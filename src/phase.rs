pub mod haplograph;
pub mod haplotree;
pub mod haplotype;
pub mod phasedblock;

use std::fs;
use std::collections::VecDeque;
use std::sync::mpsc;
use std::thread;
use std::path::{Path,PathBuf};

use anyhow::{Context, Result};
use itertools::Itertools;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use ahash::AHashMap as HashMap;
use ahash::AHashSet as HashSet;

use crate::cli::Options;
use crate::bam::BamRecordId;
use crate::seq::{SeqDatabase, SeqInterval, SuccinctSeq};
use crate::variant::{Var,VarDict};

use phasedblock::PhasedBlock;

use self::haplograph::HaploGraph;
use self::haplotype::{Haplotype, HaplotypeId, HaplotypeHit};

// input: 
//   - read-alignment (bam)
//   - work-dir
//   - target sequences
//   - variant positions
//   - other options
// output:
//   - phased haplotypes (with read files)
//   - filtered variant positions

#[derive(Default)]
pub struct PhaseResult {
    pub haplotypes: HashMap<HaplotypeId,Haplotype>,
    pub seq2haplo: HashMap<BamRecordId,Vec<HaplotypeHit>>,
}

pub struct Phaser<'a> {
    bam_path: &'a Path,
    ref_db: &'a SeqDatabase,
    read_db: &'a SeqDatabase,
    ref_intervals: &'a [SeqInterval],
    bam_index: Vec<usize>,
    fragments_dir: PathBuf,
    dots_dir: PathBuf,
    opts: &'a Options,
}

impl<'a> Phaser<'a> {

    pub fn new(bam_path: &'a Path, ref_db: &'a SeqDatabase, read_db: &'a SeqDatabase, ref_intervals: &'a[SeqInterval], work_dir: PathBuf, opts: &'a Options) -> Result<Phaser<'a>> {

        let fragments_dir = work_dir.join("fragments");
        fs::create_dir_all(&fragments_dir)
            .with_context(|| format!("Cannot create output directory: \"{}\"", fragments_dir.display()))?;

        let dots_dir = work_dir.join("dots");
        fs::create_dir_all(&dots_dir)
            .with_context(|| format!("Cannot create output directory: \"{}\"", dots_dir.display()) )?;

        let bam_index = crate::bam::build_target_index(bam_path, ref_db);

        Ok(Phaser {
            bam_path,
            ref_db,
            read_db,
            ref_intervals,
            bam_index,
            fragments_dir,
            dots_dir,
            opts
        })
    }

    pub fn fragments_dir(&self) -> &Path {
        self.fragments_dir.as_path()
    }

    pub fn dots_dir(&self) -> &Path {
        self.fragments_dir.as_path()
    }
    
    pub fn phase(&self, variants: &VarDict) -> PhaseResult {

        let mut interval_variants = vec![];
        for siv in self.ref_intervals {
            if let Some(vars) = variants.get(&siv.tid) {
                let left = vars.partition_point(|v| v.pos < siv.beg);
                let right = vars.partition_point(|v: &Var| v.pos < siv.end);
                if right - left >= self.opts.min_snv {
                    interval_variants.push((siv, &vars[left..right]));
                }
            }
        }
        interval_variants.sort_unstable_by_key(|(_,vars)| vars.len() as isize);

        let nb_threads = std::cmp::min(self.opts.nb_threads, interval_variants.len());
        let (tx, rx) = mpsc::channel();
        thread::scope(|scope| {
            for thread_id in 0..nb_threads {
                let sender = tx.clone();
                let inteval_variants_ref = &interval_variants;
                scope.spawn(move || {
                    for &(target_interval, variants) in inteval_variants_ref[thread_id..].iter().step_by(nb_threads) {
                        let result = self.phase_interval(target_interval, variants);
                        sender.send(result).unwrap();
                    }
                });
            }
        });

        let mut result = PhaseResult::default();
        for _ in 0..interval_variants.len() {
            let PhaseResult { haplotypes, seq2haplo } = rx.recv().unwrap();
            result.haplotypes.extend(haplotypes);
            result.seq2haplo.extend(seq2haplo);
        }
        result
    }

    fn phase_interval(&self, ref_interval:&SeqInterval, variants: &[Var]) -> PhaseResult {

        // spdlog::debug!("----- Phasing {target_interval} ({} variants) -----", variants.len());
        let (haplotypes, succinct_seqs) = self.phase_variants(ref_interval,variants);

        // for interval in haplotypes.keys().sorted_unstable() {
        //     if let Some(haplotypes) = haplotypes.get(interval) {
        //         if haplotypes.len() != 3 && haplotypes.first().is_some_and(|ht| ht.raw_size() > 1) {
        //             spdlog::trace!("Local Haplotypes @ {}:{}-{}", self.ref_db.names[interval.tid].as_str(), interval.beg, interval.end);
        //             haplotypes.iter().for_each(|ht| spdlog::trace!("  * {ht}"));
        //         }
        //     }
        // }
        
        let haplotypes: HashMap<HaplotypeId, Haplotype> = haplotypes.into_values().flatten().map(|ht| (ht.uid(),ht)).collect();
        let variant_positions = self.variant_positions(&haplotypes);
        let haplotypes: HashMap<HaplotypeId,Haplotype> = self.remove_false_haplotypes(haplotypes, &variant_positions);

        let sread_haplotypes = self::separate_reads(&succinct_seqs, &haplotypes, 1);

        let mut haplograph = HaploGraph::new(haplotypes, sread_haplotypes);
        let dot_file = format!("{}_{}-{}.raw.dot", self.ref_db.names[ref_interval.tid], ref_interval.beg, ref_interval.end);
        haplograph.write_dot(&self.dots_dir.join(dot_file)).unwrap();
        
        // Merge contiguous haplotypes when it is not ambiguous to do so
        let haplotypes: HashMap<HaplotypeId,Haplotype> = haplograph.scaffold_haplotypes(&variant_positions, 5, 0.8, self.opts.min_snv);
        let dot_file = format!("{}_{}-{}.final.dot", self.ref_db.names[ref_interval.tid], ref_interval.beg, ref_interval.end);
        haplograph.write_dot(&self.dots_dir.join(dot_file)).unwrap();

        let seq2haplo = self::separate_reads(&succinct_seqs, &haplotypes, 1);
        self.write_reads(&haplotypes, &seq2haplo).unwrap();

        // if self.opts.trace {
        //     let mut haplo_bases: HashMap<HaplotypeId, usize> = HashMap::new();
        //     for hits in seq2haplo.values() {
        //         for h in hits.iter().filter(|h| !h.is_ambiguous()) {
        //             haplo_bases.entry(h.hid).and_modify(|e| *e += h.nb_pos).or_insert(h.nb_pos);
        //         }
        //     }
        //     for (hid, bases) in haplo_bases {
        //         if bases / haplotypes[&hid].raw_size() < self.opts.min_alt_count {
        //             spdlog::warn!("Haplotype {hid} has low coverage: {}", bases / haplotypes[&hid].raw_size());
        //         }
        //     }
        // }

        PhaseResult {
            haplotypes,
            seq2haplo,
        }
    }


    fn write_reads(&self, haplotypes: &HashMap<HaplotypeId,Haplotype>, sread_haplotypes: &HashMap<BamRecordId,Vec<HaplotypeHit>>) -> std::io::Result<()> {
        
        let mut hap_to_reads: HashMap<HaplotypeId, Vec<usize>> = HashMap::new();
        for (record_id,hits) in sread_haplotypes.iter() {
            for hid in hits.iter().filter(|hit| hit.nb_alt == 0 && hit.dist != hit.nb_pos).map(|hit| hit.hid) {
                hap_to_reads.entry(hid).or_default().push(record_id.index);
            }
        }

        for ht_id in haplotypes.keys() {
            let ht_file_gz = format!("{}_{}-{}_h{}.fa.gz", self.ref_db.names[ht_id.tid], ht_id.beg, ht_id.end, ht_id.hid);
            let ht_file_path = self.fragments_dir().join(ht_file_gz.as_str());
            let mut writer = crate::utils::get_file_writer(&ht_file_path);
            if let Some(read_indices) = hap_to_reads.get(ht_id) {
                for read_idx in read_indices {
                    writer.write_all(format!(">{read_idx}\n").as_bytes())?;
                    writer.write_all(&self.read_db.sequences[*read_idx].as_bytes())?;
                    writer.write_all(b"\n")?;
                }
            }
        }

        // let mut read_files: HashMap<HaplotypeId,_> = HashMap::new();
        
        // for ht_id in haplotypes.keys() {
        //     let ht_file_gz = format!("{}_{}-{}_h{}.fa.gz", self.ref_db.names[ht_id.tid], ht_id.beg, ht_id.end, ht_id.hid);
        //     let ht_file_path = self.fragments_dir().join(ht_file_gz.as_str());
        //     read_files.insert(*ht_id, crate::utils::get_file_writer(&ht_file_path));
        // }
        
        // for (record_id, hits) in sread_haplotypes.iter() {
        //     for hit in hits.iter().filter(|hit| hit.nb_alt == 0) {
        //         let out = read_files.get_mut(&hit.hid).unwrap();
        //         out.write_all(format!(">{}\n", &record_id.index).as_bytes())?;
        //         out.write_all(&self.read_db.sequences[record_id.index].as_bytes())?;
        //         out.write_all(b"\n")?;
        //     }
        // }

        Ok(())
    }

    fn remove_false_haplotypes(&self, mut haplotypes: HashMap<HaplotypeId,Haplotype>, positions: &HashSet<(usize,usize)>) -> HashMap<HaplotypeId,Haplotype> {
        haplotypes.retain(|_, ht| {
            ht.variants()
                .iter()
                .map(|snv| (ht.tid(),snv.pos))
                .any(|pos| positions.contains(&pos))
        });
        haplotypes
    }

    fn variant_positions(&self, haplotypes: &HashMap<HaplotypeId,Haplotype>) -> HashSet<(usize,usize)> {

        let mut counter: HashMap<_,usize> = HashMap::default();
        for (tid,snv_pos) in haplotypes.values().flat_map(|ht| ht.variants().iter().map(|snv| (ht.tid(),snv.pos))) {
            counter.entry((tid,snv_pos))
                .and_modify(|cnt| *cnt+=1)
                .or_insert(1);
        }

        counter.into_iter()
            .filter(|&(_,cnt)| cnt > 1)
            .map(|(pos,_)| pos)
            .collect()
    }

    pub fn phase_variants(&self, ref_interval:&SeqInterval, variants: &[Var]) -> (HashMap<SeqInterval,Vec<Haplotype>>,Vec<SuccinctSeq>)  {

        if variants.is_empty() {
            return (HashMap::default(),vec![]);
        }

        // spdlog::trace!("phasing variants @ {ref_interval}");

        let mut haplotypes: HashMap<SeqInterval,Vec<Haplotype>> = HashMap::new();
        let mut succinct_records: HashMap<BamRecordId,SuccinctSeq> = HashMap::new();

        let mut supporting_sreads = HashSet::new();
        let mut lookback_positions: VecDeque<usize> = VecDeque::new();
        let mut var_iter = variants.iter();
        let mut phasedblock = PhasedBlock::new(ref_interval.tid);

        let mut bam_reader = bam::IndexedReader::from_path(self.bam_path).unwrap();
        let bam_region = (self.bam_index[ref_interval.tid] as i32, ref_interval.beg as i64, ref_interval.end as i64);
        bam_reader.fetch(bam_region).unwrap_or_else(|_| panic!("Cannot fetch bam interval: {:?}", bam_region));

        let mut var = var_iter.next();
        for pileup in bam_reader.pileup().flatten() {
            // possibly update snv_position
            if var.is_none() {
                break 
            };

            let target_pos = pileup.pos() as usize;
            let mut var_position = var.unwrap().pos;
            let mut var_nucleotides = var.unwrap().alleles.iter().map(|x| x.0 as u8).collect_vec();
            if var_position < target_pos {
                var = var_iter.next();
                if var.is_none() { break; }
                var_position = var.unwrap().pos;
                var_nucleotides = var.unwrap().alleles.iter().map(|x| x.0 as u8).collect_vec();
            }
            if target_pos < var_position {
                continue;
            }
            debug_assert_eq!(target_pos, var_position);

            // load edges from pileup and previous variant position
            let edges = self.process_pileup(&pileup, ref_interval.tid, &mut succinct_records, &mut supporting_sreads, &var_nucleotides);
            lookback_positions.push_back(var_position);

            // if this is the start of a new phased block
            if lookback_positions.len() == 1 {
                phasedblock.init(var_position, var_nucleotides);
                continue
            }

            // if no available edges
            if edges.is_empty() {
                if var_position - phasedblock.begin() > self.opts.lookback {
                    let interval = phasedblock.interval().unwrap();
                    haplotypes.insert(interval, phasedblock.drain());
                }
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue
            }

            // identify haplotypes that cannot be extended with current edges
            let unsupported_haplotypes = phasedblock.haplotypes().iter()
                .filter(|&(_, ht)| edges.iter().all(|&(s,_)| s != ht.last_nuc() ))
                .map(|(&hid,_)| hid)
                .collect_vec();

            // possibly save haplotypes and start a new phased block
            if !unsupported_haplotypes.is_empty() && (var_position - phasedblock.begin() > self.opts.lookback) {
                // for hid in unsupported_haplotypes {
                //     spdlog::trace!("Cannot edge-extend for haplotypes @ {}:{target_pos}", self.ref_db.names[ref_interval.tid].as_str());
                //     spdlog::trace!("  * {}", phasedblock.haplotypes().get(&hid).unwrap());
                // }
                let phased_interval = phasedblock.interval().unwrap();
                haplotypes.insert(phased_interval, phasedblock.drain());
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue
            }
            for hid in unsupported_haplotypes {
                phasedblock.remove_haplotype(hid)
            }

            if phasedblock.haplotypes().len() < 2 {
                phasedblock.init(var_position, var_nucleotides);
                continue
            }

            let mut _is_ambiguous = phasedblock.extend(var_position, edges);

            // discard haplotypes unsupported by the reads
            let back_pos = (var_position+1).saturating_sub(self.opts.lookback.min(1000));
            let back_i = lookback_positions.partition_point(|&pos| pos < back_pos);
            let min_position = lookback_positions[back_i];
            supporting_sreads.retain(|sr_id| succinct_records[sr_id].positions().last().unwrap() == &var_position);
            let candidate_records = supporting_sreads.iter()
                .filter(|&sr_id| succinct_records[sr_id].positions()[0] <= min_position)
                .collect_vec();

            // spdlog::trace!("Current haplotypes @ {}:{target_pos}", self.ref_db.names[ref_interval.tid].as_str());
            // for hid in phasedblock.haplotypes().keys() {
            //     spdlog::trace!("  * {}", phasedblock.haplotypes().get(hid).unwrap());
            // }

            let (unsupported_haplotypes, ambiguous_haplotypes) = self.validate_haplotypes(&succinct_records, &candidate_records, &phasedblock);

            if !unsupported_haplotypes.is_empty() && (var_position - phasedblock.begin() + 1 > self.opts.lookback) {
                phasedblock.split_and_init(0);
                let phased_interval = phasedblock.interval().unwrap();
                haplotypes.insert(phased_interval, phasedblock.drain());
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue;
            }
            if !unsupported_haplotypes.is_empty() {
                for hid in unsupported_haplotypes {
                    phasedblock.remove_haplotype(hid);
                }
            }

            let is_ambiguous = !ambiguous_haplotypes.is_empty(); // is_ambiguous |= !ambiguous_haplotypes.is_empty();

            // if ambiguous extension, create a new phaseset
            if is_ambiguous && (var_position - phasedblock.begin() > self.opts.lookback) {
                phasedblock.split_and_init(0);
                let phased_interval = phasedblock.interval().unwrap();
                haplotypes.insert(phased_interval, phasedblock.drain());
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue;
                // let mut new_phasedblock = phasedblock.split_and_init(self.opts.lookback);
                // let phased_interval = phasedblock.interval().unwrap();
                // haplotypes.insert(phased_interval, phasedblock.drain());
                // std::mem::swap(&mut phasedblock, &mut new_phasedblock);
                // continue;
            }
            // discard ambiguous haplotypes if phased region was too short
            for hid in ambiguous_haplotypes {
                phasedblock.remove_haplotype(hid);
            }

            if phasedblock.is_empty() {
                phasedblock.init(var_position, var_nucleotides);
                lookback_positions.clear();
                lookback_positions.push_back(var_position);
                continue;
            }
        }

        let phased_interval = phasedblock.interval().unwrap();
        haplotypes.insert(phased_interval, phasedblock.drain());
        
        (haplotypes, succinct_records.into_values().collect_vec())
    }

    fn validate_haplotypes(&self, succinct_records: &HashMap<BamRecordId,SuccinctSeq>, supporting_records: &[&BamRecordId], phasedblock: &PhasedBlock) -> (Vec<usize>,Vec<usize>) {

        fn get_best_haplotypes(sr: &SuccinctSeq, haplotypes:&[&Haplotype]) -> Vec<(usize,usize)> {
            let mut sr_distances = Vec::with_capacity(haplotypes.len());
            for ht in haplotypes {
                let back_i = ht.raw_variants().partition_point(|snv| snv.pos < sr.positions()[0]);
                let back_pos = ht.raw_variants().get(back_i).unwrap().pos;
                let sr_left = sr.positions().partition_point(|pos| *pos < back_pos);
                let sr_right = sr.positions().partition_point(|pos| *pos <= ht.last_pos());
                let sr_nucleotides = sr.nucleotides().get(sr_left..sr_right).unwrap();
                let sr_qualities = sr.qualities().get(sr_left..sr_right).unwrap();
                let ht_nucleotides = ht.raw_variants().get(back_i..).unwrap();

                let nb_pos = sr_nucleotides.len();
                let dist = itertools::izip!(ht_nucleotides, sr_nucleotides, sr_qualities)
                    .filter(|&(a, b, qual)| *qual < 10 || a.nuc != *b)
                    .count();

                assert!(nb_pos > 0);
                sr_distances.push((ht.hid(), dist, nb_pos));
            }
            if sr_distances.is_empty() {
                return Vec::new()
            }
            let min_dist = sr_distances.iter().map(|(_,d,_)| *d).min().unwrap();
            sr_distances.into_iter()
                .filter(|(_,d,_)| *d == min_dist)
                .map(|(hid,_,nb_pos)| (hid,nb_pos))
                .collect_vec()
        }
        
        let mut supporting: HashMap<usize,usize> = HashMap::new();
        let mut unambiguous: HashMap<usize,usize> = HashMap::new();
        let haplotypes = phasedblock.haplotypes().values().collect_vec();
        let nb_snvs = phasedblock.haplotypes().values().next().unwrap().raw_size();

        for sr_id in supporting_records {
            let sr = &succinct_records[sr_id];
            let mut best_haplotypes = get_best_haplotypes(sr, &haplotypes);
            if best_haplotypes.len() == 1 {
                let (hid,nb_pos) = best_haplotypes.pop().unwrap();
                let nb_pos = nb_pos.min(nb_snvs);
                unambiguous.entry(hid).and_modify(|cnt| *cnt+=nb_pos).or_insert(nb_pos);
                supporting.entry(hid).and_modify(|cnt| *cnt+=nb_pos).or_insert(nb_pos);
            }
            while let Some((hid,nb_pos)) = best_haplotypes.pop() {
                supporting.entry(hid).and_modify(|cnt| *cnt+=nb_pos).or_insert(nb_pos);
            }
        }

        // if self.ref_db.names[phasedblock.interval().unwrap().tid].as_str() == "edge_482" && [269889].contains(&phasedblock.begin()) {
        //     spdlog::trace!("Supporting:");
        //     for (hid,count) in supporting.iter() {
        //         spdlog::trace!("  * h{hid} => {:.1}", *count as f64 / nb_snvs as f64);
        //     }
        //     spdlog::trace!("Unambiguous:");
        //     for (hid,count) in unambiguous.iter() {
        //         spdlog::trace!("  * h{hid} => {:.1}", *count as f64 / nb_snvs as f64);
        //     }
        // }

        let unsupported_haplotypes = phasedblock.haplotypes().keys()
            .filter(|hid| supporting.get(hid).is_none_or(|cnt| *cnt/nb_snvs < self.opts.min_alt_count))
            .cloned()
            .collect_vec();

        let ambiguous_haplotypes = phasedblock.haplotypes().keys()
            .filter(|hid| supporting.get(hid).is_some_and(|cnt| *cnt/nb_snvs >= self.opts.min_alt_count) && unambiguous.get(hid).is_none_or(|cnt| *cnt/nb_snvs < self.opts.min_alt_count))
            .cloned()
            .collect_vec();

        (unsupported_haplotypes, ambiguous_haplotypes)
    }

    fn process_pileup(&self, pileup: &bam::pileup::Pileup, ref_idx: usize, succinct_records: &mut HashMap<BamRecordId,SuccinctSeq>, supporting_sreads: &mut HashSet<BamRecordId>, var_nucleotides: &[u8]) -> Vec<(u8,u8)> {

        let position = pileup.pos() as usize;

        // let mut edge_total_obs: usize = 0;
        let mut edge_counter: HashMap<(u8,u8),usize> = HashMap::new();

        for alignment in pileup.alignments() {
            let record = alignment.record();

            if record.mapq() < self.opts.min_mapq
                || record.is_unmapped() 
                || record.is_secondary() 
                || record.is_quality_check_failed() 
                || record.is_duplicate()
            {
                continue;
            }
            
            let record_id = BamRecordId::from_record(&record, self.read_db);
            
            let record_nuc = if !alignment.is_del() && !alignment.is_refskip() { 
                alignment.record().seq()[alignment.qpos().unwrap()]
            } else { 
                b'-'
            };
            
            let record_qual = if !alignment.is_del() && !alignment.is_refskip() { 
                alignment.record().qual()[alignment.qpos().unwrap()]
            } else {
                0
            };

            if !succinct_records.contains_key(&record_id) {
                let srec = SuccinctSeq::new(record_id, ref_idx);
                succinct_records.insert(record_id, srec);
                supporting_sreads.insert(record_id);
            }

            let srec = succinct_records.get_mut(&record_id).unwrap();
            srec.push(position, record_nuc, record_qual);

            let srec_positions = srec.positions();
            let srec_qualities = srec.qualities();
            let srec_len = srec.len();
            if srec_len > 1 && (srec_positions[srec_len-1] - srec_positions[srec_len-2] < self.opts.lookback) && srec_qualities[srec_len-1] >= 10 && srec_qualities[srec_len-2] >= 10 {
                let srec_nucleotides = srec.nucleotides();
                let edge = (srec_nucleotides[srec_len-2], srec_nucleotides[srec_len-1]);
                if var_nucleotides.contains(&edge.1) {
                    // edge_total_obs += 1;
                    edge_counter.entry(edge)
                        .and_modify(|cnt| *cnt += 1)
                        .or_insert(1);
                }
            }
        }

        edge_counter.iter()
            .filter(|&(_edge,&cnt)| cnt >= self.opts.min_alt_count) // && (cnt as f64) >= self.opts.min_alt_frac * (edge_total_obs as f64))
            .map(|(&edge,_cnt)| edge)
            .collect_vec()
    }

}


pub fn separate_reads(succinct_records: &[SuccinctSeq], haplotypes: &HashMap<HaplotypeId,Haplotype>, min_shared_pos: usize)
    ->  HashMap<BamRecordId,Vec<HaplotypeHit>> {

    let haplotypes = haplotypes.values().sorted_unstable_by_key(|ht| ht.uid()).collect_vec();

    let mut sread_haplotypes: HashMap<BamRecordId, Vec<HaplotypeHit>> = HashMap::new();
    for sread in succinct_records {
        let best_hits = best_sread_haplotypes(sread, &haplotypes, min_shared_pos);
        sread_haplotypes.entry(sread.record_id()).or_default()
            .extend(best_hits);
    }
    sread_haplotypes
}

// Prerequisite: haplotypes are sorted by target-id and start/end coordinate
fn best_sread_haplotypes(sread: &SuccinctSeq, haplotypes: &[&Haplotype], min_shared_pos: usize)
    -> Vec<HaplotypeHit> {
    
    let mut candidates = HashMap::new();

    let idx = haplotypes.partition_point(|ht| ht.tid() < sread.tid() || (ht.tid() == sread.tid() && ht.end() <= sread.beg()));
    for ht in haplotypes[idx..].iter().take_while(|ht| ht.beg() < sread.end()) {
        if let Some(hit) = sread_haplotype_distance(sread, ht) {
            if hit.nb_pos > 0 && hit.nb_pos >= min_shared_pos && hit.dist < hit.nb_pos {
                let range = ht.beg()..ht.end();
                candidates.entry(range)
                    .or_insert(vec![])
                    .push(hit);
            }
        }
    }

    let mut best_hits = vec![];
    for hits in candidates.into_values() {
        assert!(!hits.is_empty());

        let mut hits_iter = hits.into_iter();
        let mut best_hit = hits_iter.next().unwrap();
        for hit in hits_iter {
            match hit.dist.cmp(&best_hit.dist) {
                std::cmp::Ordering::Equal => { best_hit.nb_alt += 1; },
                std::cmp::Ordering::Less  => { best_hit = hit; },
                _ => {}
            }
        }
        best_hits.push(best_hit);
    }

    best_hits
}


// returns hamming distance and number of shared positions (None if no shared position)
fn sread_haplotype_distance(sread: &SuccinctSeq, ht: &Haplotype) -> Option<HaplotypeHit> {

    let sread_positions: HashSet<usize> = sread.positions().iter().cloned().collect();
    let ht_variants = ht.raw_variants().as_slice();

    let ht_beg = ht_variants.partition_point(|snv| snv.pos < *sread.positions().first().unwrap());
    let ht_end = ht_variants.partition_point(|snv| snv.pos <= *sread.positions().last().unwrap());
    if ht_beg >= ht_end {
        return None
    }

    let ht_positions: HashSet<usize> = ht_variants[ht_beg..ht_end].iter().map(|snv| snv.pos).collect();
    let shared_positions: HashSet<usize> = ht_positions.intersection(&sread_positions).cloned().collect();

    let sread_nucleotides = sread.positions().iter().enumerate()
        .filter_map(|(i,pos)| if shared_positions.contains(pos) { Some(sread.nucleotides()[i]) } else { None })
        .collect_vec();

    let sread_qualities = sread.positions().iter().enumerate()
        .filter_map(|(i,pos)| if shared_positions.contains(pos) { Some(sread.qualities()[i]) } else { None })
        .collect_vec();

    let ht_nucleotides = ht_variants.iter()
        .filter_map(|&snv| if shared_positions.contains(&snv.pos) { Some(snv.nuc) } else { None })
        .collect_vec();

    let nb_pos = sread_nucleotides.len();
    let dist = itertools::izip!(ht_nucleotides, sread_nucleotides, sread_qualities)
        .filter(|&(a, b, qual)| qual < 10 || a != b)
        .count();

    let ht_hit = HaplotypeHit::new(ht.uid(), dist, nb_pos, 0);
    Some(ht_hit)
}
