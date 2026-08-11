use std::fmt;

use itertools::Itertools;
use ahash::AHashMap as HashMap;

use crate::phase::haplotype::{HaplotypeId, HaplotypeHit, Haplotype};
use crate::seq::SeqInterval;
use crate::alignment::{MappingType, SeqAlignment};
use crate::bam::BamRecordId;


#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Haplotype(usize),
    #[default] Unphased,
    Read, // read strand
}

impl fmt::Display for SeqType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeqType::Haplotype(hid) => { write!(f, "Haplotype({hid})") }
            SeqType::Unphased => { write!(f, "Unphased") }
            SeqType::Read => { write!(f, "Read") }
        }
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct AwareContig {
    contig_type: SeqType,
    interval: SeqInterval,
    strand: u8,
    aligned_bases: usize,
}

impl AwareContig {

    pub fn new(contig_type:SeqType, interval:SeqInterval, strand:u8, aligned_bases:usize) -> Self {
        Self {
            contig_type,
            interval,
            strand,
            aligned_bases
        }
    }

    pub fn tid(&self) -> usize { self.interval.tid }
    pub fn beg(&self) -> usize { self.interval.beg }
    pub fn end(&self) -> usize { self.interval.end }
    pub fn length(&self) -> usize { self.interval.length() }
    pub fn strand(&self) -> u8 { self.strand }

    pub fn contig_type(&self) -> SeqType { self.contig_type }
    pub fn is_phased(&self) -> bool { matches!(self.contig_type, SeqType::Haplotype(_)) }
    
    pub fn hid(&self) -> Option<usize> {
        match self.contig_type {
            SeqType::Haplotype(hid) => Some(hid),
            _ => None
        }
    }

    pub fn haplotype_id(&self) -> Option<HaplotypeId> {
        Some(HaplotypeId{
            tid: self.interval.tid,
            beg: self.interval.beg,
            end: self.interval.end,
            hid: self.hid()?
        })
    }

    pub fn interval(&self) -> SeqInterval { self.interval }

    // pub fn phaseset_id(&self) -> SeqInterval { self.interval() }
    // pub fn region(&self) -> (usize,usize) { (self.beg(), self.end()) }

    pub fn depth(&self) -> f64 {
        let contig_length = self.length();
        if contig_length > 0 {
            (self.aligned_bases as f64) / (contig_length as f64)
        } else {
            0.0
        }
    }

}

#[derive(Debug, Clone)]
pub struct AwareAlignment {
    pub aware_id: usize,
    pub query_idx: usize,
    pub query_len: usize,
    pub query_beg: usize,
    pub query_end: usize,
    pub strand: u8,
    pub target_idx: usize,
    pub target_beg: usize,
    pub target_end: usize,
    pub query_aln_beg: usize,
    pub query_aln_end: usize,
    pub mapping_type: MappingType,
    pub dist: usize,
    pub nb_shared_snvs: usize,
    pub nb_alt_matches: usize,
}

impl AwareAlignment {

    pub fn is_ambiguous(&self) -> bool {
        self.nb_alt_matches > 0
    }

    pub fn target_interval(&self) -> SeqInterval {
        SeqInterval { tid: self.target_idx, beg: self.target_beg, end: self.target_end }
    }

}


fn retrieve_unphased_intervals(ref_intervals:&[SeqInterval], haplotypes:&HashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<SeqInterval> {

    let mut unphased_intervals: Vec<SeqInterval> = vec![];
    let phased_intervals = haplotypes.values()
        .map(|ht| ht.region())
        .unique()
        .sorted_unstable()
        .collect_vec();

    for iv in ref_intervals {

        let first = phased_intervals.partition_point(|x| x.tid < iv.tid || (x.tid == iv.tid && x.beg < iv.beg));
        
        let mut beg = iv.beg;
        for &piv in phased_intervals[first..].iter().take_while(|piv| piv.tid == iv.tid && piv.end <= iv.end) {
            let uiv = SeqInterval{ tid:iv.tid, beg, end:piv.beg };
            if !uiv.empty() && uiv.length() >= min_length {
                unphased_intervals.push(uiv);
            }
            beg = piv.end;
        }

        let uiv = SeqInterval{ tid:iv.tid, beg, end:iv.end };
        if !uiv.empty() && uiv.length() >= min_length {
            unphased_intervals.push(uiv);
        }
    }

    unphased_intervals
}


pub fn build_aware_contigs(ref_intervals:&[SeqInterval], haplotypes:&HashMap<HaplotypeId,Haplotype>, min_length:usize) -> Vec<AwareContig> {
    
    let mut aware_contigs = retrieve_unphased_intervals(ref_intervals, haplotypes, min_length)
        .into_iter()
        .map(|interval| AwareContig{
            contig_type: SeqType::Unphased,
            interval,
            strand: b'+',
            aligned_bases: 0
        }).collect_vec();

    aware_contigs.extend(haplotypes.values()
        .map(|ht| {
            AwareContig {
                contig_type: SeqType::Haplotype(ht.hid()),
                interval: ht.region(),
                strand: b'+',
                aligned_bases: 0
            }
        })
    );

    aware_contigs.sort_by_key(|ctg| (ctg.interval(), ctg.hid()));
    spdlog::debug!("{} strain-aware contigs built", aware_contigs.len());

    aware_contigs
}


pub fn build_phased_contigs(haplotypes:&HashMap<HaplotypeId,Haplotype>) -> Vec<AwareContig> {
    
    let mut aware_contigs = haplotypes.values()
        .map(|ht| {
            AwareContig {
                contig_type: SeqType::Haplotype(ht.hid()),
                interval: ht.region(),
                strand: b'+',
                aligned_bases: 0
            }
        }).collect_vec();

    aware_contigs.sort_by_key(|ctg| (ctg.interval(), ctg.hid()));
    spdlog::debug!("{} strain-aware contigs built", aware_contigs.len());

    aware_contigs
}


pub fn map_sequences_to_aware_contigs(seq_alignments: &HashMap<usize,Vec<SeqAlignment>>, aware_contigs: &mut[AwareContig], seq2haplo: &HashMap<BamRecordId,Vec<HaplotypeHit>>) -> HashMap<usize,Vec<AwareAlignment>> {

    let mut read2aware = HashMap::new();
    for (read_idx, alignments) in seq_alignments {
        let mut aware_alignments = map_alignments_to_aware_contigs(alignments, aware_contigs, seq2haplo);
        aware_alignments.sort_unstable_by_key(|a| (a.query_aln_beg,a.query_beg));
        let aware_alignments = merge_aware_alignments(aware_alignments, aware_contigs);
        let aware_alignments = filter_aware_alignments(aware_alignments, aware_contigs);
        read2aware.insert(*read_idx, aware_alignments);
    }

    read2aware
}


// aware_contigs should be sorted by interval
pub fn map_alignments_to_aware_contigs(alignments: &[SeqAlignment], aware_contigs: &[AwareContig], seq2haplo: &HashMap<BamRecordId,Vec<HaplotypeHit>>) -> Vec<AwareAlignment> {
    
    let mut aware_alignments: Vec<AwareAlignment> = vec![];

    for sa in alignments {
        let sa_id = sa.bam_record_id();

        let idx = aware_contigs.partition_point(|ctg| ctg.tid() < sa.tid() || (ctg.tid() == sa.tid() && ctg.end() <= sa.target_beg()));
        for (i,ctg) in aware_contigs[idx..].iter().take_while(|ctg| ctg.tid() == sa.tid() && ctg.beg() < sa.target_end()).enumerate() {

            if ctg.is_phased() && (!seq2haplo.contains_key(&sa_id) || !seq2haplo[&sa_id].iter().any(|hit| hit.hid == ctg.haplotype_id().unwrap())) {
                continue
            }

            let aware_id = idx+i;
            
            let aware_target_beg = std::cmp::max(sa.target_beg(), ctg.beg());
            let aware_target_end = std::cmp::min(sa.target_end(), ctg.end());
            
            let aware_ctg_beg = aware_target_beg - ctg.beg();
            let aware_ctg_end = aware_target_end - ctg.beg();
            let aware_ctg_range = (aware_ctg_beg, aware_ctg_end, ctg.length());

            let query_beg = if sa.is_forward() { sa.query_beg() } else { sa.query_length() - sa.query_end() };
            let (mut aware_query_beg, aware_target_beg) = crate::alignment::first_match_from(ctg.beg(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let (mut aware_query_end, aware_target_end) = crate::alignment::last_match_until(ctg.end(), query_beg, sa.target_beg(), sa.cigar()).unwrap();
            let aware_query_range = (aware_query_beg, aware_query_end, sa.query_length());

            if aware_query_end <= aware_query_beg || aware_target_end <= aware_target_beg {
                continue
            }

            let maptype = crate::alignment::classify_mapping(aware_query_range, aware_ctg_range, 100, 0.5);

            if sa.is_reverse() {
                (aware_query_beg, aware_query_end) = (sa.query_length() - aware_query_end, sa.query_length() - aware_query_beg);
            }

            let mut nb_alt_matches = 0;
            let mut nb_shared_snvs = 0;
            let mut dist = 0;

            if ctg.is_phased() {
                debug_assert!(seq2haplo.contains_key(&sa_id));
                if let Some(hit) = seq2haplo[&sa_id].iter().find(|hit| hit.hid == ctg.haplotype_id().unwrap()) {
                    nb_alt_matches = hit.nb_alt;
                    nb_shared_snvs = hit.nb_pos;
                    dist = hit.dist;
                }
            }

            let awa = AwareAlignment{
                aware_id,
                query_idx: sa.query_index(),
                query_len: sa.query_length(),
                query_beg: aware_query_beg,
                query_end: aware_query_end,
                strand: sa.strand(),
                target_idx: sa.target_index(),
                target_beg: aware_target_beg,
                target_end: aware_target_end,
                query_aln_beg: sa.query_beg(),
                query_aln_end: sa.query_end(),
                mapping_type: maptype,
                dist,
                nb_shared_snvs,
                nb_alt_matches,
            };

            aware_alignments.push(awa);
        }
    }

    aware_alignments
}


pub fn merge_aware_alignments(aware_alignments: Vec<AwareAlignment>, aware_contigs: &[AwareContig]) -> Vec<AwareAlignment> {

    let mut merged_alignments: Vec<AwareAlignment> = Vec::with_capacity(aware_alignments.len());
    for a in aware_alignments {
        if let Some(last) = merged_alignments.last_mut()
            && a.aware_id == last.aware_id && a.strand == last.strand && ((a.strand == b'+' && last.target_beg < a.target_beg) || (a.strand != b'+' && a.target_beg < last.target_beg)) {
                let ctg = aware_contigs[a.aware_id];
                last.query_end = last.query_end.max(a.query_end);
                last.target_beg = last.target_beg.min(a.target_beg);
                last.target_end = last.target_end.max(a.target_end);
                let query_range = if a.strand == b'+' { (a.query_beg, a.query_end, a.query_len) } else { (a.query_len-a.query_end, a.query_len-a.query_beg, a.query_len) };
                let target_range = (last.target_beg-ctg.beg(), last.target_end-ctg.beg(), ctg.length());
                // spdlog::trace!("ctg={}  qry={query_range:?} trg={target_range:?}", ctg.interval());
                last.mapping_type = crate::alignment::classify_mapping(query_range, target_range, 100, 0.5);
                last.dist += a.dist;
                last.nb_shared_snvs += a.nb_shared_snvs;
                last.nb_alt_matches = last.nb_alt_matches.max(a.nb_alt_matches);
                continue
        }
        merged_alignments.push(a);
    }

    merged_alignments
}


pub fn filter_aware_alignments(aware_alignments: Vec<AwareAlignment>, aware_contigs: &mut [AwareContig]) -> Vec<AwareAlignment> {
    let nb_alignments = aware_alignments.len();
    aware_alignments.into_iter()
        .enumerate()
        .filter(|(_,awa)| awa.mapping_type.is_containment() || awa.mapping_type.is_dovetail())
        .filter(|(_,awa)| awa.query_end - awa.query_beg >= 1)
        .filter(|(i,awa)| (1..nb_alignments).contains(i) || 10 * awa.dist <= awa.nb_shared_snvs)
        .map(|(_,awa)| {
            let ctg = &mut aware_contigs[awa.aware_id];
            ctg.aligned_bases += awa.target_end - awa.target_beg;
            awa
        }).collect()
}