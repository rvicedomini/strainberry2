pub mod alignment;
pub mod polisher;
pub mod window;

use std::{path::Path, str::FromStr};

use ahash::AHashMap as HashMap;
use anyhow::{bail,Context, Result};
use itertools::Itertools;

use crate::bitseq::BitSeq;
use crate::polish::alignment::Alignment;

const ERROR_THRESHOLD: f64 = 0.3;


pub fn compute_alignments(target_path: &Path, read_path: &Path, opts: &crate::cli::Options) -> Result<Vec<Alignment>> {

    use std::io::{BufRead,BufReader};
    use std::process::{Command, Stdio};

    let mm2_preset = match opts.mode {
        crate::cli::Mode::Hifi => "-xmap-hifi",
        crate::cli::Mode::Nano => "-xlr:hq",
    };
    
    let mm2_args = ["-t", &opts.nb_threads.to_string(), "-c", mm2_preset, target_path.to_str().unwrap(), read_path.to_str().unwrap()];
    let mut minimap2 = Command::new("minimap2")
        .args(mm2_args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .context("cannot run minimap2")?;

    let mm2_stdout = minimap2.stdout.take()
        .context("Could not capture minimap2 standard output.")?;

    let reader = BufReader::new(mm2_stdout);
    let alignments = reader.lines()
        .map_while(Result::ok)
        .filter_map(|line| {
            let line = line.trim();
            if !line.is_empty() { Some(line.to_string()) } else { None }
        })
        .map(|line| {
            Alignment::from_str(&line)
                .with_context(||format!("error parsing paf record:\n{line}"))
                .unwrap()
        })
        .collect_vec();

    let mm2_status = minimap2.wait().context("minimap2 did not run successfully")?;
    if !mm2_status.success() {
        bail!("minimap2 exited with {mm2_status}");
    }

    // filter alignments

    let mut retained: HashMap<usize, Alignment> = HashMap::new();

    for a in alignments {
        let a_err = 1.0 - std::cmp::min(a.query_end-a.query_beg, a.target_end-a.target_beg) as f64 / a.length  as f64;
        if a_err > ERROR_THRESHOLD {
            continue;
        }
        
        retained.entry(a.query_idx)
            .and_modify(|e| {
                if a.length > e.length { *e = a.clone(); }
            }).or_insert(a);
    }

    let alignments = retained.into_values().collect();

    Ok(alignments)
}


pub fn polish(ref_sequences: &[BitSeq], read_sequences: &[BitSeq], mut alignments: Vec<Alignment>) -> Vec<BitSeq> {

    let mut polisher = polisher::Polisher::new(ref_sequences, read_sequences, &mut alignments);

    // spdlog::debug!("initializing polisher");
    polisher.initialize();

    // spdlog::debug!("polishing windows");
    let polished_sequences = polisher.polish();

    polished_sequences.into_iter()
        .map(|seq| BitSeq::from_utf8(&seq))
        .collect_vec()
}


// pub fn build_haplotigs(ref_db: &SeqDatabase, read_db: &SeqDatabase, haplotypes: &HashMap<HaplotypeId,Haplotype>, aware_contigs: &[AwareContig], read2aware: &HashMap<usize,Vec<AwareAlignment>>) -> HashMap<HaplotypeId,Vec<u8>> {

//     let haplotypes = haplotypes.keys().cloned().sorted().collect_vec();
//     let mut hap_index = HashMap::new();
//     let hap_sequences = haplotypes.iter().enumerate().map(|(idx, hid)| {
//         hap_index.insert(*hid, idx);
//         ref_db.sequences[hid.tid][hid.beg..hid.end].to_vec()
//     }).collect_vec();

//     let mut alignments = Vec::new();
//     for read_alignments in read2aware.values() {
//         for a in read_alignments {
//             let ctg = &aware_contigs[a.aware_id];
//             if a.is_ambiguous || !ctg.is_phased() {
//                 continue
//             }

//             let hid = aware_contigs[a.aware_id].haplotype_id().unwrap();
//             let hid_idx = hap_index[&hid];
//             let hid_sequence = &hap_sequences[hid_idx];

//             assert!(a.target_beg <= a.target_end);
//             assert!(hid.beg <= a.target_beg && a.target_end <= hid.end);

//             let length = std::cmp::max(a.query_end - a.query_beg, a.target_end-a.target_beg);

//             alignments.push( alignment::Alignment {
//                 query_idx: a.query_idx,
//                 query_len: a.query_len,
//                 query_beg: a.query_beg,
//                 query_end: a.query_end,
//                 strand: a.strand,
//                 target_idx: hid_idx,
//                 target_len: hid_sequence.len(),
//                 target_beg: a.target_beg - ctg.beg(),
//                 target_end: a.target_end - ctg.beg(),
//                 matches: 0,         // dummy, not supposed to be used
//                 mapping_length: 0,  // dummy, not supposed to be used
//                 mapq: 60,           // dummy, not supposed to be used
//                 length,
//                 identity: 0.0,      // dummy, not supposed to be used
//                 cigar: CigarString(Vec::new()),
//                 breaking_points: Vec::new()
//             });
//         }
//     }

//     let mut polisher = polish::Polisher::new(&hap_sequences, &read_db.sequences, &mut alignments);

//     spdlog::debug!("initializing polisher");
//     polisher.initialize();

//     spdlog::debug!("polishing windows");
//     let mut polished_sequences = polisher.polish();

//     hap_index.drain()
//         .map(|(hid, hap_idx)| (hid, std::mem::take(&mut polished_sequences[hap_idx])))
//         .collect()
// }
