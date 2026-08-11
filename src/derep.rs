// use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::str::FromStr;

use anyhow::{bail,Context, Result};
use ahash::{AHashMap as HashMap, AHashSet as HashSet};
use itertools::Itertools;
// use flate2::read::MultiGzDecoder;
// use itertools::Itertools;
use needletail::Sequence;

use crate::cli::Options;
use crate::alignment::PafAlignment;


type MatchInterval = (usize,usize,String);

fn compute_matching_intervals(fasta_path: &Path, opts: &Options) -> Result<HashMap<String,Vec<MatchInterval>>> {

    let fasta_path_str = fasta_path.to_str().unwrap();

    // let args = ["-t", &opts.nb_threads.to_string(), "-cDP", "--dual=no", "--no-long-join", "-r85", fasta_path_str, fasta_path_str];
    // let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", fasta_path_str, fasta_path_str];
    let mm2_args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];

    let mut minimap2 = Command::new("minimap2")
        .args(mm2_args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .context("cannot run minimap2")?;

    let mm2_stdout = minimap2.stdout.take()
        .context("Could not capture minimap2 standard output.")?;

    let reader = BufReader::new(mm2_stdout);
    let alignments: Vec<PafAlignment> = reader.lines()
        .map_while(Result::ok)
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .map(|line| PafAlignment::from_str(&line).with_context(||format!("error parsing line:\n{line}")).unwrap())
        .collect();

    let mm2_status = minimap2.wait().context("minimap2 did not run successfully")?;
    if !mm2_status.success() {
        bail!("minimap2 exited with {mm2_status}");
    }

    let mut matching_intervals: HashMap<String,Vec<(usize,usize,String)>> = HashMap::new();
    for a in alignments.into_iter().filter(|a| a.query_name != a.target_name) {
        
        // let map_type = a.map_type(100, 0.5);
        assert!(a.mapping_length > 0);
        
        // if a.mapping_length < 1000 && !map_type.is_containment() {
        //     continue
        // }

        if (a.query_end-a.query_beg).min(a.target_end-a.target_beg) < 2000 || a.identity < 90.0 {
            continue
        }

        // let (query_beg, query_end) = if a.strand == b'+' { (a.query_beg,a.query_end) } else { (a.query_beg,a.query_length) };
        matching_intervals.entry(a.query_name.clone()).or_default().push((a.query_beg, a.query_end, a.target_name.clone()));
        matching_intervals.entry(a.target_name).or_default().push((a.target_beg, a.target_end, a.query_name));

        // match map_type {
        //     MappingType::QueryContained => {
        //         matching_intervals.entry(a.query_name).or_default().push((0,a.query_length,a.target_name));
        //     },
        //     MappingType::QueryPrefix => {
        //         let (query_beg, query_end) = if a.strand == b'+' { (0,a.query_end) } else { (a.query_beg,a.query_length) };
        //         // println!("{} => {map_type:?} | {} | ({},{},{})", a.query_name, a.strand as char, query_beg, query_end, a.target_name);
        //         matching_intervals.entry(a.query_name).or_default().push((query_beg, query_end,a.target_name));
        //     },
        //     MappingType::QuerySuffix => {
        //         let (query_beg, query_end) = if a.strand == b'+' { (a.query_beg,a.query_length) } else { (0,a.query_end) };
        //         // println!("{} => {map_type:?} | {} | ({},{},{})", a.query_name, a.strand as char, query_beg, query_end, a.target_name);
        //         matching_intervals.entry(a.query_name).or_default().push((query_beg, query_end,a.target_name));
        //     },
        //     MappingType::ReferenceContained => {
        //         matching_intervals.entry(a.target_name).or_default().push((0,a.target_length,a.query_name));
        //     },
        //     MappingType::ReferencePrefix => {
        //         // println!("{} => {map_type:?} | {} | ({},{},{})", a.target_name, a.strand as char, 0,a.target_end, a.query_name);
        //         matching_intervals.entry(a.target_name).or_default().push((0,a.target_end,a.query_name));
        //     },
        //     MappingType::ReferenceSuffix => {
        //         // println!("{} => {map_type:?} | {} | ({},{},{})", a.target_name, a.strand as char, a.target_beg,a.target_length, a.query_name);
        //         matching_intervals.entry(a.target_name).or_default().push((a.target_beg,a.target_length,a.query_name));
        //     },
        //     MappingType::DovetailPrefix => {
        //         let (query_beg, query_end) = if a.strand == b'+' { (0,a.query_end) } else { (a.query_beg,a.query_length) };
        //         matching_intervals.entry(a.query_name.clone()).or_default().push((query_beg, query_end,a.target_name.clone()));
        //         matching_intervals.entry(a.target_name).or_default().push((a.target_beg, a.target_length, a.query_name));
        //     },
        //     MappingType::DovetailSuffix => {
        //         let (query_beg, query_end) = if a.strand == b'+' { (a.query_beg,a.query_length) } else { (0,a.query_end) };
        //         matching_intervals.entry(a.query_name.clone()).or_default().push((query_beg, query_end,a.target_name.clone()));
        //         matching_intervals.entry(a.target_name).or_default().push((0, a.target_end, a.query_name));
        //     },
        //     MappingType::Internal => {}
        // }
    }

    Ok(matching_intervals)
}


pub fn derep_assembly(fasta_path: &Path, work_dir: &Path, opts: &Options) -> Result<PathBuf> {
    assert!(fasta_path.is_file());

    let tmp_dir = work_dir.join("tmp");
    std::fs::create_dir_all(&tmp_dir)
        .with_context(|| format!("Cannot create output directory: \"{}\"", tmp_dir.display()))?;

    let output_path = work_dir.join("assembly_derep.fasta");
    let mut asm_path = fasta_path.to_path_buf();
    
    for it in 1usize.. {

        // println!("  derep iteration #{it} on \"{}\"", asm_path.display());

        let mut modified = false;
        let matching_intervals = self::compute_matching_intervals(&asm_path, opts)?;

        let mut seq_dict: HashMap<String,String> = HashMap::new();
        let mut asm_reader = needletail::parse_fastx_file(asm_path).context("Cannot open fasta file")?;
        while let Some(record) = asm_reader.next() {
            let record = record.unwrap();
            let name = record.id().split(|b| b.is_ascii_whitespace()).next().unwrap();
            let name = std::str::from_utf8(name).unwrap().to_string();
            let sequence = std::str::from_utf8(&record.normalize(false)).unwrap().to_string();
            seq_dict.insert(name, sequence);
        }

        let mut derep_ctg_id: usize = 0;
        let derep_path = tmp_dir.join(format!("derep_{it}.fasta"));
        let mut derep_writer = crate::utils::get_file_writer(&derep_path);
        let mut processed: HashSet<&str> = HashSet::new();
        for (seq_id, sequence) in seq_dict.iter().sorted_by_key(|(_,seq)| seq.len()) {
            
            let mut seq_intervals = Vec::new();
            
            if let Some(intervals) = matching_intervals.get(seq_id.as_str()) {
                let sorted_intervals = intervals.iter()
                    .filter_map(|(beg,end,name)| {
                        if !processed.contains(name.as_str()) { Some((*beg,*end,name.clone())) } else { None }
                    }).sorted().collect_vec();
                
                let mut merged_intervals = Vec::with_capacity(sorted_intervals.len());
                for (beg,end,_) in sorted_intervals {
                    if let Some((_last_beg,last_end)) = merged_intervals.last_mut() {
                        if beg <= *last_end {
                            *last_end = std::cmp::max(*last_end, end);
                            continue
                        }
                    }
                    merged_intervals.push((beg,end));
                }

                modified = modified || !merged_intervals.is_empty();

                let mut beg = 0;
                for &(mi_beg, mi_end) in &merged_intervals {
                    seq_intervals.push((beg, mi_beg));
                    beg = mi_end;
                }
                seq_intervals.push((beg, sequence.len()));
                processed.insert(seq_id.as_str());
            } else { // no matching interval: output the full sequence
                seq_intervals.push((0,sequence.len()));
            }

            // Sequences are processed in increasing length and only alignments involving the suffix/prefix are retained
            // This means I cannot get more than 1 interval
            // assert!(seq_intervals.len() <= 1);
            
            for &(beg,end) in seq_intervals.iter().filter(|(beg,end)| end-beg >= 2000) {
                let header = format!(">{derep_ctg_id}\n");
                derep_writer.write_all(header.as_bytes())?;
                let sequence = crate::utils::insert_newlines(&sequence[beg..end], 120);
                derep_writer.write_all(sequence.as_bytes())?;
                derep_writer.write_all(b"\n")?;
                derep_ctg_id += 1;
            }
        }

        if !modified {
            std::fs::rename(&derep_path, &output_path)
                .with_context(|| format!("cannot move file \"{}\" to \"{}\"", derep_path.display(), output_path.display()))?;
            break
        }

        asm_path = derep_path;
    }

    if !opts.keep_temp {
        std::fs::remove_dir_all(&tmp_dir)?;
    }

    Ok(output_path)
}

