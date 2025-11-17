#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

use ahash::AHashMap as HashMap;
use anyhow::{bail,Context};
use clap::{crate_name, crate_version, Parser};
use itertools::Itertools;

use strainberry::alignment;
use strainberry::cli;
use strainberry::graph;
use strainberry::misassembly;
use strainberry::phase;
use strainberry::seq;
use strainberry::utils;
use strainberry::variant;


fn main() -> anyhow::Result<(), anyhow::Error> {
    
    let t_start = Instant::now();

    let opts = cli::Options::parse();

    if opts.trace {
        spdlog::default_logger().set_level_filter(spdlog::LevelFilter::All);
    } else if opts.debug {
        spdlog::default_logger().set_level_filter(spdlog::LevelFilter::MoreSevereEqual(spdlog::Level::Debug));
    }

    utils::check_dependencies(&["minimap2", "samtools"])?;

    rayon::ThreadPoolBuilder::new().num_threads(opts.nb_threads).build_global().unwrap();

    run_pipeline(opts)?;

    spdlog::info!("Time: {:.2}s | MaxRSS: {:.2}GB", t_start.elapsed().as_secs_f64(), utils::get_maxrss());

    Ok(())
}


fn run_pipeline(mut opts: cli::Options) -> anyhow::Result<(), anyhow::Error> {
    
    let reads_path = match (opts.in_hifi.as_ref(), opts.in_ont.as_ref()) {
        (None, None) => { bail!("Either --in-hifi or --in-ont is required. For more information, try '--help'.") },
        (Some(hifi_path), None) => {
            opts.mode = cli::Mode::Hifi;
            opts.strand_bias_pvalue = 0.005;
            Path::new(hifi_path)
        },
        (None, Some(ont_path)) => {
            opts.mode = cli::Mode::Nano;
            opts.strand_bias_pvalue = 0.01;
            Path::new(ont_path)
        },
        _ => unreachable!()
    };

    spdlog::info!("{} {}", crate_name!(), crate_version!());

    let mut reference_path = PathBuf::from_str(&opts.reference)?;
    let output_dir = Path::new(&opts.output_dir);

    // TODO: consider creating output directory in cli module (after option validation)
    // TODO: consider adding "--force" option and terminate if directory already exists
    fs::create_dir_all(output_dir).with_context(|| format!("Cannot create output directory: \"{}\"", output_dir.display()))?;

    // --------------
    // PREPROCESSING
    // --------------
    
    let preprocess_dir = output_dir.join("00-preprocess");
    std::fs::create_dir_all(&preprocess_dir)
        .with_context(|| format!("Cannot create output directory: \"{}\"", preprocess_dir.display()))?;

    if let Some("gfa") = reference_path.extension().and_then(std::ffi::OsStr::to_str) {
        spdlog::info!("Converting input GFA to FASTA");
        let gfa_graph = strainberry::graph::GfaGraph::from_gfa(&reference_path)?;
        let mut fasta_name = reference_path.file_stem().unwrap().to_os_string();
        fasta_name.push(".fasta");
        reference_path = preprocess_dir.join(fasta_name);
        gfa_graph.write_fasta(&reference_path, 0)?;
    }

    if !opts.no_derep {
        spdlog::info!("Purging duplications from input reference");
        reference_path = strainberry::derep::derep_assembly(&reference_path, &preprocess_dir, &opts)?;
        spdlog::info!("Dereplicated reference written to {}", reference_path.display());
    }

    let bam_path = if !opts.no_derep || opts.bam.is_none() {
        spdlog::info!("Mapping reads to dereplicated reference");
        let bam_path = preprocess_dir.join("alignment.bam");
        utils::run_minimap2(&reference_path, reads_path, &bam_path, &opts)?;
        bam_path
    } else {
        PathBuf::from_str(opts.bam.as_ref().unwrap())?
    };

    spdlog::info!("Loading sequences from: {}", reference_path.display());
    let ref_db = seq::SeqDatabase::build(&reference_path, true)?;
    spdlog::info!("{} sequences processed", ref_db.size());

    spdlog::info!("Loading reads");
    let read_db = seq::SeqDatabase::build(reads_path, true)?;
    let read_n75 = read_db.compute_nx(75).unwrap();
    // opts.lookback = (9 * read_n75) / 10;
    spdlog::debug!("{} sequences loaded / N75: {read_n75}", read_db.size());

    // ----------------
    // VARIANT CALLING
    // ----------------

    spdlog::info!("Calling variants from pileup");
    let variants = variant::load_variants_from_bam(&bam_path, &ref_db, &opts);
    let variant_path = preprocess_dir.join("variants.vcf");
    variant::write_variants_to_file(&variant_path, &variants, &ref_db)?;
    spdlog::debug!("{} variants identified", variants.values().map(|vars| vars.len()).sum::<usize>());
    let variant_info = preprocess_dir.join("variants.info.txt");
    variant::write_variants_info(&variant_info, &variants, &ref_db)?;

    // spdlog::info!("Filtering variants");
    // let variants = variant::filter_variants_by_density(variants, &ref_db, opts.min_snv_density);
    // let variant_path = preprocess_dir.join("variants.filtered.vcf");
    // variant::write_variants_to_file(&variant_path, &variants, &ref_db)?;
    // let variants = variant::filter_variants_hp(variants, &ref_db, 3);
    // let variant_path = preprocess_dir.join("variants.filtered.vcf");
    // variant::write_variants_to_file(&variant_path, &variants, &ref_db)?;
    // spdlog::debug!("{} variants retained", variants.values().map(|vars| vars.len()).sum::<usize>());

    spdlog::info!("Splitting reference at putative misjoins");
    let ref_intervals = misassembly::partition_reference(&bam_path, &ref_db, &read_db, &opts);
    spdlog::debug!("{} sequences after split", ref_intervals.len());

    spdlog::info!("Loading read alignments");
    let read_alignments = alignment::load_bam_alignments(&bam_path, &ref_db, &read_db, &opts);

    if opts.no_phase {
        return Ok(());
    }

    let ref_intervals = {
        spdlog::info!("Filtering low-coverage sequences");
        let mut ref_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &HashMap::new(), 0);
        strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut ref_contigs, &HashMap::new());
        ref_contigs.retain(|ctg| ctg.depth() >= opts.min_alt_count as f64);
        ref_contigs.into_iter().map(|ctg| ctg.interval()).collect_vec()
    };

    // ------------------
    // HAPLOTYPE PHASING
    // ------------------

    let phased_dir = output_dir.join("20-phased");
    spdlog::info!("Phasing strain haplotypes");
    let phaser = phase::Phaser::new(&bam_path, &ref_db, &read_db, &ref_intervals, phased_dir, &opts).unwrap();
    let phaser_result = phaser.phase(&variants);
    spdlog::info!("{} haplotypes phased", phaser_result.haplotypes.len());

    // ---------------------
    // STRAIN-AWARE CONTIGS
    // ---------------------

    spdlog::info!("Building strain-aware contigs");
    let mut aware_contigs = strainberry::awarecontig::build_aware_contigs(&ref_intervals, &phaser_result.haplotypes, opts.min_aware_ctg_len);

    spdlog::info!("Building succinct reads");
    let succinct_reads = seq::build_succinct_sequences(&bam_path, &ref_db, &read_db, &variants, &opts);

    spdlog::info!("Read realignment to haplotypes");
    let seq2haplo = strainberry::phase::separate_reads(&succinct_reads, &phaser_result.haplotypes, opts.min_shared_snv);

    spdlog::info!("Mapping reads to strain-aware contigs");
    let read2aware = strainberry::awarecontig::map_sequences_to_aware_contigs(&read_alignments, &mut aware_contigs, &seq2haplo);

    if opts.debug {

        let haplotypes_path = output_dir.join("haplotypes.txt");
        let mut writer = crate::utils::get_file_writer(&haplotypes_path);
        for ht_id in phaser_result.haplotypes.keys().sorted_unstable() {
            let ht = phaser_result.haplotypes.get(ht_id).unwrap();
            let ref_name = ref_db.names[ht_id.tid].as_str();
            let ref_beg = ht_id.beg;
            let ref_end = ht_id.end;
            let hid = ht_id.hid;
            writer.write_all(format!("{ref_name}_{ref_beg}-{ref_end}_h{hid}\t{}\n", ht.seq_string()).as_bytes())?;
        }

        let read2aware_path = output_dir.join("read2aware.txt");
        let mut writer = crate::utils::get_file_writer(&read2aware_path);
        for (read_idx, read_alignments) in read2aware.iter() {
            writer.write_all(format!("{}/{read_idx}", read_db.names[*read_idx]).as_bytes())?;
            for a in read_alignments {
                let ctg = &aware_contigs[a.aware_id];
                let ref_name = ref_db.names[ctg.tid()].as_str();
                let ctg_beg = ctg.beg();
                let ctg_end = ctg.end();
                let ctg_suf = if let Some(hid) = ctg.hid() { format!("_h{hid}") } else { String::new() };
                let ctg_name = format!("{ref_name}_{ctg_beg}-{ctg_end}{ctg_suf}");
                let qbeg = a.query_beg;
                let qend = a.query_end;
                let qlen = a.query_len;
                let dist = a.dist;
                let nb_shared_snvs = a.nb_shared_snvs;
                let is_ambiguous = a.is_ambiguous();
                writer.write_all(format!("\t({ctg_name},{qbeg},{qend},{qlen},{dist}/{nb_shared_snvs},{is_ambiguous})").as_bytes())?;
            }
            writer.write_all(b"\n")?;
        }
    }

    drop(phaser_result);

    // -------------------
    // STRAIN-AWARE GRAPH
    // -------------------

    let graphs_dir = output_dir.join("40-graphs");
    fs::create_dir_all(graphs_dir.as_path()).with_context(|| format!("Cannot create graphs directory: \"{}\"", graphs_dir.display()))?;
    
    spdlog::info!("Building strain-aware graph");
    let mut aware_graph = graph::AwareGraph::build(&aware_contigs);
    aware_graph.add_edges_from_aware_alignments(&read2aware);
    aware_graph.write_gfa(graphs_dir.join("aware_graph.raw.gfa"), &ref_db)?;

    spdlog::info!("Strain-aware graph resolution");
    aware_graph.remove_self_loops();
    aware_graph.remove_weak_edges(opts.min_alt_count);
    aware_graph.add_bridges_2(&read2aware);
    aware_graph.remove_weak_transitive_biedges();
    aware_graph.write_gfa(graphs_dir.join("aware_graph.gfa"), &ref_db)?;
    
    aware_graph.resolve_junctions(opts.min_alt_count);
    aware_graph.remove_weak_transitive_biedges();
    aware_graph.write_gfa(graphs_dir.join("aware_graph.resolved.gfa"), &ref_db)?;

    if opts.no_asm {
        return Ok(());
    }

    // ---------------
    // ASSEMBLY GRAPH
    // ---------------

    spdlog::info!("Building assembly graph");
    let unitig_dir = output_dir.join("50-unitigs");
    let unitig_graph = aware_graph.build_assembly_graph(&ref_db, &read_db, phaser.fragments_dir(), &unitig_dir, &opts)?;
    unitig_graph.write_gfa(&output_dir.join("assembly.gfa"))?;
    let assembly_fasta_path = output_dir.join("assembly.fasta");
    unitig_graph.write_fasta(&assembly_fasta_path, 100)?;
    spdlog::info!("Final assembly written to {}", assembly_fasta_path.display());
    
    Ok(())
}
