use clap::{Parser, ValueEnum};

const CLI_HH_OTHERS: &str = "Other options";
const CLI_HH_ADVANCED: &str = "Advanced options";

#[derive(Parser)]
#[command(version)]
#[command(about = "Strainberry2: Strain-aware metagenome assembly with accurate long reads", long_about = None)]
pub struct Options {

    /// Input HiFi reads
    #[arg(long = "in-hifi", value_name = "PATH", conflicts_with = "in_ont")]
    pub in_hifi: Option<String>,
    
    /// Input ONT R10.4 reads
    #[arg(long = "in-ont", value_name = "PATH")]
    pub in_ont: Option<String>,

    /// Input assembly in FASTA format
    #[arg(short = 'r', long = "reference", value_name = "PATH")]
    pub reference: String,

    /// Output directory
    #[arg(short = 'o', long = "out-dir", value_name = "PATH")]
    pub output_dir: String,

    /// Minimum MAPQ value to consider a read alignment
    #[arg(short = 'q', long = "min-mapq", value_name = "NUM", default_value_t = 20)]
    pub min_mapq: u8,

    /// Maximum number of theads
    #[arg(short = 't', long = "threads", value_name = "NUM", default_value_t = 1)]
    pub nb_threads: usize,

    // ----------------
    // OTHER OPTIONS
    // ----------------

    /// Disable input de-replication
    #[arg(long = "no-derep", help_heading = CLI_HH_OTHERS)]
    pub no_derep: bool,

    /// Do not assembly, only phase and build strain-aware graph
    #[arg(long = "no-asm", help_heading = CLI_HH_OTHERS)]
    pub no_asm: bool,

    /// Do not delete temporary files
    #[arg(long = "keep-temp", help_heading = CLI_HH_OTHERS)]
    pub keep_temp: bool,

    /// Print debug information
    #[arg(long = "debug", help_heading = CLI_HH_OTHERS)]
    pub debug: bool,

    // ----------------
    // ADVANCED OPTIONS
    // ----------------

    /// Minimum number of phased variants to retain a haplotype
    #[arg(long = "min-snv", value_name = "NUM", default_value_t = 3, help_heading = CLI_HH_ADVANCED)]
    pub min_snv: usize,

    /// Minimum number of shared SNV positions between a read and a haplotype to consider a match
    #[arg(long = "min-shared-snv", value_name = "NUM", default_value_t = 1, help_heading = CLI_HH_ADVANCED)]
    pub min_shared_snv: usize,

    /// Minimum number of alternative-allele observations
    #[arg(long = "min-alt-count", value_name = "NUM", default_value_t = 5, help_heading = CLI_HH_ADVANCED)]
    pub min_alt_count: usize,

    /// Minimum fraction of alternative-allele observations
    #[arg(long = "min-alt-frac", value_name = "FLOAT", default_value_t = 0.125, help_heading = CLI_HH_ADVANCED)]
    pub min_alt_frac: f64,

    /// Minimum indel length
    #[arg(long = "min-indel", value_name = "NUM", default_value_t = 100, help_heading = CLI_HH_ADVANCED)]
    pub min_indel: usize,

    /// Minimum overhang length
    #[arg(long = "min-overhang", value_name = "NUM", default_value_t = 500, help_heading = CLI_HH_ADVANCED)]
    pub min_overhang: usize,

    /// Minimum aware-contig length
    #[arg(long = "min-aware-ctg-len", value_name = "NUM", default_value_t = 2000, help_heading = CLI_HH_ADVANCED)]
    pub min_aware_ctg_len: usize,    

    // HIDDEN ARGUMENTS

    /// Sequencing read technology
    #[arg(value_enum, long="mode", value_name="STR", default_value_t = Mode::Hifi, hide = true)]
    pub mode: Mode,

    /// Long-read alignment in BAM format
    #[arg(short = 'b', long = "bam", value_name = "PATH", hide = true)]
    pub bam: Option<String>,

    /// Lookback distance
    #[arg(short = 'l', long = "lookback", value_name = "NUM", default_value_t = 3000, hide = true)]
    pub lookback: usize,

    /// p-value for strain-bias computation
    #[arg(long = "p-value", value_name = "FLOAT", default_value_t = 0.005, hide = true)]
    pub strand_bias_pvalue: f64,

    /// Stop before phasing
    #[arg(long = "no-phase", hide = true)]
    pub no_phase: bool,

    /// Print more verbose debug information
    #[arg(long = "trace", hide = true)]
    pub trace: bool,

    // OLD DEPRECATED ARGUMENTS

    // /// User provided vcf file of SNV positions to consider
    // #[arg(short = 'v', long = "vcf", value_name = "PATH")]
    // pub vcf: Option<String>,

    // /// Do not split at putative misassembly events
    // #[arg(long = "no-split")]
    // pub no_split: bool,

    // /// Disable post-assembly polishing
    // #[arg(long = "no-polish")]
    // pub no_polish: bool,

    // /// Minimum QUAL value for loaded variants (effective only with --vcf)
    // #[arg(long = "min-var-qual", value_name = "NUM", default_value_t = 0)]
    // pub min_var_qual: usize,

    // /// Minimum SNV fraction to phase haplotypes in a sequence
    // #[arg(long = "min-snv-density", value_name = "NUM", default_value_t = 0.001)]
    // pub min_snv_density: f64,
}


#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Mode {
    Hifi,
    Nano,
}


// TODO: validate parsed options and possibly estimate other parameters
pub struct Config;

impl Config {

    pub fn from_options(_opts:Options) -> Config {
        todo!()
    }
}
