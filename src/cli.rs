use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
pub(crate) struct RustyNuc {
    /// Number of threads
    #[structopt(short, long, default_value = "4")]
    pub threads: usize,
    /// Determines verbosity of the processing, can be specified multiple times -vvv
    #[structopt(short, long, parse(from_occurrences))]
    pub verbosity: u8,
    #[structopt(
        short,
        long,
        help = "Minimum number of aligned reads in ff or fr orientation for a position to be considered",
        default_value = "4"
    )]
    pub(crate) min_reads: u32,
    #[structopt(
        short,
        long,
        help = "Minimum base quality to consider",
        default_value = "20"
    )]
    pub(crate) quality: u8,
    #[structopt(
        long,
        help = "AF above this cutoff in EITHER read orientation will be excluded from OxoAF filter",
        default_value = "0.25"
    )]
    pub(crate) oxo_af_ceiling: f32,
    #[structopt(long, help = "Maximum pileup depth to use", default_value = "1000")]
    pub(crate) max_depth: u32,
    #[structopt(
        long,
        help = "Significance threshold for Fisher's test",
        default_value = "0.05"
    )]
    pub(crate) fishers_sig: f64,
    #[structopt(long, help = "FDR threshold", default_value = "0.2")]
    pub(crate) alpha: f32,
    #[structopt(
        long,
        help = "AF on both the ff and fr at which point to not apply Fisher's exact pval filter in VCF",
        default_value = "0.1"
    )]
    pub(crate) fishers_af: f32,
    #[structopt(
        short,
        long,
        help = "Optional reference which will be used to determine sequence context and type of change"
    )]
    pub(crate) reference: Option<PathBuf>,
    #[structopt(
        short,
        long,
        help = "Whether to process and print information for every position in the BAM file"
    )]
    pub(crate) all: bool,
    #[structopt(
        help = "Alignments to investigate for possible 8-oxoG damage",
        required = true,
        parse(from_os_str)
    )]
    pub(crate) bam: PathBuf,
    #[structopt(
        short,
        long,
        help = "Whether to use pseudocounts (increments all counts by 1) when calculating statistics"
    )]
    pub(crate) pseudocount: bool,
    #[structopt(
        short,
        long,
        help = "Include track line (for correct visualization with IGV)"
    )]
    pub(crate) with_track_line: bool,
    #[structopt(long, help = "A BED file to restrict analysis to specific regions")]
    pub(crate) bed: Option<PathBuf>,
    #[structopt(short, long, help = "Skip calculating qvalue")]
    pub(crate) no_qval: bool,
    #[structopt(
        short,
        long,
        help = "BCF/VCF for variants called on the supplied alignment file"
    )]
    pub(crate) bcf: Option<PathBuf>,
    #[structopt(
        long,
        help = "Do not count overlapping mates when calculating total depth"
    )]
    pub(crate) no_overlapping: bool,
}

impl RustyNuc {
    pub fn set_logging(&self) {
        use log::LevelFilter::*;

        let log_level = match self.verbosity {
            level if level == 1 => Info,
            level if level == 2 => Debug,
            level if level > 2 => Trace,
            _ => Warn,
        };

        env_logger::builder()
            .format_module_path(false)
            .filter_module("rustynuc", log_level)
            .init();
    }
}
