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
    #[structopt(long, help = "Maximum pileup depth to use", default_value = "400")]
    pub(crate) max_depth: u32,
    #[structopt(
        long,
        help = "Significance threshold for Fisher's test",
        default_value = "0.05"
    )]
    pub(crate) fisher_sig: f64,
    #[structopt(long, help = "FDR threshold", default_value = "0.2")]
    pub(crate) alpha: f32,
    #[structopt(
        short,
        long,
        help = "Optional reference which will be used to determine sequence context and type of change"
    )]
    pub(crate) reference: Option<PathBuf>,
    #[structopt(short, long, help = "Whether to just print results for all positions")]
    pub(crate) all: bool,
    #[structopt(
        help = "Alignments to correct for possible 8-oxoG damage",
        required = true,
        parse(from_os_str)
    )]
    pub(crate) bam: PathBuf,
    #[structopt(
        short,
        long,
        help = "Whether to use pseudocounts (adds +1 to all counts) when calculating statistics"
    )]
    pub(crate) pseudocount: bool,
    #[structopt(
        short,
        long,
        help = "Include track line (for correct visualization with IGV)"
    )]
    pub(crate) with_track_line: bool,
    #[structopt(
        short,
        long,
        help = "A BED file to restrict analysis to specific regions"
    )]
    pub(crate) bed: Option<PathBuf>,
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
