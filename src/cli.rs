use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
pub(crate) struct RustyNuc {
    #[structopt(
        short,
        long,
        help = "Alignments to correct for possible oxo-g damage",
        default_value = "4"
    )]
    pub(crate) min_number_variant: usize,
    #[structopt(short, long, help = "Minimum base quality to consider")]
    pub(crate) base_quality: Option<usize>,
    #[structopt(
        short,
        long,
        help = "Maximum pileup depth to use",
        default_value = "200"
    )]
    pub(crate) max_depth: u32,
    #[structopt(short, long, help = "Specific position to investigate")]
    pub(crate) positions: Option<Vec<u32>>,
    #[structopt(short, long, help = "Optional reference that might increase accuracy")]
    pub(crate) reference: Option<PathBuf>,
    #[structopt(
        help = "Alignments to correct for possible oxo-g damage",
        required = true,
        parse(from_os_str)
    )]
    pub(crate) bam: PathBuf,
}
