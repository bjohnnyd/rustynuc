use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
pub(crate) struct RustyNuc {
    #[structopt(
        short,
        long,
        help = "Alignments to correct for possible oxo-g damage",
        default_value = "2"
    )]
    pub(crate) min_number_variant: u32,
    #[structopt(short, long, help = "Minimum base quality to consider")]
    pub(crate) quality: Option<u8>,
    #[structopt(long, help = "Maximum pileup depth to use", default_value = "200")]
    pub(crate) max_depth: u32,
    #[structopt(long, help = "P-value signficance threshold", default_value = "0.05")]
    pub(crate) significance: f32,
    #[structopt(short, long, help = "Specific position to investigate")]
    pub(crate) positions: Option<Vec<u32>>,
    #[structopt(short, long, help = "Optional reference that might increase accuracy")]
    pub(crate) reference: Option<PathBuf>,
    #[structopt(short, long, help = "Whether to just print results for all positions")]
    pub(crate) all: bool,
    #[structopt(
        help = "Alignments to correct for possible oxo-g damage",
        required = true,
        parse(from_os_str)
    )]
    pub(crate) bam: PathBuf,
}
