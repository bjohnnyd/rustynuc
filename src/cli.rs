use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
pub(crate) struct Opt {
    #[structopt(
        short,
        long,
        help = "Represents paired end forward reads or all reads in case of SE",
        required = true,
        parse(from_os_str)
    )]
    fwd: Vec<PathBuf>,
    #[structopt(
        short,
        long,
        help = "Represents paired end reverse reads",
        parse(from_os_str)
    )]
    rev: Option<Vec<PathBuf>>,
}
