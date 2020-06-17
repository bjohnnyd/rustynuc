#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

//! QC tool for assesment of likelihood of oxo-G related variation in reads.
//! Currently not applicable for multigenomic experiments due to assumptions made.

mod cli;

use structopt::StructOpt;

fn main() {
    let opt = cli::RustyNuc::from_args();
    println!("{:?}", opt);
}
