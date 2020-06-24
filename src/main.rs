#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
//TODO: need to trim and map Sub-0-Rep-5 from run3 and compare to current
//! QC tool for assesment of likelihood of oxo-G related variation in reads.
//! Currently not applicable for multigenomic experiments due to assumptions made.

mod cli;

use oxo::alignment::OxoPileup;
use rust_htslib::{bam, bam::Read};
use structopt::StructOpt;

type Result<T> = std::result::Result<T, oxo::error::Error>;

static TEST_ISSUE_POS: [u32; 4] = [2977104, 2981956, 2982320, 2982432];

static TEST_BAM: &str = "/home/johnny/data/klaus_winzer/alignments/Sub-0-Rep-5_S1_LMERGED.bam";

fn main() -> Result<()> {
    let opt = cli::RustyNuc::from_args();
    println!("{:?}", opt);

    let mut bam = bam::Reader::from_path(opt.bam)?;
    let mut pileups = bam.pileup();
    pileups.set_max_depth(opt.max_depth);

    for p in pileups {
        let pileup = p?;
        let oxo_pileup = OxoPileup::new(pileup);

        if TEST_ISSUE_POS.contains(&oxo_pileup.ref_pos) {
            dbg!(&oxo_pileup);
        }
    }

    Ok(())
}
