#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
//TODO: need to trim and map Sub-0-Rep-5 from run3 and compare to current
//TODO: Need to probably keep fishers test but deal differently and instead of FDR
//      make OR or similar, like GATK
//! QC tool for assesment of likelihood of oxo-G related variation in reads.

mod alignment;
mod cli;
mod error;

use crate::alignment::OxoPileup;
use rust_htslib::{bam, bam::Read};
use structopt::StructOpt;

/// Nucleotide alphabet used
pub const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

type Result<T> = std::result::Result<T, crate::error::Error>;

fn main() -> Result<()> {
    let opt = cli::RustyNuc::from_args();

    let fdr = opt.fdr;
    let mut bam = bam::Reader::from_path(opt.bam)?;
    let mut pileups = bam.pileup();
    pileups.set_max_depth(opt.max_depth);

    let mut oxo_pileups = Vec::new();
    for p in pileups {
        let pileup = p?;

        match &opt.positions {
            Some(positions) => {
                if positions.contains(&pileup.pos()) {
                    let oxo_pileup =
                        OxoPileup::new(pileup, Some(opt.min_number_variant), opt.quality);
                    println!("{}", &oxo_pileup);
                    if oxo_pileup.is_imbalanced(opt.fdr)? {
                        println!("Significantly Different");
                    }
                }
            }
            _ => {
                let oxo_pileup = OxoPileup::new(pileup, Some(opt.min_number_variant), opt.quality);
                if opt.all {
                    println!("{}", &oxo_pileup);
                } else {
                    oxo_pileups.push(oxo_pileup);
                }
            }
        }
    }

    if !oxo_pileups.is_empty() {
        let mut oxo_pileups = oxo_pileups
            .iter()
            .map(|oxo_pileup| (oxo_pileup, oxo_pileup.min_p_value()))
            .collect::<Vec<(&OxoPileup, f64)>>();
        oxo_pileups.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let tests = oxo_pileups.len() as f64;

        oxo_pileups
            .iter()
            .enumerate()
            .rev()
            .skip_while(|(i, (oxo_pileup, pval))| *pval > (fdr * (*i as f64 / tests)))
            .for_each(|(_, (oxo_pileup, _))| println!("{}", &oxo_pileup));
    }

    Ok(())
}
