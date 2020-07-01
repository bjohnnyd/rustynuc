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
use bio::io::fasta::Record;
use log::warn;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use structopt::StructOpt;

/// Nucleotide alphabet used
pub const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

type Result<T> = std::result::Result<T, crate::error::Error>;

fn main() -> Result<()> {
    let opt = cli::RustyNuc::from_args();
    env_logger::init();

    let fdr = opt.alpha;
    let mut bam = bam::Reader::from_path(opt.bam)?;
    let header = bam.header().clone();

    let records = match opt.reference {
        Some(fpath) => Some(read_reference(fpath)?),
        None => None,
    };

    let records = records
        .into_iter()
        .collect::<Vec<Record<Box<dyn std::io::Read>>>>();

    let seq_map = match records {
        Some(records) => {
            let mut seq_map = HashMap::<&[u8], &[u8]>::new();
            for record in records {
                let record = record?;
                seq_map.insert(record.id().as_bytes(), record.seq());
            }
            Some(seq_map)
        }
        None => None,
    };

    let mut pileups = bam.pileup();
    pileups.set_max_depth(opt.max_depth);

    let mut oxo_pileups = Vec::new();

    for p in pileups {
        let pileup = p?;
        oxo_pileups.push(OxoPileup::new(
            pileup,
            Some(opt.min_number_variant),
            opt.quality,
        ));
    }
    let mut oxo_pileups = oxo_pileups
        .iter()
        .map(|oxo_pileup| (oxo_pileup, oxo_pileup.min_p_value()))
        .collect::<Vec<(&OxoPileup, f64)>>();
    oxo_pileups.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let tests = oxo_pileups.len();

    let significant = tests
        - oxo_pileups
            .iter()
            .enumerate()
            .rev()
            .take_while(|(i, (oxo_pileup, pval))| *pval > (fdr * (*i as f64 / tests as f64)))
            .count();

    for (i, (p, pval)) in oxo_pileups.iter().enumerate() {
        let fdr_sig = i < significant;
        if let Some(ref reference) = seq_map {
            let seq_name = header.tid2name(p.ref_id);
            let seq = reference.get(&seq_name);
            match &seq {
                Some(seq) => {
                    let idx = (seq.len() - 1) as u32;
                    let context = match p.ref_pos {
                        0 => format!("X_{}", String::from_utf8(seq[0..0].to_vec())?),
                        last_idx if last_idx == idx => {
                            let last_idx = last_idx as usize;
                            format!("{}_X", String::from_utf8(seq[last_idx..last_idx].to_vec())?)
                        }
                        i => {
                            let left_idx = (i - 1) as usize;
                            let right_idx = (i + 1) as usize;
                            let left_nuc = String::from_utf8(seq[left_idx..left_idx].to_vec())?;
                            let right_nuc = String::from_utf8(seq[right_idx..right_idx].to_vec())?;
                            format!("{}_{}", left_nuc, right_nuc)
                        }
                    };

                    println!("{}\t{}\t{}\t{}", p, pval, fdr_sig, context);
                }
                None => {
                    warn!(
                        "The reference provided does not have record present in the bam file, {}",
                        String::from_utf8(seq_name.to_vec())?
                    );
                    println!("{}\t{}\t{}", p, pval, fdr_sig);
                }
            }
        } else {
            println!("{}\t{}\t{}", p, pval, fdr_sig);
        }
    }

    Ok(())
}

fn read_reference(
    path: std::path::PathBuf,
) -> Result<bio::io::fasta::Records<Box<dyn std::io::Read>>> {
    let (rdr, _) = niffler::from_path(path)?;
    let fasta = bio::io::fasta::Reader::new(rdr);
    Ok(fasta.records())
}
