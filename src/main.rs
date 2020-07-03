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
use log::warn;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use structopt::StructOpt;

/// Nucleotide alphabet used
pub const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

type Result<T> = std::result::Result<T, crate::error::Error>;

fn main() -> Result<()> {
    let opt = cli::RustyNuc::from_args();
    opt.set_logging();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .or_else(|_| Err(crate::error::Error::ThreadError))?;

    let mut bam = bam::Reader::from_path(opt.bam)?;
    let header = bam.header().clone();
    let mut pileups = bam.pileup();
    pileups.set_max_depth(opt.max_depth);
    let mut oxo_pileups = Vec::new();

    let fdr = opt.alpha;
    let mut seq_map = None;
    let mut id_map = HashMap::<u32, String>::new();

    if let Some(path) = opt.reference {
        let (rdr, _) = niffler::from_path(path)?;
        let fasta_rdr = bio::io::fasta::Reader::new(rdr);
        seq_map = Some(create_seq_map(fasta_rdr)?);
    }

    for p in pileups {
        let pileup = p?;
        let oxo = OxoPileup::new(pileup, Some(opt.min_number_variant), opt.quality);
        id_map.insert(
            oxo.ref_id,
            String::from_utf8(header.tid2name(oxo.ref_id).to_vec())?,
        );
        oxo_pileups.push(oxo);
    }

    let mut oxo_pileups = oxo_pileups
        .iter()
        .enumerate()
        .map(|(i, oxo_pileup)| (i, oxo_pileup, oxo_pileup.min_p_value()))
        .collect::<Vec<(usize, &OxoPileup, f64)>>();
    oxo_pileups.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    let tests = oxo_pileups.len();

    let significant = tests
        - oxo_pileups
            .iter()
            .enumerate()
            .rev()
            .take_while(|(i, (_, oxo_pileup, pval))| *pval > (fdr * (*i as f64 / tests as f64)))
            .count();

    let _ = oxo_pileups
        .par_iter()
        .map(|(i, p, pval)| {
            let fdr_sig = *i < significant;
            if let Some(ref reference) = seq_map {
                let seq_name = id_map.get(&p.ref_id).unwrap();
                // let seq_name = String::from_utf8(id_name.to_vec())?;
                let seq = reference.get(seq_name);
                // reference.iter().for_each(|(k, v)| println!("{}", k));
                match &seq {
                    Some(seq) => {
                        let idx = (seq.len() - 1) as u32;
                        let context = match p.ref_pos {
                            0 => format!("X{}", &seq[0..2]),
                            last_idx if last_idx == idx => {
                                let last_idx = last_idx as usize;
                                format!("{}X", &seq[(last_idx - 1)..(last_idx + 1)])
                            }
                            i => format!("{}", &seq[(i - 1) as usize..(i + 2) as usize]),
                        };

                        println!("{}\t{}\t{}\t{}", &seq_name, p, fdr_sig, context);
                    }
                    None => {
                        warn!(
                        "The reference provided does not have record present in the bam file, {}",
                        &seq_name
                    );
                        println!("{}\t{}\t{}", p.ref_id, p, fdr_sig);
                    }
                }
            } else {
                println!("{}\t{}\t{}", p.ref_id, p, fdr_sig);
            }
            Ok(())
        })
        .collect::<Result<Vec<()>>>();

    Ok(())
}

// TODO: Cannot return fasta::Record to share between threads
fn create_seq_map<T: std::io::Read>(
    rdr: bio::io::fasta::Reader<T>,
) -> Result<HashMap<String, String>> {
    let mut seq_map = HashMap::<String, String>::new();
    for record in rdr.records() {
        let record = record?;
        seq_map.insert(
            record.id().to_string(),
            String::from_utf8(record.seq().to_vec())?,
        );
    }
    Ok(seq_map)
}
