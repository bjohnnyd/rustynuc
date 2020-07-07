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
/// Defualt visualization settings for track line
pub const TRACK_LINE: &str = "#coords 0";

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

    let alpha = opt.alpha;
    let mut seq_map = None;
    let mut id_map = HashMap::<u32, String>::new();

    if let Some(path) = opt.reference {
        let (rdr, _) = niffler::from_path(path)?;
        let fasta_rdr = bio::io::fasta::Reader::new(rdr);
        seq_map = Some(create_seq_map(fasta_rdr)?);
    }

    for p in pileups {
        let pileup = p?;
        let seq_name = String::from_utf8(header.tid2name(pileup.tid()).to_vec())?;
        let seq = match seq_map {
            Some(ref seq_map) => seq_map.get(&seq_name),
            None => None,
        };

        let oxo = OxoPileup::new(
            pileup,
            Some(opt.min_number_variant),
            opt.quality,
            opt.pseudocount,
            seq,
        );

        if !oxo.is_monomorphic() {
            match oxo.is_ref_g_or_c() {
                Some(true) | None => {
                    id_map.insert(oxo.ref_id, seq_name);
                    oxo_pileups.push(oxo);
                }
                _ => {}
            }
        }
    }

    // NOTE: If wanting ot use custom error will need to collect into tuple of (i, i.pval) but
    // could do this in parallel
    oxo_pileups.sort_by(|a, b| {
        a.min_p_value()
            .expect("Could not calculate min pval")
            .partial_cmp(&b.min_p_value().expect("Could not calculate min pval"))
            .expect("Could not compare the flot values")
    });

    let m = oxo_pileups.len();

    if opt.with_track_line {
        println!("{}", TRACK_LINE);
    }
    let _ = oxo_pileups
        .par_iter()
        .enumerate()
        .map(|(rank, p)| {
            let pval = p.min_p_value()?;
            let qval = (alpha * rank as f32) / m as f32;
            let sig = if pval < qval as f64 { 1 } else { 0 };
            let corrected = format!("{}\t{}", qval, sig);

            if let Some(ref reference) = seq_map {
                let seq_name = id_map.get(&p.ref_id).unwrap();
                let seq = reference.get(seq_name);
                match &seq {
                    Some(seq) => {
                        println!(
                            "{}\t{}\t{}",
                            p.to_bed_entry(Some(&id_map))?,
                            get_seq_context(seq, p.ref_pos),
                            corrected
                        );
                    }
                    None => {
                        warn!(
                        "The reference provided does not have record present in the bam file, {}",
                        &seq_name
                    );
                        println!("{}\t{}", p.to_bed_entry(Some(&id_map))?, corrected);
                    }
                }
            } else {
                println!("{}\t{}", p.to_bed_entry(Some(&id_map))?, corrected);
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

fn get_seq_context(seq: &str, pos: u32) -> String {
    let idx = (seq.len() - 1) as u32;
    match pos {
        0 => format!("X{}", &seq[0..2]),
        last_idx if last_idx == idx => {
            let last_idx = last_idx as usize;
            format!("{}X", &seq[(last_idx - 1)..(last_idx + 1)])
        }
        i => format!("{}", &seq[(i - 1) as usize..(i + 2) as usize]),
    }
}
