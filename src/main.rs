#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

//! QC tool for assesment of likelihood of 8-oxoG related variation.  
mod alignment;
mod cli;
mod error;

use crate::alignment::OxoPileup;
use bio::utils::Interval;
use log::error;
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
    let alpha = opt.alpha;
    let calc_qval = !opt.no_qval;

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .or_else(|_| Err(crate::error::Error::ThreadError))?;

    let mut bam = bam::Reader::from_path(opt.bam)?;
    let header = bam.header().clone();
    let pileups = create_pileups(&mut bam, opt.threads, opt.max_depth)?;
    let mut oxo_pileups = Vec::new();

    let mut seq_map = HashMap::new();
    let mut bed_map = HashMap::new();

    if let Some(ref path) = opt.reference {
        let (rdr, _) = niffler::from_path(path)?;
        let fasta_rdr = bio::io::fasta::Reader::new(rdr);
        seq_map = create_seq_map(fasta_rdr)?;
    }

    if let Some(ref path) = opt.bed {
        let (rdr, _) = niffler::from_path(path)?;
        let bed_rdr = bio::io::bed::Reader::new(rdr);
        create_bed_map(bed_rdr, &mut bed_map)?;
    }

    'pileups: for p in pileups {
        let pileup = p?;
        let pos = pileup.pos() as usize;
        let seq_name = String::from_utf8(header.tid2name(pileup.tid()).to_vec())?;
        if let Some(regions) = bed_map.get(&seq_name) {
            if !regions
                .par_iter()
                .any(|region| region.contains(&(pos as u64)))
            {
                continue 'pileups;
            }
        }
        let seq = seq_map.get(&seq_name);
        match seq {
            Some(seq) if !is_s(seq.as_bytes(), pos) => {}
            _ => {
                let oxo = OxoPileup::new(
                    seq_name,
                    pileup,
                    Some(opt.min_reads),
                    opt.quality,
                    opt.pseudocount,
                    seq,
                );

                if opt.all {
                    oxo_pileups.push(oxo);
                } else if !oxo.is_monomorphic() && oxo.occurence_sufficient(opt.min_reads) {
                    oxo_pileups.push(oxo);
                }
            }
        };
    }

    if calc_qval {
        oxo_pileups.sort_by(|a, b| {
            a.min_p_value()
                .expect("Could not calculate min pval")
                .partial_cmp(&b.min_p_value().expect("Could not calculate min pval"))
                .expect("Could not compare the float values")
        });
    }

    let m = if opt.reference.is_some() {
        oxo_pileups.len() * 2
    } else {
        oxo_pileups.len()
    };

    if opt.with_track_line {
        println!("{}", TRACK_LINE);
    }
    let _ = oxo_pileups
        .par_iter()
        .enumerate()
        .map(|(rank, p)| {
            let pval = p.min_p_value()?;
            let mut corrected = String::from("");
            if calc_qval {
                let qval = (alpha * rank as f32) / m as f32;
                let sig = if pval < qval as f64 { 1 } else { 0 };
                corrected = format!("{}\t{}", qval, sig);
            }
            println!("{}\t{}", p.to_bed_entry()?, corrected);
            Ok(())
        })
        .collect::<Result<Vec<()>>>();

    Ok(())
}

fn create_pileups<'a>(
    bam: &'a mut bam::Reader,
    threads: usize,
    max_depth: u32,
) -> Result<bam::pileup::Pileups<'a, bam::Reader>> {
    bam.set_threads(threads)?;

    let mut pileups = bam.pileup();
    pileups.set_max_depth(max_depth);
    Ok(pileups)
}

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

fn create_bed_map<T: std::io::Read>(
    mut rdr: bio::io::bed::Reader<T>,
    map: &mut HashMap<String, Vec<Interval<u64>>>,
) -> Result<()> {
    for (i, record) in rdr.records().enumerate() {
        let record = record.or_else(|_| Err(crate::error::Error::BedRecordError(i + 1)))?;
        let interval = Interval::new(std::ops::Range {
            start: record.start(),
            end: record.end(),
        })
        .or_else(|_| {
            Err(crate::error::Error::IncorrectInterval(
                i,
                record.start(),
                record.end(),
            ))
        })?;
        let regions = map.entry(record.chrom().to_string()).or_default();
        regions.push(interval);
    }
    Ok(())
}

fn is_s(seq: &[u8], pos: usize) -> bool {
    if pos >= seq.len() {
        error!("Reference sequence is shorter than BAM alignment positions");
        std::process::exit(1)
    } else {
        match seq[pos] {
            b'G' | b'C' | b'g' | b'c' => true,
            _ => false,
        }
    }
}
