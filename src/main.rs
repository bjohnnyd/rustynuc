#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]

//! QC tool for assesment of likelihood of 8-oxoG related variation.
mod alignment;
mod cli;
mod error;
mod genomic;
mod vcf;

use crate::alignment::{is_pos_iupac_s, OxoPileup};
use genomic::{create_fasta_records, create_pileups, create_regions};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use rust_htslib::{bcf, bcf::Read as VcfRead};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::{io, sync::mpsc::channel};
use structopt::StructOpt;
use vcf::{
    apply_af_too_low_filter, apply_fishers_filter, update_vcf_record, AF_MIN, HEADER_RECORDS,
    INSUFFICIENT_FILTER,
};

/// Nucleotide alphabet used
pub const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
/// Default visualization settings for track line
pub const TRACK_LINE: &str = "#coords 0";

type Result<T> = std::result::Result<T, crate::error::Error>;

fn main() -> Result<()> {
    match main_try() {
        Ok(()) => Ok(()),
        Err(e) => {
            error!("{}", e);
            Err(e)
        }
    }
}

fn main_try() -> Result<()> {
    let opt = cli::RustyNuc::from_args();
    opt.set_logging();
    let alpha = opt.alpha;
    let calc_qval = !opt.no_qval;
    let stdout = io::stdout();
    let lock = stdout.lock();
    let mut wstdout = io::BufWriter::new(lock);

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .or_else(|_| Err(crate::error::Error::ThreadError))?;

    let mut bam = bam::Reader::from_path(opt.bam)?;
    let header = bam.header().clone();
    let pileups = create_pileups(&mut bam, opt.threads, opt.max_depth)?;
    let mut oxo_pileups = Vec::new();

    let mut fasta_records = HashMap::new();
    let mut regions = HashMap::new();
    let mut vcf_positions = HashSet::new();
    let mut oxo_check_locations = HashMap::new();

    if let Some(ref path) = opt.reference {
        info!("Reading reference...");
        let (rdr, _) = niffler::from_path(path)?;
        let fasta_rdr = bio::io::fasta::Reader::new(rdr);
        fasta_records = create_fasta_records(fasta_rdr)?;
        debug!(
            "Read in the following chromosomes: {}",
            fasta_records
                .keys()
                .cloned()
                .collect::<Vec<String>>()
                .join(";")
        );
    }

    if let Some(ref path) = opt.bed {
        info!("Reading regions...");
        let (rdr, _) = niffler::from_path(path)?;
        let bed_rdr = bio::io::bed::Reader::new(rdr);
        create_regions(bed_rdr, &mut regions)?;
    }

    if let Some(ref path) = opt.bcf {
        info!("Reading VCF...");
        let mut vcf = bcf::Reader::from_path(path)?;
        for (i, record) in vcf.records().enumerate() {
            debug!("Accesing {} record in the VCF...", i + 1);
            let record = record?;
            let reference = record.alleles()[0];
            debug!(
                "Succesfully obtained description {} with reference allele {}",
                record.desc(),
                String::from_utf8(reference.to_vec())?
            );
            if reference.len() == 1 && is_pos_iupac_s(reference, 0) {
                debug!(
                    "Record {} in the VCF will be processed for potential 8-oxoG",
                    record.desc()
                );
                vcf_positions.insert(record.desc());
            }
        }
    }

    info!("Started pileup analysis...");
    'pileups: for p in pileups {
        let pileup = p?;
        let pos = pileup.pos() as usize;
        let chrom = String::from_utf8(header.tid2name(pileup.tid()).to_vec())?;
        debug!(
            "Processing pileup on chromosome {} at position {}...",
            &chrom, pos
        );
        if let Some(chrom_regions) = regions.get(&chrom) {
            if !chrom_regions
                .par_iter()
                .any(|region| region.contains(&(pos as u64)))
            {
                debug!(
                    "Pileup {}:{} is outside of the provided BED regions and will be skipped",
                    &chrom, pos
                );
                continue 'pileups;
            }
        }

        if !vcf_positions.is_empty() && !vcf_positions.contains(&format!("{}:{}", &chrom, pos)) {
            debug!(
                "Pileup on chromosome {} at position {} is not amongst the variants in the VCF and will be skipped",
                &chrom, pos
            );

            continue 'pileups;
        }

        let seq = fasta_records.get(&chrom);

        if let Some(seq) = seq {
            if !is_pos_iupac_s(seq.as_bytes(), pos) {
                debug!("Pileup {}:{} is not G/C and will be skipped", &chrom, pos);
                continue 'pileups;
            }
        }

        let oxo = OxoPileup::new(
            chrom,
            pileup,
            opt.min_reads,
            opt.quality,
            opt.pseudocount,
            seq,
            opt.no_overlapping,
        );

        if let Some(_) = opt.bcf {
            oxo_check_locations.insert(oxo.desc(), oxo);
        } else if opt.all {
            debug!("Pileup {}:{} will be analysed for 8-oxoG", &oxo.chrom, pos);
            oxo_pileups.push(oxo);
        } else if !oxo.is_monomorphic() && oxo.occurence_sufficient() {
            debug!("Pileup {}:{} will be analysed for 8-oxoG", &oxo.chrom, pos);
            oxo_pileups.push(oxo);
        }
    }

    if let Some(ref path) = opt.bcf {
        let mut vcf = bcf::Reader::from_path(path)?;

        let header = HEADER_RECORDS.iter().fold(
            bcf::header::Header::from_template(&vcf.header()),
            |mut header, header_record| {
                header.push_record(header_record);
                header
            },
        );

        let mut wrt = bcf::Writer::from_stdout(&header, true, bcf::Format::VCF)?;

        wrt.set_threads(opt.threads)?;

        info!("Started processing VCF file variants...");
        for record in vcf.records() {
            let mut record = record?;
            wrt.translate(&mut record);
            let ref_allele = record.alleles()[0].to_ascii_uppercase()[0];
            debug!(
                "Checking if record {} is among the total of {} records labelled as to be processed for OxoG...",
                &record.desc(),
                oxo_check_locations.len(),
            );
            match oxo_check_locations.get(&record.desc()) {
                Some(oxo)
                    if !oxo.is_monomorphic()
                        && oxo.nuc_sufficient(ref_allele).iter().all(|x| *x) =>
                {
                    debug!("Updating VCF information for {}", record.desc());
                    if let Some(g2t_or_c2a_max_af) = update_vcf_record(&mut record, &oxo)? {
                        if g2t_or_c2a_max_af > opt.af_either_pass {
                            debug!("Record {} has a max orientation AF of {} which is above the frequency of {} and will not be considered for OxoG filtering", record.desc(), g2t_or_c2a_max_af, opt.af_either_pass)
                        } else {
                            if !oxo.occurence_sufficient() {
                                debug!("Record {} has insufficient count", record.desc());
                                let insufficient_filter =
                                    wrt.header().name_to_id(INSUFFICIENT_FILTER)?;
                                record.push_filter(insufficient_filter);
                            } else {
                                if !opt.skip_fishers {
                                    debug!(
                                        "Checking if record {} passes Fisher's OxoG filter...",
                                        record.desc()
                                    );
                                    apply_fishers_filter(
                                        &mut record,
                                        opt.fishers_sig as f32,
                                        opt.af_both_pass,
                                    )?;
                                }
                                apply_af_too_low_filter(&mut record, AF_MIN)?;
                            }
                        }
                    }
                }
                _ => {}
            }
            wrt.write(&record)?
        }
    } else {
        if calc_qval {
            info!("Ranking p-values...");
            oxo_pileups.sort_by(|a, b| {
                a.pval()
                    .expect("Could not calculate min pval")
                    .partial_cmp(&b.pval().expect("Could not calculate min pval"))
                    .expect("Could not compare the float values")
            });
        }

        let m = if opt.reference.is_some() {
            oxo_pileups.len() * 2
        } else {
            oxo_pileups.len()
        };

        if opt.with_track_line {
            writeln!(wstdout, "{}", TRACK_LINE)?;
        }

        info!("Started adj.p-value calculation and output processing...");
        let (tx, rx) = channel();
        let _ = oxo_pileups
            .par_iter()
            .enumerate()
            .map_with(tx, |s, (rank, p)| {
                let pval = p.pval()?;
                let mut corrected = String::from("");
                if calc_qval {
                    let qval = (alpha * rank as f32) / m as f32;
                    let sig = if pval < qval as f64 { 1 } else { 0 };
                    corrected = format!("{}\t{}", qval, sig);
                }
                s.send(format!("{}\t{}", p.to_bed_row()?, corrected))?;
                Ok(())
            })
            .collect::<Result<Vec<()>>>();
        for _ in 0..oxo_pileups.len() {
            writeln!(wstdout, "{}", rx.recv()?)?;
        }
    }

    info!("Done!");
    Ok(())
}
