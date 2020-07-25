#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

//! QC tool for assesment of likelihood of 8-oxoG related variation.  
mod alignment;
mod cli;
mod error;

use crate::alignment::OxoPileup;
use alignment::ReadType;
use bio::utils::Interval;
use log::error;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use rust_htslib::{bcf, bcf::header::Id, bcf::Read as VcfRead};
use std::collections::{HashMap, HashSet};
use structopt::StructOpt;

// TODO:
// 1. Examples like 581774 where the primary alt is oxo-g but the secondary is real (might be
//    better to label with a different filter like MultiOxoG)

/// Nucleotide alphabet used
pub const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
/// Defualt visualization settings for track line
pub const TRACK_LINE: &str = "#coords 0";
/// Header information added to VCF
pub const HEADER_RECORDS: [&[u8]; 11] = [
    br#"##FILTER=<ID=OxoG,Description="OxoG two-sided p-value < 0.05">"#,
    br#"##FILTER=<ID=OccurenceInsufficient,Description="There is not a sufficient number of reads aligning in the FF and FR orientation">"#,
    br#"##INFO=<ID=OXO_DEPTH,Number=1,Type=Integer,Description="OxoG Pileup Depth">"#,
    br#"##INFO=<ID=ADENINE_FF_FR,Number=2,Type=Integer,Description="Adenine counts at FF and FR">"#,
    br#"##INFO=<ID=CYTOSINE_FF_FR,Number=2,Type=Integer,Description="Cytosine counts at FF and FR">"#,
    br#"##INFO=<ID=GUANINE_FF_FR,Number=2,Type=Integer,Description="Guanine counts at FF and FR">"#,
    br#"##INFO=<ID=THYMINE_FF_FR,Number=2,Type=Integer,Description="Thymine counts at FF and FR">"#,
    br#"##INFO=<ID=A_C_PVAL,Number=1,Type=Float,Description="A/C two-sided p-value">"#,
    br#"##INFO=<ID=G_T_PVAL,Number=1,Type=Float,Description="G/T two-sided p-value">"#,
    br#"##INFO=<ID=AF_FF_FR,Number=A,Type=Float,Description="Alternate frequency calculations on the FF and FR">"#,
    br#"##INFO=<ID=OXO_CONTEXT,Number=1,Type=String,Description="3mer Context of reference">"#,
];

/// Name for significant filter
pub const OXO_FILTER: &[u8] = b"OxoG";
/// Name for insufficient counts filter
pub const INSUFFICIENT_FILTER: &[u8] = b"OccurenceInsufficient";

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
    let mut variant_pos = HashSet::new();

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

    if let Some(ref path) = opt.bcf {
        let mut vcf = bcf::Reader::from_path(path)?;
        for record in vcf.records() {
            let record = record?;
            let reference = record.alleles()[0];
            if reference.len() == 1 && is_s(reference, 0) {
                variant_pos.insert(record.desc());
            }
        }
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

        if !variant_pos.is_empty() && !variant_pos.contains(&format!("{}:{}", &seq_name, pos)) {
            continue 'pileups;
        }
        let seq = seq_map.get(&seq_name);
        match seq {
            Some(seq) if !is_s(seq.as_bytes(), pos) => {}
            _ => {
                let oxo = OxoPileup::new(
                    seq_name,
                    pileup,
                    opt.min_reads,
                    opt.quality,
                    opt.pseudocount,
                    seq,
                );

                if opt.all {
                    oxo_pileups.push(oxo);
                } else if !oxo.is_monomorphic() && oxo.occurence_sufficient() {
                    oxo_pileups.push(oxo);
                }
            }
        };
    }

    if let Some(ref path) = opt.bcf {
        // TODO: Implement methods to
        // 1. Update VCF record if Oxo (filter and info fields)
        // 2. Update header with the Oxo Info
        let mut oxo_dict = HashMap::new();
        for oxo in oxo_pileups.into_iter() {
            oxo_dict.insert(oxo.desc(), oxo);
        }

        let mut vcf = bcf::Reader::from_path(path)?;
        let header = HEADER_RECORDS.iter().fold(
            bcf::header::Header::from_template(&vcf.header()),
            |mut header, header_record| {
                header.push_record(header_record);
                header
            },
        );

        let mut wrt = bcf::Writer::from_stdout(&header, true, bcf::Format::VCF)?;
        // let mut wrt_record = wrt.empty_record();
        let oxo_id = wrt.header().name_to_id(OXO_FILTER)?;
        let insufficient_filter = wrt.header().name_to_id(INSUFFICIENT_FILTER)?;

        wrt.set_threads(opt.threads)?;

        for record in vcf.records() {
            let mut record = record?;
            wrt.translate(&mut record);
            if let Some(oxo) = oxo_dict.get(&record.desc()) {
                update_vcf_record(
                    &mut record,
                    &oxo,
                    (oxo_id, insufficient_filter),
                    opt.fisher_sig,
                )?;
            }
            wrt.write(&record)?;
        }
    } else {
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
    }

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

fn update_vcf_record(
    record: &mut bcf::Record,
    oxo: &OxoPileup,
    filter_id: (Id, Id),
    sig_threshold: f64,
) -> Result<()> {
    let counts = oxo
        .nuc_counts(ReadType::FF)
        .iter()
        .zip(oxo.nuc_counts(ReadType::FR))
        .enumerate()
        .map(|(i, (ff_count, fr_count))| {
            let counts = [
                (ff_count - oxo.pseudocount) as i32,
                (fr_count - oxo.pseudocount) as i32,
            ];
            (NUCLEOTIDES[i], counts)
        })
        .collect::<HashMap<u8, [i32; 2]>>();

    let ac_pval = oxo.get_pval(b'C')?;
    let gt_pval = oxo.get_pval(b'G')?;
    let is_oxog = ac_pval < sig_threshold || gt_pval < sig_threshold;

    let context = match oxo.context {
        Some(ref seq) => seq,
        _ => "".as_bytes(),
    };

    let alt_alleles = record
        .alleles()
        .into_iter()
        .skip(1)
        .map(|allele| allele.to_owned())
        .collect::<Vec<Vec<u8>>>();

    let ref_allele = record.alleles()[0].to_owned();

    record.push_info_integer(b"OXO_DEPTH", &[oxo.ref_depth as i32])?;
    record.push_info_integer(
        b"ADENINE_FF_FR",
        counts.get(&b'A').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"CYTOSINE_FF_FR",
        counts.get(&b'C').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"GUANINE_FF_FR",
        counts.get(&b'G').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"THYMINE_FF_FR",
        counts.get(&b'T').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_float(b"A_C_PVAL", &[ac_pval as f32])?;
    record.push_info_float(b"G_T_PVAL", &[gt_pval as f32])?;
    record.push_info_string(b"OXO_CONTEXT", &[context])?;

    match ref_allele.as_slice() {
        b"G" | b"C" => {
            if is_oxog {
                record.push_filter(filter_id.0)
            };

            let oxo_af =
                record
                    .alleles()
                    .into_iter()
                    .skip(1)
                    .fold(Vec::new(), |mut oxo_af, alt_allele| {
                        let alt_counts = counts.get(&alt_allele[0]).unwrap_or_else(|| &[0, 0]);
                        let ref_counts = counts.get(&ref_allele[0]).unwrap_or_else(|| &[1, 1]);
                        let ff_af = alt_counts[0] as f32 / ref_counts[0] as f32;
                        let fr_af = alt_counts[1] as f32 / ref_counts[1] as f32;
                        oxo_af.push(ff_af);
                        oxo_af.push(fr_af);
                        oxo_af
                    });

            record.push_info_float(b"AF_FF_FR", oxo_af.as_slice())?;
        }
        _ => {}
    }

    if !oxo.occurence_sufficient() {
        record.push_filter(filter_id.1)
    }

    Ok(())
}
