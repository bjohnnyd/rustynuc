use crate::Result;
use bio::utils::Interval;
use log::{debug, error, info};
use rust_htslib::bcf::Read as VcfRead;
use rust_htslib::{bam, bam::Read, bcf};
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
};

pub fn create_pileups(
    bam: &mut bam::Reader,
    threads: usize,
    max_depth: u32,
) -> Result<bam::pileup::Pileups<'_, bam::Reader>> {
    bam.set_threads(threads)?;

    let mut pileups = bam.pileup();
    pileups.set_max_depth(max_depth);
    Ok(pileups)
}

pub fn create_fasta_records<T: std::io::Read>(
    rdr: bio::io::fasta::Reader<T>,
) -> Result<HashMap<String, String>> {
    let mut fasta_records = HashMap::<String, String>::new();
    for record in rdr.records() {
        let record = record?;
        fasta_records.insert(
            record.id().to_string(),
            String::from_utf8(record.seq().to_vec())?,
        );
    }
    Ok(fasta_records)
}

pub fn create_regions<T: std::io::Read>(
    mut rdr: bio::io::bed::Reader<T>,
    map: &mut HashMap<String, Vec<Interval<u64>>>,
) -> Result<()> {
    for (i, record) in rdr.records().enumerate() {
        let record = record.map_err(|_| crate::error::Error::BedRecordError(i + 1))?;
        let interval = Interval::new(std::ops::Range {
            start: record.start(),
            end: record.end(),
        })
        .map_err(|_| crate::error::Error::IncorrectInterval(i, record.start(), record.end()))?;
        let regions = map.entry(record.chrom().to_string()).or_default();
        regions.push(interval);
    }
    Ok(())
}

pub fn create_vcf_positions(path: &PathBuf) -> Result<HashSet<String>> {
    info!("Reading VCF...");
    let mut vcf = bcf::Reader::from_path(path)?;

    vcf.records()
        .enumerate()
        .fold(Ok(HashSet::new()), |vcf_positions, (i, record)| {
            debug!("Accesing {} record in the VCF...", i + 1);
            let record = record?;
            let reference = record.alleles()[0];
            let mut vcf_positions = vcf_positions?;
            debug!(
                "Succesfully obtained description {} with reference allele {}",
                record.desc(),
                String::from_utf8(reference.to_vec())?
            );
            if is_any_iupac_s(reference, 0, reference.len()) {
                debug!(
                    "Record {} in the VCF will be processed for potential 8-oxoG",
                    record.desc()
                );
                vcf_positions.insert(record.desc());
            }
            Ok(vcf_positions)
        })
}

/// Checks if nucleotide at specifc positions is IUPAC S
pub fn is_any_iupac_s(seq: &[u8], start: usize, end: usize) -> bool {
    if start >= seq.len() || end > seq.len() {
        error!("Reference sequence is shorter than BAM alignment positions");
        std::process::exit(1)
    } else {
        seq[start..end]
            .iter()
            .any(|nuc| matches!(nuc, b'G' | b'C' | b'g' | b'c'))
    }
}

/// Checks if a ref and alt allele are G --> T or A --> C
pub fn is_any_g2t_or_c2a(ref_allele: &[u8], alt_allele: &[u8]) -> bool {
    ref_allele
        .iter()
        .zip(alt_allele.iter())
        .any(|alleles| matches!(alleles, (b'G', b'T') | (b'C', b'A')))
}
