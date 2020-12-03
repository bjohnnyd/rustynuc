use crate::Result;
use bio::utils::Interval;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;

pub fn create_pileups<'a>(
    bam: &'a mut bam::Reader,
    threads: usize,
    max_depth: u32,
) -> Result<bam::pileup::Pileups<'a, bam::Reader>> {
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
