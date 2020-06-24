use crate::error::Error;
use rust_htslib::bam::pileup::Pileup;
use std::collections::HashMap;

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
/// Contains summary of a pileup split into
/// `first in pair` and `second in pair`.
/// Contains helper functions for determining likelihood of
/// reads containing oxo damage
#[allow(missing_docs)]
pub struct OxoPileup {
    pub ref_id: u32,
    pub ref_depth: u32,
    pub ref_pos: u32,
    // ref_base: Option<u8>,
    pub first: FirstReadCount,
    pub second: SecondReadCount,
}

#[derive(Debug)]
/// Contains nucleotide counts for this position
/// for reads labeled as `first in pair (0x64)`
pub struct FirstReadCount(pub HashMap<u8, u32>);
/// Contains nucleotide counts for this position
/// for reads labeled as `second in pair (0x128)`
#[derive(Debug)]
pub struct SecondReadCount(pub HashMap<u8, u32>);

// TODO: Need to decide if want to implement the reference base
impl OxoPileup {
    pub fn new(pileup: Pileup) -> Self {
        let ref_id = pileup.tid();
        let ref_depth = pileup.depth();
        let ref_pos = pileup.pos();

        let mut first = FirstReadCount(HashMap::new());
        let mut second = SecondReadCount(HashMap::new());

        for alignment in pileup.alignments() {
            if let Some(pos) = alignment.qpos() {
                let record = alignment.record();
                let base = record.seq()[pos];

                if record.is_first_in_template() {
                    let count = first.0.entry(base).or_insert(0);
                    *count += 1;
                } else if record.is_last_in_template() {
                    let count = second.0.entry(base).or_insert(0);
                    *count += 1;
                }
            }
        }

        Self {
            ref_id,
            ref_depth,
            ref_pos,
            first,
            second,
        }
    }
}
