use crate::error::Error;
use rust_htslib::bam::pileup::Pileup;
use std::collections::HashMap;

type Result<T> = std::result::Result<T, Error>;

// TODO: Need to test if just doing a fishers exact test on:
// 1  A+C ff v A+C fr
// 2. G+T ff v G+T fr
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
    pub first_on_forward: NucleotideCount,
    pub second_on_forward: NucleotideCount,
}

#[derive(Debug)]
/// Contains nucleotide counts for a specific position
pub struct NucleotideCount(pub HashMap<u8, u32>);

// TODO: Need to decide if want to implement the reference base
impl OxoPileup {
    /// Creates a pileup summary in terms of `F1R2` and `F2R1`
    pub fn new(pileup: Pileup) -> Self {
        let ref_id = pileup.tid();
        let ref_depth = pileup.depth();
        let ref_pos = pileup.pos();

        let mut first_on_forward = NucleotideCount(HashMap::new());
        let mut second_on_forward = NucleotideCount(HashMap::new());

        for alignment in pileup.alignments() {
            if let Some(pos) = alignment.qpos() {
                let record = alignment.record();
                let base = record.seq()[pos];

                if record.is_first_in_template() {
                    if record.is_reverse() {
                        let count = second_on_forward.0.entry(base).or_insert(0);
                        *count += 1;
                    } else {
                        let count = first_on_forward.0.entry(base).or_insert(0);
                        *count += 1;
                    }
                }
            }
        }

        Self {
            ref_id,
            ref_depth,
            ref_pos,
            first_on_forward,
            second_on_forward,
        }
    }
}

impl std::fmt::Display for OxoPileup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut summary = format!("{}\t{}", self.ref_pos, self.ref_depth);
        for nuc in crate::NUCLEOTIDES.iter() {
            let ff_count = self.first_on_forward.0.get(nuc).unwrap_or(&0_u32);
            let fr_count = self.second_on_forward.0.get(nuc).unwrap_or(&0_u32);

            summary.push_str(&format!("\t{}:{}", ff_count, fr_count));
        }
        write!(f, "{}", summary)
    }
}
