use crate::error::Error;
use fishers_exact::fishers_exact;
use rust_htslib::bam::pileup::Pileup;
use std::collections::HashMap;

type Result<T> = std::result::Result<T, Error>;

// TODO: Need to test if just doing a fishers exact test on:
// 1  A+C ff v A+C fr
// 2. G+T ff v G+T fr
// 3. First Pass get positional summaries
// 4. Identify issue positions and forward or reverse
// 5. Second pass collect fragment locations
// 6. Either return fragment IDs
// or if fastq provided remove reads from those fragments
#[derive(Debug, Eq, PartialEq)]
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
    pub ff_count: NucleotideCount,
    pub fr_count: NucleotideCount,
    pub sf_count: NucleotideCount,
    pub sr_count: NucleotideCount,
    pub min_count: Option<u32>,
}

#[derive(Debug, Eq, PartialEq)]
/// Contains nucleotide counts for a specific position
pub struct NucleotideCount(pub HashMap<u8, u32>);

// TODO: Need to decide if want to implement the reference base
// TODO: Needs refactoring kind of repetitive
impl OxoPileup {
    /// Creates a pileup summary in terms of `F1R2` and `F2R1`
    pub fn new(pileup: Pileup, min_count: Option<u32>, min_qual: u8) -> Self {
        let ref_id = pileup.tid();
        let ref_depth = pileup.depth();
        let ref_pos = pileup.pos();

        let mut ff_count = NucleotideCount(HashMap::new());
        let mut fr_count = NucleotideCount(HashMap::new());

        let mut sf_count = NucleotideCount(HashMap::new());
        let mut sr_count = NucleotideCount(HashMap::new());

        for alignment in pileup.alignments() {
            if let Some(pos) = alignment.qpos() {
                let record = alignment.record();
                let base = record.seq()[pos];
                let qual = record.qual()[pos];

                if record.is_first_in_template() && qual > min_qual {
                    if record.is_reverse() {
                        let count = fr_count.0.entry(base).or_insert(0);
                        *count += 1;
                    } else {
                        let count = ff_count.0.entry(base).or_insert(0);
                        *count += 1;
                    }
                }

                if record.is_last_in_template() && qual > min_qual {
                    if record.is_reverse() {
                        let count = sr_count.0.entry(base).or_insert(0);
                        *count += 1;
                    } else {
                        let count = sf_count.0.entry(base).or_insert(0);
                        *count += 1;
                    }
                }
            }
        }

        Self {
            ref_id,
            ref_depth,
            ref_pos,
            ff_count,
            fr_count,
            sf_count,
            sr_count,
            min_count,
        }
    }

    /// Checks if there are significant differences between nucleotides
    /// on ff v fr
    // TODO: This should return p-value and then needs to be corrected by FDR(order first and then just threshold) or Bonferroni
    // (stricter)
    pub fn is_imbalanced(&self, sig: f64) -> Result<bool> {
        let ff_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.ff_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let fr_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.fr_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let ac_counts = [
            ff_nuc_counts[0],
            ff_nuc_counts[1],
            fr_nuc_counts[0],
            fr_nuc_counts[1],
        ];

        let gt_counts = [
            ff_nuc_counts[2],
            ff_nuc_counts[3],
            fr_nuc_counts[2],
            fr_nuc_counts[3],
        ];

        let p_ac = fishers_exact(&ac_counts)?;
        let p_gt = fishers_exact(&gt_counts)?;

        let (ac_sufficient, gt_sufficient) = match self.min_count {
            // NOTE: I think only count for alt is necessary (i.e A or T)
            Some(count) => {
                let ac_sufficient = (ac_counts[0] > count || ac_counts[2] > count)
                    && (ac_counts[1] > count || ac_counts[3] > count);

                let gt_sufficient = (gt_counts[0] > count || gt_counts[2] > count)
                    && (gt_counts[1] > count || gt_counts[3] > count);
                (ac_sufficient, gt_sufficient)
            }
            _ => (true, true),
        };

        Ok((p_ac.two_tail_pvalue < sig && ac_sufficient)
            || (p_gt.two_tail_pvalue < sig && gt_sufficient))
    }

    /// Returns the two-tailed p.value for the A/C comparison
    pub fn p_ac(&self) -> f64 {
        let ff_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.ff_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let fr_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.fr_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let ac_counts = [
            ff_nuc_counts[0],
            fr_nuc_counts[0],
            ff_nuc_counts[1],
            fr_nuc_counts[1],
        ];

        let p_ac = fishers_exact(&ac_counts).unwrap();
        p_ac.two_tail_pvalue
    }

    /// Returns the two-tailed p.value for the G/T comparison
    pub fn p_gt(&self) -> f64 {
        let ff_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.ff_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let fr_nuc_counts = crate::NUCLEOTIDES
            .iter()
            .map(|nuc| *self.fr_count.0.get(nuc).unwrap_or(&0))
            .collect::<Vec<u32>>();

        let gt_counts = [
            ff_nuc_counts[2],
            fr_nuc_counts[2],
            ff_nuc_counts[3],
            fr_nuc_counts[3],
        ];

        let p_gt = fishers_exact(&gt_counts).unwrap();
        p_gt.two_tail_pvalue
    }

    /// Returns the smaller of the A/C and G/T two-tailed p.values
    pub fn min_p_value(&self) -> f64 {
        let ac = self.p_ac();
        let gt = self.p_gt();

        if ac < gt {
            ac
        } else {
            gt
        }
    }
}

impl std::fmt::Display for OxoPileup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut summary = format!("{}\t{}", self.ref_pos, self.ref_depth);
        for nuc in crate::NUCLEOTIDES.iter() {
            let ff_count = self.ff_count.0.get(nuc).unwrap_or(&0_u32);
            let fr_count = self.fr_count.0.get(nuc).unwrap_or(&0_u32);

            summary.push_str(&format!("\t{}:{}", ff_count, fr_count,));
        }
        write!(f, "{}\t{}\t{}", summary, self.p_ac(), self.p_gt())
    }
}
