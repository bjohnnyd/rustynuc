use crate::error::Error;
use fishers_exact::fishers_exact;
use rust_htslib::bam::pileup::Pileup;
use std::collections::HashMap;

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Eq, PartialEq)]
/// Contains summary of a pileup split into
/// `first in pair` and `second in pair`.
/// Contains helper functions for determining likelihood of reads containing oxo damage
#[allow(missing_docs)]
pub struct OxoPileup {
    pub chrom: String,
    pub ref_depth: u32,
    pub ref_pos: u32,
    pub context: Option<[u8; 3]>,
    pub ff_count: NucleotideCount,
    pub fr_count: NucleotideCount,
    pub min_count: Option<u32>,
    pub pseudocount: u32,
}

#[derive(Debug, Eq, PartialEq)]
/// Contains nucleotide counts for a specific position
pub struct NucleotideCount(pub HashMap<u8, u32>);

impl OxoPileup {
    /// Creates a pileup summary in terms of `F1R2` and `F2R1`
    pub fn new<T: AsRef<[u8]>>(
        chrom: String,
        pileup: Pileup,
        min_count: Option<u32>,
        min_qual: u8,
        pseudocount: bool,
        seq: Option<T>,
    ) -> Self {
        let ref_depth = pileup.depth();
        let ref_pos = pileup.pos();
        let pseudocount = if pseudocount { 1 } else { 0 };
        let context = match seq {
            Some(seq) => Some(get_seq_context(seq.as_ref(), ref_pos)),
            None => None,
        };

        let mut ff_count = NucleotideCount(HashMap::new());
        let mut fr_count = NucleotideCount(HashMap::new());

        for alignment in pileup.alignments() {
            if let Some(pos) = alignment.qpos() {
                let record = alignment.record();
                let base = record.seq()[pos];
                let qual = record.qual()[pos];

                if record.is_first_in_template() && qual > min_qual {
                    if record.is_reverse() {
                        let count = fr_count.0.entry(base).or_insert(pseudocount);
                        *count += 1;
                    } else {
                        let count = ff_count.0.entry(base).or_insert(pseudocount);
                        *count += 1;
                    }
                }
            }
        }

        Self {
            chrom,
            ref_depth,
            ref_pos,
            context,
            ff_count,
            fr_count,
            min_count,
            pseudocount,
        }
    }

    /// Returns the counts of A, C, G and T
    pub fn nuc_counts(&self, read_type: ReadType) -> Vec<u32> {
        crate::NUCLEOTIDES
            .iter()
            .map(|nuc| match read_type {
                ReadType::FF => *self.ff_count.0.get(nuc).unwrap_or(&self.pseudocount),
                ReadType::FR => *self.fr_count.0.get(nuc).unwrap_or(&self.pseudocount),
            })
            .collect::<Vec<u32>>()
    }

    /// Checks if the base in the reference, if available, is G/C (IUPAC S)
    pub fn is_iupac_s(&self) -> Option<bool> {
        match self.context {
            Some(context) => match context[1] {
                b'G' | b'C' => Some(true),
                _ => Some(false),
            },
            None => None,
        }
    }

    /// Checks if the base in the reference, if available, is G/C (IUPAC S)
    pub fn occurence_sufficient(&self, min_count: u32) -> bool {
        self.nuc_counts(ReadType::FF)
            .into_iter()
            .zip(self.nuc_counts(ReadType::FR).into_iter())
            .filter(|(ff_count, fr_count)| *ff_count > min_count || *fr_count > min_count)
            .count()
            > 1
    }

    /// Checks if there are more than one nucleotides appearing at 0 or more counts
    pub fn is_monomorphic(&self) -> bool {
        let ff_monomorphic = self
            .nuc_counts(ReadType::FF)
            .into_iter()
            .filter(|count| *count > 0)
            .count()
            <= 1;

        if ff_monomorphic {
            self.nuc_counts(ReadType::FR)
                .into_iter()
                .filter(|count| *count > 0)
                .count()
                <= 1
        } else {
            false
        }
    }

    /// Checks if there are significant differences between nucleotides
    /// on ff v fr
    // TODO: This should return p-value and then needs to be corrected by FDR(order first and then just threshold) or Bonferroni
    // (stricter)
    pub fn is_imbalanced(&self, sig: f64) -> Result<bool> {
        let ff_nuc_counts = self.nuc_counts(ReadType::FF);
        let fr_nuc_counts = self.nuc_counts(ReadType::FR);

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

    /// Returns the two-tailed p.value for the A/C or G/T comparison
    pub fn get_pval(&self, nuc: u8) -> Result<f64> {
        let idx = match nuc {
            b'A' | b'C' => 0,
            b'G' | b'T' => 2,
            _ => return Err(Error::IncorrectNuc(nuc.to_string())),
        };

        let ff_nuc_counts = self.nuc_counts(ReadType::FF);
        let fr_nuc_counts = self.nuc_counts(ReadType::FR);

        let counts = [
            ff_nuc_counts[idx],
            fr_nuc_counts[idx],
            ff_nuc_counts[idx + 1],
            fr_nuc_counts[idx + 1],
        ];

        let fisher_result = fishers_exact(&counts).unwrap();
        Ok(fisher_result.two_tail_pvalue)
    }

    /// Returns the smaller of the A/C and G/T two-tailed p.values
    pub fn min_p_value(&self) -> Result<f64> {
        let ac = self.get_pval(b'A')?;
        let gt = self.get_pval(b'G')?;

        if ac < gt {
            Ok(ac)
        } else {
            Ok(gt)
        }
    }

    /// Crates a string representation of the OxoPileup summary
    fn to_string(&self) -> Result<String> {
        let summary = self
            .nuc_counts(ReadType::FF)
            .iter()
            .zip(self.nuc_counts(ReadType::FR))
            .fold(
                format!("{}", self.ref_depth),
                |mut summary, (ff_count, fr_count)| {
                    summary.push_str(&format!(
                        "\t{}:{}",
                        ff_count - self.pseudocount,
                        fr_count - self.pseudocount,
                    ));
                    summary
                },
            );

        Ok(format!(
            "{}\t{}\t{}",
            summary,
            self.get_pval(b'A')?,
            self.get_pval(b'G')?
        ))
    }

    pub fn to_bed_entry(&self) -> Result<String> {
        let mut seq = String::from("");

        let name = match self.context {
            Some(context) => {
                seq = format!("\t{}", String::from_utf8(context.to_vec())?);

                format!(
                    "{}_{}_{}_{}",
                    &self.chrom,
                    context[1] as char,
                    self.ref_pos,
                    self.ref_pos + 1
                )
            }
            _ => format!("{}_{}_{}", &self.chrom, self.ref_pos, self.ref_pos + 1),
        };

        Ok(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            &self.chrom,
            self.ref_pos,
            self.ref_pos + 1,
            name,
            -(self.min_p_value()?.log10()),
            "*",
            self.to_string()?,
            seq,
        ))
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum ReadType {
    FR,
    FF,
}

fn get_seq_context(seq: &[u8], pos: u32) -> [u8; 3] {
    let mut context = [0; 3];
    let idx = (seq.len() - 1) as u32;
    match pos {
        0 => {
            context[0] = b'X';
            context[1..3].copy_from_slice(&seq[0..2])
        }
        last_idx if last_idx == idx => {
            let last_idx = last_idx as usize;
            context[2] = b'X';
            context[0..2].copy_from_slice(&seq[(last_idx - 1)..(last_idx + 1)]);
        }
        i => context.copy_from_slice(&seq[(i - 1) as usize..(i + 2) as usize]),
    }

    context
}
