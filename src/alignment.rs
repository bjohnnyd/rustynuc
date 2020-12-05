use crate::error::Error;
use fishers_exact::fishers_exact;
use log::{debug, error};
use rust_htslib::bam::pileup::Pileup;
use std::collections::HashMap;

type Result<T> = std::result::Result<T, Error>;

/// Contains summary of a pileup split into `first in pair` and `second in pair`.
/// Contains helper functions for determining likelihood of reads containing oxo damage.
#[allow(missing_docs)]
pub struct OxoPileup {
    pub chrom: String,
    pub depth: u32,
    pub pos: u32,
    pub context: Option<[u8; 3]>,
    pub ff_count: [u32; 256],
    pub fr_count: [u32; 256],
    pub min_count: u32,
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
        min_count: u32,
        min_qual: u8,
        pseudocount: bool,
        seq: Option<T>,
        no_overlapping: bool,
    ) -> Self {
        let mut depth = 0;
        let pos = pileup.pos();
        let pseudocount = if pseudocount { 1 } else { 0 };
        let context = match seq {
            Some(seq) => Some(get_seq_context(seq.as_ref(), pos)),
            None => None,
        };

        let mut ff_count = [0; 256];
        let mut fr_count = [0; 256];

        for alignment in pileup.alignments() {
            if let Some(pos) = alignment.qpos() {
                let record = alignment.record();
                let base = record.seq()[pos];
                let qual = record.qual()[pos];

                if record.is_first_in_template() && qual > min_qual {
                    if record.is_reverse() {
                        fr_count[base as usize] += 1;
                    } else {
                        ff_count[base as usize] += 1;
                    }
                }

                if no_overlapping {
                    if record.insert_size() <= 0 && !record.is_first_in_template() {
                        continue;
                    } else {
                        depth += 1;
                    }
                }
            }
        }

        if !no_overlapping {
            depth = pileup.depth()
        }

        Self {
            chrom,
            depth,
            pos,
            context,
            ff_count,
            fr_count,
            min_count,
            pseudocount,
        }
    }

    /// Checks if specific nucleotide is at sufficient numbers `[ff, fr]`
    pub fn nuc_sufficient(&self, nuc: u8) -> [bool; 2] {
        debug!("Counting nucleotide {}", nuc as char);
        let nuc = nuc.to_ascii_uppercase() as usize;
        let ff_count = self.ff_count[nuc] + self.ff_count[nuc + 32];
        let fr_count = self.fr_count[nuc] + self.fr_count[nuc + 32];
        debug!("{} in ff and {} in fr", ff_count, fr_count);
        [ff_count >= self.min_count, fr_count >= self.min_count]
    }

    /// Returns the counts of A, C, G and T
    pub fn nuc_counts(&self, orientation: Orientation) -> Vec<u32> {
        crate::NUCLEOTIDES
            .iter()
            .map(|nuc| match orientation {
                Orientation::FF => {
                    self.ff_count[*nuc as usize]
                        + self.ff_count[(*nuc + 32) as usize]
                        + self.pseudocount
                }
                Orientation::FR => {
                    self.fr_count[*nuc as usize]
                        + self.fr_count[(*nuc + 32) as usize]
                        + self.pseudocount
                }
            })
            .collect::<Vec<u32>>()
    }

    /// Returns a description of the pileup in the form `chrom:pos`
    pub fn desc(&self) -> String {
        format!("{}:{}", self.chrom, self.pos)
    }

    /// Checks if the coverage is above the supplied minimum
    pub fn occurence_sufficient(&self) -> bool {
        self.nuc_counts(Orientation::FF)
            .into_iter()
            .zip(self.nuc_counts(Orientation::FR).into_iter())
            .filter(|(ff_count, fr_count)| {
                *ff_count >= self.min_count || *fr_count >= self.min_count
            })
            .count()
            > 1
    }

    /// Checks if there is more than one nucleotide appearing at 0 or more counts
    pub fn is_monomorphic(&self) -> bool {
        let ff_monomorphic = self
            .nuc_counts(Orientation::FF)
            .into_iter()
            .filter(|count| *count > 0)
            .count()
            <= 1;

        if ff_monomorphic {
            self.nuc_counts(Orientation::FR)
                .into_iter()
                .filter(|count| *count > 0)
                .count()
                <= 1
        } else {
            false
        }
    }

    /// Returns the two-tailed p.value for the A/C or G/T comparison
    pub fn get_nuc_pval(&self, nuc: u8) -> Result<f64> {
        let idx = match nuc {
            b'A' | b'C' | b'c' | b'a' => 0,
            b'G' | b'T' | b'g' | b't' => 2,
            _ => return Err(Error::IncorrectNuc(nuc.to_string())),
        };

        let ff_nuc_counts = self.nuc_counts(Orientation::FF);
        let fr_nuc_counts = self.nuc_counts(Orientation::FR);

        let counts = [
            ff_nuc_counts[idx],
            fr_nuc_counts[idx],
            ff_nuc_counts[idx + 1],
            fr_nuc_counts[idx + 1],
        ];

        let fisher_result = fishers_exact(&counts).unwrap();
        Ok(fisher_result.two_tail_pvalue)
    }

    /// Returns if context sequence is known the p-value based on the reference nucleotide
    /// or if the reference is unknown the smaller of the A/C and G/T two-tailed p.values
    pub fn pval(&self) -> Result<f64> {
        match self.context {
            Some(context) if is_pos_iupac_s(&context[..], 1) => self.get_nuc_pval(context[1]),
            _ => {
                let ac = self.get_nuc_pval(b'A')?;
                let gt = self.get_nuc_pval(b'G')?;

                if ac < gt {
                    Ok(ac)
                } else {
                    Ok(gt)
                }
            }
        }
    }

    /// Crates a string representation of the OxoPileup summary
    fn to_string(&self) -> Result<String> {
        let summary = self
            .nuc_counts(Orientation::FF)
            .iter()
            .zip(self.nuc_counts(Orientation::FR))
            .fold(
                format!("{}", self.depth),
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
            self.get_nuc_pval(b'A')?,
            self.get_nuc_pval(b'G')?
        ))
    }

    /// Convertes the pileups summary into a bed entry
    pub fn to_bed_row(&self) -> Result<String> {
        let mut seq = String::from("");

        let name = match self.context {
            Some(context) => {
                seq = format!("\t{}", String::from_utf8(context.to_vec())?);

                format!(
                    "{}_{}_{}_{}",
                    &self.chrom,
                    context[1] as char,
                    self.pos,
                    self.pos + 1
                )
            }
            _ => format!("{}_{}_{}", &self.chrom, self.pos, self.pos + 1),
        };

        Ok(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            &self.chrom,
            self.pos,
            self.pos + 1,
            name,
            -(self.pval()?.log10()),
            "*",
            self.to_string()?,
            seq,
        ))
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum Orientation {
    FR,
    FF,
}

/// Returns the 3-mer sequence context with the pileup in the middle. For pileup at zero
/// the first base is labeled `X` and for pileup at end the last
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

/// Checks if nucleotide at specifc positions is IUPAC S
pub fn is_pos_iupac_s(seq: &[u8], pos: usize) -> bool {
    if pos >= seq.len() {
        error!("Reference sequence is shorter than BAM alignment positions. Got pileup at position {} but reference is only {} long", pos + 1, seq.len());
        std::process::exit(1)
    } else {
        match seq[pos] {
            b'G' | b'C' | b'g' | b'c' => true,
            _ => false,
        }
    }
}
