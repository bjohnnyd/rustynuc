use bio::io::fastq;
use std::collections::HashMap;

// Basic idea could be that take kmers across fwd reads and filter for ones containing `G` in
// forward reads and ones containing `C` in reverse reads. On forward reads there should be
// decrease in G but not a decrease in C while on reverse there should be decrease in C but not a
// decrease in G for reverse complement.  Is this possible without reference present, example with
// foward would be take kmers with single T and substitute with G reverse complement both and count
// reverse complement of G should not exceed T similar while reverse complement of.
// Need to index the reads for this to work!!

#[derive(Debug)]
pub(crate) struct KmerRead<'a> {
    read: &'a fastq::Record,
    kmer_size: usize,
    motif_counts: HashMap<String, usize>,
    read_type: Orientation,
}

#[derive(Debug)]
pub enum Orientation {
    FWD,
    REV,
}

impl<'a> KmerRead<'a> {
    fn new(read: &'a fastq::Record, kmer_size: usize, read_type: Orientation) -> Self {
        let motif_counts = read
            .seq()
            .chunks(kmer_size)
            .filter(|chunk| {
                chunk.iter().any(|c| *c == b'G' || *c == b'T')
                    && chunk.len() == kmer_size
                    && !chunk.iter().all(|c| *c == b'G')
            })
            .fold(HashMap::new(), |mut motif_counts, chunk| {
                let motif = String::from_utf8(chunk.to_vec()).unwrap();
                let count = motif_counts.entry(motif).or_default();
                *count += 1;
                motif_counts
            });

        Self {
            read,
            kmer_size,
            motif_counts,
            read_type,
        }
    }
}

// TODO: Still returns original could change by:
// 1. Passing an additional variable original_text
// 2. Wrapping in another function that will filter out original
// 3. Filter from result
// NOTE: Might be best to create struct for this so I can also keep track of original kmer
// as might be needed to get diff to original in order to create the expected variant on reverse
// complement too if this one is
/// "Oxidizes" the specified string depending on orientation.
/// FWD produces all combinations of G --> C mutations
/// while REV produces all C-->T mutations for an input kmer
pub fn oxidize<T: AsRef<str>>(text: T, orientation: &Orientation) -> Vec<String> {
    let text = text.as_ref();
    let (ref_oxo_nuc, alt_oxo_nuc) = match orientation {
        Orientation::FWD => ('G', 'T'),
        _ => ('C', 'A'),
    };

    if let Some(i) = text.find(ref_oxo_nuc) {
        let prefix = &text[0..i];
        let suffix = &text[(i + 1)..];

        let comb = oxidize(&suffix, &orientation);
        let mods = vec![ref_oxo_nuc, alt_oxo_nuc];
        return mods.iter().fold(Vec::new(), |mut oxidized_set, nuc| {
            for suffix_result in &comb {
                let seq = format!("{}{}{}", &prefix, nuc, &suffix_result);
                oxidized_set.push(seq);
            }
            oxidized_set
        });
    }

    return vec![text.to_string()];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oxidize() {
        let orientation = Orientation::FWD;
        let kmer = "GCATTGAGTT";

        let mut fwd_oxo = oxidize(kmer, &Orientation::FWD);
        let mut rev_oxo = oxidize(kmer, &Orientation::REV);

        fwd_oxo.remove(0);
        rev_oxo.remove(0);

        println!(
            "For Kmer '{}' forward oxidized versions are {:?}",
            kmer, &fwd_oxo
        );

        println!(
            "For Kmer '{}' reverse oxidized versions are {:?}",
            kmer, &rev_oxo
        );

        assert_eq!(fwd_oxo.len(), 7);
        assert_eq!(rev_oxo.len(), 1);
    }
}
