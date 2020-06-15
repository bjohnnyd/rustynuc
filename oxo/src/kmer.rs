use bio::alignment::sparse::{hash_kmers, HashMapFx};
use bio::io::fastq;

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
    kmers: HashMapFx<&'a [u8], Vec<u32>>,
    read_type: Orientation,
}

#[derive(Debug)]
pub enum Orientation {
    FWD,
    REV,
}

impl<'a> KmerRead<'a> {
    fn new(read: &'a fastq::Record, kmer_size: usize, read_type: Orientation) -> Self {
        let kmers = hash_kmers(read.seq(), kmer_size);

        Self {
            read,
            kmer_size,
            kmers,
            read_type,
        }
    }
}

pub struct OxidizedReads {
    original_read: String,
    oxidized_reads: Vec<String>,
    read_diffs: Vec<usize>,
}

impl OxidizedReads {
    pub fn new<T: AsRef<str>>(text: T, orientation: Orientation) -> Self {
        let mut oxidized_reads = Self::oxidize(text, &orientation);
        let original_read = oxidized_reads.remove(0);
        let read_diffs = oxidized_reads
            .iter()
            .map(|oxo_read| {
                oxo_read
                    .chars()
                    .zip(original_read.chars())
                    .filter(|(oxo_nuc, orig_nuc)| oxo_nuc != orig_nuc)
                    .count()
            })
            .collect();

        Self {
            original_read,
            oxidized_reads,
            read_diffs,
        }
    }

    /// "Oxidizes" the specified string depending on orientation.
    /// FWD produces all combinations of G --> C mutations
    /// while REV produces all C-->T mutations for an input kmer
    fn oxidize<T: AsRef<str>>(text: T, orientation: &Orientation) -> Vec<String> {
        let text = text.as_ref();
        let (ref_oxo_nuc, alt_oxo_nuc) = match orientation {
            Orientation::FWD => ('G', 'T'),
            _ => ('C', 'A'),
        };

        if let Some(i) = text.find(ref_oxo_nuc) {
            let prefix = &text[0..i];
            let suffix = &text[(i + 1)..];

            let comb = OxidizedReads::oxidize(&suffix, &orientation);
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oxidize() {
        let kmer = "GCATTGAGTT";

        let fwd_oxo = OxidizedReads::new(kmer, Orientation::FWD);
        let rev_oxo = OxidizedReads::new(kmer, Orientation::REV);

        assert_eq!(fwd_oxo.original_read, rev_oxo.original_read);

        assert_eq!(fwd_oxo.oxidized_reads.len(), 7);
        assert_eq!(rev_oxo.oxidized_reads.len(), 1);

        assert_eq!(fwd_oxo.read_diffs.iter().max(), Some(&3usize));
        assert_eq!(rev_oxo.read_diffs.iter().max(), Some(&1usize));
    }

    #[test]
    fn test_kmers() {
        let seq = b"GCATGGFG";
        let qual = vec![19, 22, 15, 20, 19, 22, 15, 20];
        let fastq = fastq::Record::with_attrs("test_read", None, seq, &qual);

        let kmer_read = KmerRead::new(&fastq, 3, Orientation::FWD);

        let expected_first_kmer = b"GCA";
        assert_eq!(
            kmer_read.kmers.get(&expected_first_kmer[..]),
            Some(&vec![0])
        );
    }
}
