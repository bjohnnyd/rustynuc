use crate::kmer::Orientation;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use std::collections::HashMap;

// TODO: Maybe better us FMDIndex and concatenate fmindex for forward and reverse as
// https://docs.rs/bio/0.31.0/src/bio/data_structures/fmindex.rs.html#417
pub struct Reference<'a> {
    pub(crate) fmindex: FMIndex<BWT, Less, Occ>,
    pub(crate) sa: RawSuffixArray,
    pub(crate) seq: &'a [u8],
    pub(crate) revcomp: Vec<u8>,
    pub(crate) ref_len: usize,
}

// TODO: in order to hash reverse comp need to create it might be better to own the key in hash
// then in some cases
impl<'a> Reference<'a> {
    pub fn new(seq: &'a [u8]) -> Self {
        let ref_len = seq.len();
        let revcomp = dna::revcomp(seq);
        let text_builder: Vec<&[u8]> = vec![seq, b"$", &revcomp, b"$"];
        let text = text_builder.concat();

        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&text);
        let bwt = bwt(&text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, crate::DEFAULT_SAMPLING_RATE_OCC, &alphabet);

        let fmindex = FMIndex::new(bwt, less, occ);

        Self {
            fmindex,
            sa,
            seq,
            revcomp,
            ref_len,
        }
    }

    pub fn search_pattern(&self, pattern: &[u8]) -> (Vec<usize>, Vec<usize>) {
        let sai = self.fmindex.backward_search(pattern.iter());
        let positions = sai.occ(&self.sa);
        positions.iter().partition(|pos| **pos < self.ref_len)
    }

    pub fn hash_kmers(&self, kmer_size: usize) -> HashMap<&[u8], usize> {
        let mut hashes = hash_kmers(self.seq, kmer_size)
            .into_iter()
            .map(|(seq, pos)| (seq, pos.len()))
            .collect::<HashMap<&[u8], usize>>();

        hash_kmers(&self.revcomp[..], kmer_size)
            .into_iter()
            .for_each(|(seq, pos)| {
                let count = hashes.entry(seq).or_insert(0);
                *count += pos.len();
            });

        hashes
    }
}

pub struct OxidizedReads {
    original_read: String,
    oxidized_reads: Vec<String>,
    read_diffs: Vec<usize>,
}

impl OxidizedReads {
    pub fn new<T: AsRef<str>>(text: T, orientation: Orientation) -> Self {
        let (ref_oxo_nuc, alt_oxo_nuc) = match orientation {
            Orientation::FWD => ('G', 'T'),
            Orientation::REV => ('C', 'A'),
        };

        let mut oxidized_reads = Self::oxidize(text, ref_oxo_nuc, alt_oxo_nuc);
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
    fn oxidize<T: AsRef<str>>(text: T, ref_oxo_nuc: char, alt_oxo_nuc: char) -> Vec<String> {
        let text = text.as_ref();

        if let Some(i) = text.find(ref_oxo_nuc) {
            let prefix = &text[0..i];
            let suffix = &text[(i + 1)..];

            let comb = OxidizedReads::oxidize(&suffix, ref_oxo_nuc, alt_oxo_nuc);
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
    fn test_ref_creation() {
        let text = b"GCCTTAACATTATTACGCCTA";
        let reference = Reference::new(text);

        let pattern = b"GGC";

        let (fwd_matches, mut rev_matches) = reference.search_pattern(pattern);

        rev_matches.remove(0);

        assert!(fwd_matches.is_empty());
        assert_eq!(rev_matches.len(), 1);
    }

    #[test]
    fn test_ref_real() {
        let pattern = b"GCTAAGGT";

        let (rdr, _) = niffler::from_path("tests/refs/ref.fa.gz").unwrap();
        let reference = bio::io::fasta::Reader::new(rdr);

        for record in reference.records() {
            let record = record.unwrap();
            let record_reference = Reference::new(record.seq());
            let (fwd_matches, rev_matches) = record_reference.search_pattern(pattern);

            assert_eq!(fwd_matches.len(), 2);
            assert_eq!(rev_matches.len(), 1);

            // NOTE: Need to check that the kmer occurs only once in reference
            // to ensure that the closest oxo sequence is most likely source from this kmer
            for (seq, occurence) in
                record_reference
                    .hash_kmers(8)
                    .iter()
                    .filter(|(seq, occurences)| {
                        (seq.iter().any(|c| *c == b'C' || *c == b'G')) && **occurences == 1usize
                    })
            {
                let seq = unsafe { String::from_utf8_unchecked(seq.to_vec()) };

                println!("Kmer: '{}' occurs {} times.", seq, occurence)
            }
        }
    }

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
}
