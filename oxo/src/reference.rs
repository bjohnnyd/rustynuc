use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};

// TODO: Maybe better us FMDIndex and concatenate fmindex for forward and reverse as
// https://docs.rs/bio/0.31.0/src/bio/data_structures/fmindex.rs.html#417
pub struct Reference {
    pub(crate) fmindex: FMIndex<BWT, Less, Occ>,
    pub(crate) sa: RawSuffixArray,
    pub(crate) seq: Vec<u8>,
    pub(crate) ref_len: usize,
}

impl Reference {
    pub fn new<T: AsRef<str>>(text: T) -> Self {
        let text = text.as_ref().as_bytes();
        let ref_len = text.len();
        let revcomp = dna::revcomp(text);
        let text_builder: Vec<&[u8]> = vec![text, b"$", &revcomp, b"$"];
        let text = text_builder.concat();

        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&text);
        let bwt = bwt(&text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fmindex = FMIndex::new(bwt, less, occ);

        Self {
            fmindex,
            sa,
            seq: text,
            ref_len,
        }
    }
}

impl Reference {
    pub fn search_pattern(&self, pattern: &[u8]) -> (Vec<usize>, Vec<usize>) {
        let sai = self.fmindex.backward_search(pattern.iter());
        let positions = sai.occ(&self.sa);
        positions.iter().partition(|pos| **pos < self.ref_len)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_ref_creation() {
        let text = "GCCTTAACATTATTACGCCTA";
        let reference = Reference::new(&text);

        let pattern = b"GGC";

        let (fwd_matches, mut rev_matches) = reference.search_pattern(pattern);

        rev_matches.remove(0);

        assert!(fwd_matches.is_empty());
        assert_eq!(rev_matches.len(), 1);
    }
}
