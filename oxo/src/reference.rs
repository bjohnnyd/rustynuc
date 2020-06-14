use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;

// TODO: Maybe better us FMDIndex and concatenate fmindex for forward and reverse as
// https://docs.rs/bio/0.31.0/src/bio/data_structures/fmindex.rs.html#417
pub struct Reference {
    pub(crate) forward: FMIndex<BWT, Less, Occ>,
    pub(crate) reverse: FMIndex<BWT, Less, Occ>,
}

impl Reference {
    pub fn new<T: AsRef<str>>(text: T) -> Self {
        let mut complement = dna::revcomp(text.as_ref().as_bytes());
        complement.push(b'$');
        let text = format!("{}$", text.as_ref());
        let forward = Self::create_fm_index(text.as_bytes());
        let reverse = Self::create_fm_index(&complement);

        Self { forward, reverse }
    }

    fn create_fm_index(text: &[u8]) -> FMIndex<BWT, Less, Occ> {
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        FMIndex::new(bwt, less, occ)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_ref_creation() {
        let text = "GCCTTAACATTATTACGCCTA";
        let texts = "GCCTTAACATTATTACGCCTA$";
        let pattern = "TTT";
        let reference = Reference::new(&text);

        let mut revcomp = dna::revcomp(text.as_bytes());
        revcomp.push(b'$');
        let sa = suffix_array(texts.as_bytes());
        let suffix_indices = reference.forward.backward_search(pattern.as_bytes().iter());
        dbg!(&suffix_indices);
        let pos = suffix_indices.occ(&sa);
        println!(
            "Reference '{}' has pattern '{}' occuring at '{:?}'",
            &text, &pattern, pos
        );
    }
}
