use crate::error::Error;
use crate::kmer::Orientation;
use bio::alignment::sparse::hash_kmers;
use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::io::fasta;
use std::collections::HashMap;

// TODO: Maybe better us FMDIndex and concatenate fmindex for forward and reverse as
// https://docs.rs/bio/0.31.0/src/bio/data_structures/fmindex.rs.html#417
pub struct Reference {
    pub(crate) fmindex: FMIndex<BWT, Less, Occ>,
    pub(crate) sa: RawSuffixArray,
}

type Result<T> = std::result::Result<T, Error>;

// TODO: Need to create singletons using full sequences from record

impl Reference {
    pub fn new(seq: &[u8]) -> Self {
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

        Self { fmindex, sa }
    }

    pub fn from_reader<T: std::io::Read>(rdr: fasta::Reader<T>) -> Result<Self> {
        let mut text_builder = Vec::new();
        for record in rdr.records() {
            let record = record?;
            let seq = record.seq();
            let revcomp = dna::revcomp(seq.to_vec());
            text_builder.push(vec![seq, b"$", &revcomp[..], b"$"].concat());
        }
        let text = text_builder.concat();

        let (fmindex, sa) = Self::create_fm_index(&text[..]);

        Ok(Self { fmindex, sa })
    }

    fn create_fm_index(text: &[u8]) -> (FMIndex<BWT, Less, Occ>, RawSuffixArray) {
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&text);
        let bwt = bwt(&text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, crate::DEFAULT_SAMPLING_RATE_OCC, &alphabet);

        let fmindex = FMIndex::new(bwt, less, occ);
        (fmindex, sa)
    }

    pub fn search_pattern(&self, pattern: &[u8]) -> Vec<usize> {
        let sai = self.fmindex.backward_search(pattern.iter());
        sai.occ(&self.sa)
    }
}

#[derive(Debug)]
pub struct OxidizedReads {
    original_read: String,
    oxidized_reads: Vec<String>,
    read_diffs: Vec<usize>,
}

impl OxidizedReads {
    pub fn new(text: &[u8], orientation: Orientation) -> Result<Self> {
        let (ref_oxo_nuc, alt_oxo_nuc) = match orientation {
            Orientation::FWD => (b'G', b'T'),
            Orientation::REV => (b'C', b'A'),
        };

        let mut oxidized_reads = Self::oxidize(text, ref_oxo_nuc, alt_oxo_nuc)?;
        dbg!(&oxidized_reads);
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

        Ok(Self {
            original_read,
            oxidized_reads,
            read_diffs,
        })
    }

    /// "Oxidizes" the specified string depending on orientation.
    /// FWD produces all combinations of G --> C mutations
    /// while REV produces all C-->T mutations for an input kmer
    fn oxidize(text: &[u8], ref_oxo_nuc: u8, alt_oxo_nuc: u8) -> Result<Vec<String>> {
        let mut oxidized_set = Vec::new();
        if let Some(i) = text.iter().position(|c| *c == ref_oxo_nuc) {
            let prefix = &text[0..i];
            let suffix = &text[(i + 1)..];

            let comb = OxidizedReads::oxidize(&suffix, ref_oxo_nuc, alt_oxo_nuc)?;
            let mods = vec![ref_oxo_nuc, alt_oxo_nuc];
            let mut oxidized_set = Vec::new();
            for nuc in mods.iter() {
                for suffix_result in comb.iter().map(|seq| seq.as_bytes()) {
                    let nuc = vec![*nuc];
                    let seq = vec![prefix, &nuc[..], suffix_result];
                    oxidized_set.push(String::from_utf8(seq.concat())?);
                }
            }
            return Ok(oxidized_set);
        } else {
            oxidized_set.push(String::from_utf8(text.to_vec())?);
        }

        Ok(oxidized_set)
    }
}

/// Creates unique kmers occuring in the first and second strand of seqs
/// that contain nucleotides that can be affected by oxo damage `C` and `G`
pub fn create_singletons(seq: &[u8], kmer_size: usize) -> Result<Vec<String>> {
    Ok(hash_kmers_with_counts(seq, kmer_size)?
        .into_iter()
        .filter(|(seq, occurences)| {
            (seq.chars().any(|c| c == 'C' || c == 'G')) && *occurences == 1usize
        })
        .map(|(seq, occurences)| seq)
        .collect())
}

/// Creates the hash and counts for both strands of input sequence
pub fn hash_kmers_with_counts(seq: &[u8], kmer_size: usize) -> Result<HashMap<String, usize>> {
    let mut hashes = HashMap::new();
    let revcomp = dna::revcomp(seq.to_vec());

    for (kmer, pos) in hash_kmers(seq, kmer_size).into_iter() {
        hashes.insert(String::from_utf8(seq.to_vec())?, pos.len());
    }

    for (rev_kmer, pos) in hash_kmers(&revcomp[..], kmer_size).into_iter() {
        let count = hashes
            .entry(String::from_utf8(rev_kmer.to_vec())?)
            .or_insert(0);
        *count += pos.len();
    }

    Ok(hashes)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_ref_creation() {
        let text = b"GCCTTAACATTATTACGCCTA";
        let reference = Reference::new(text);
        let pattern = b"GGC";

        assert_eq!(reference.search_pattern(pattern).len(), 2);
    }

    #[test]
    fn test_ref_real() {
        let pattern = b"GCTAAGGT";

        let (rdr, _) = niffler::from_path("tests/refs/ref.fa.gz").unwrap();
        let reference = bio::io::fasta::Reader::new(rdr);

        let records = reference
            .records()
            .into_iter()
            .flat_map(|record| record.ok())
            .collect::<Vec<fasta::Record>>();

        let seqs = records
            .iter()
            .map(|record| record.seq())
            .collect::<Vec<&[u8]>>()
            .concat();

        #[allow(unused)]
        let mut singletons = Vec::new();

        // TODO: this will not work need to create a single reference object for the fasta
        // do i even need the struct, all i need is the singletons with certain n from the fasta
        // probably need to iter over fasta and have fm index of first$first_rev$second$second_rev
        for k in crate::DEFAULT_KMER_SIZE..21 {
            let record_reference = Reference::new(&seqs[..]);
            let raw_singletons = create_singletons(&seqs[..], k).unwrap();
            println!(
                "With kmer size {} found {} singletons",
                k,
                raw_singletons.len()
            );
            if raw_singletons.len() as u32 > crate::DEFAULT_KMER_MIN {
                singletons = raw_singletons;
                break;
            }
        }
    }

    #[test]
    fn test_oxidize() {
        let kmer = b"GCATTGAGTT";

        let fwd_oxo = OxidizedReads::new(kmer, Orientation::FWD).unwrap();
        let rev_oxo = OxidizedReads::new(kmer, Orientation::REV).unwrap();

        assert_eq!(fwd_oxo.original_read, rev_oxo.original_read);

        assert_eq!(fwd_oxo.oxidized_reads.len(), 7);
        assert_eq!(rev_oxo.oxidized_reads.len(), 1);

        assert_eq!(fwd_oxo.read_diffs.iter().max(), Some(&3usize));
        assert_eq!(rev_oxo.read_diffs.iter().max(), Some(&1usize));
    }
}
