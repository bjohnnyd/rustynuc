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
pub(crate) struct OxoKmerRead<'a> {
    fwd_read: &'a fastq::Record,
    rev_read: &'a fastq::Record,
    kmer_size: usize,
    fwd_kmers: HashMapFx<&'a [u8], Vec<u32>>,
    rev_kmers: HashMapFx<&'a [u8], Vec<u32>>,
}

#[derive(Debug)]
pub enum Orientation {
    FWD,
    REV,
}

impl<'a> OxoKmerRead<'a> {
    pub fn new(fwd_read: &'a fastq::Record, rev_read: &'a fastq::Record, kmer_size: usize) -> Self {
        let fwd_kmers = Self::get_oxo_kmers(&fwd_read, Orientation::FWD, kmer_size);
        let rev_kmers = Self::get_oxo_kmers(&rev_read, Orientation::REV, kmer_size);
        Self {
            fwd_read,
            fwd_kmers,
            rev_read,
            rev_kmers,
            kmer_size,
        }
    }

    fn get_oxo_kmers(
        read: &'a fastq::Record,
        orientation: Orientation,
        kmer_size: usize,
    ) -> HashMapFx<&'a [u8], Vec<u32>> {
        let oxo_nuc = match orientation {
            Orientation::FWD => b'G',
            Orientation::REV => b'C',
        };

        hash_kmers(read.seq(), kmer_size)
            .into_iter()
            .filter(|(kmer, _)| kmer.contains(&oxo_nuc))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmers() {
        let seq = b"GCATGGTGTA";
        let rev_seq = b"TGTACCCTTA";
        let qual = vec![19, 22, 15, 20, 19, 22, 15, 20];
        let fwd_fastq = fastq::Record::with_attrs("fwd_read", None, seq, &qual);
        let rev_fastq = fastq::Record::with_attrs("rev_read", None, rev_seq, &qual);

        let kmer_read = OxoKmerRead::new(&fwd_fastq, &rev_fastq, 3);
        let expected_fwd_first_kmer = b"GCA";

        assert_eq!(
            kmer_read.fwd_kmers.get(&expected_fwd_first_kmer[..]),
            Some(&vec![0])
        );
    }
}
