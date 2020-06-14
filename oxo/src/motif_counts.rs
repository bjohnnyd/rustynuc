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

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_FQ_FIRST: &str = "tests/reads/normal/project_1_L001-ds.34e7f62f57824e8fb2e6f6f89738136a/culture-0-rep-2_S1_L001_R1_001.fastq.gz";
    static TEST_FQ_SECOND: &str = "tests/reads/normal/project_1_L001-ds.34e7f62f57824e8fb2e6f6f89738136a/culture-0-rep-2_S1_L001_R2_001.fastq.gz";

    static TEST_FQ_FIRST_OXO: &str = "tests/reads/high_oxo/project_1_L001-ds.64512d9e24c4410aa6a572eee0ed4b6b/culture-0-Rep-5_S10_L001_R1_001.fastq.gz";
    static TEST_FQ_SECOND_OXO: &str = "tests/reads/high_oxo/project_1_L001-ds.64512d9e24c4410aa6a572eee0ed4b6b/culture-0-Rep-5_S10_L001_R2_001.fastq.gz";
    #[test]
    fn test_fastq_complement() {
        use Orientation::*;
        let (fq1, _) =
            niffler::get_reader(Box::new(std::fs::File::open(TEST_FQ_FIRST_OXO).unwrap())).unwrap();
        let (fq2, _) =
            niffler::get_reader(Box::new(std::fs::File::open(TEST_FQ_SECOND_OXO).unwrap()))
                .unwrap();
        let rdr1 = fastq::Reader::new(fq1);
        let rdr2 = fastq::Reader::new(fq2);

        for (record1, record2) in rdr1.records().zip(rdr2.records()) {
            let record1 = record1.unwrap();
            let record2 = record2.unwrap();
            let fwd = KmerRead::new(&record1, 3, FWD);
            let rev = KmerRead::new(&record2, 3, FWD);

            dbg!(&fwd.motif_counts);
            dbg!(&rev.motif_counts);
        }
    }
}
