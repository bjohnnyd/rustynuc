use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Read};

// TODO: Steps:
// 1. Use pileup
// 2. Get position
// 3. From position get reference nuclotide, if G do the following
// 3. For each alignment at position
// 4. Get strand, forward/reverse, nucleotide count
// 5. Create Struct to store this information
// NOTE: In graph representation there is a nice way
// to traverse the graph and in each case check that the order of nucleotides by occurence
// can be represented as GCTA on forward and CGAT on reverse check if this is true  cross the graph
// where it is not true check that the diff ones are representative of oxo damage
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_alignment_file() {
        let home = std::env::var("HOME").unwrap();
        let f = format!("{}/data/rustynuc/oxo.bam", home);
        let mut alignment = bam::Reader::from_path(f).unwrap();

        for record in alignment.records() {
            let mut record = record.unwrap();
            if record.is_proper_pair() && !record.is_supplementary() && !record.is_duplicate() {
                let read_seq = record.seq();
                let ref_pos = record.reference_positions();
                dbg!(&read_seq.encoded[0..2]);
                dbg!(&ref_pos[0..2]);
                println!(
                    "Read '{}' is on '{}' and reverse is '{}' ",
                    record.tid(),
                    record.strand(),
                    record.is_reverse()
                )
            }
        }
    }
}
