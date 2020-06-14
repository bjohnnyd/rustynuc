#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

//! Crate for intepreting various nucleotide damage related to oxidation (e.g. `oxo-G`).
//! `oxo-G` damage results generates technical artefacts that cause `G > T` errors first
//! read and `C -> A` erros in the second read.  In addition, there is an enrichment in
//! specific sequnce contexts as `CCG -> CAG`.
//!

mod motif_counts;
mod reference;

const DEFAULT_KMER_SIZE: usize = 21;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
