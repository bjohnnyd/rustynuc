#![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

//! Crate for intepreting various nucleotide damage related to oxidation (e.g. `oxo-G`).
//! `oxo-G` damage results generates technical artefacts that cause `G > T` errors first
//! read and `C -> A` erros in the second read.  In addition, there is an enrichment in
//! specific sequnce contexts as `CCG -> CAG`.
//!

mod kmer;
mod reference;

pub(crate) const DEFAULT_KMER_SIZE: usize = 5;
pub(crate) const DEFAULT_KMER_MIN: u32 = 10_000;
pub(crate) const DEFAULT_SAMPLING_RATE_OCC: u32 = 3;

//TODO: To many combinations to produce even with 10mers
// however with clostridium exponential growth observed
// i.e. the `test.ref` 5 give to 1 unique, 6 to 243 and 7 to 2970 on forward only.
// 1. Could be also possible to iterate over kmer sizes until a significance has been reached? Maybe
// 10,000 for default or as a fraction of refrence or max number
// 2. Then maybe iterate these kmers to create an oxidized set from reverse and forward strand (7mer
// is unlikely to have > 2 mutations)
// 3. Check the oxidized set to not occur in the reference, check that the expected alt is not in
// the reference if not due to oxo damage
// 4. Create the informative kmer struct that will contain the oxidized kmer seq, the differnce and
//    location and type of change to original, the expected alt on reverse
// Might be best to hash unique reference kmers and then hash forward and reverse kmers
// Iterate over the forward and reverse kmers
//
//
