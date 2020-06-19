use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not read FASTA record")]
    FastaRecordError(#[from] std::io::Error),
}
