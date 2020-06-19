use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not read FASTA record")]
    FastaRecordError(#[from] std::io::Error),
    #[error("Could not convert bytes as it is invalid UTF-8")]
    NotUTF8(#[from] std::string::FromUtf8Error),
}
