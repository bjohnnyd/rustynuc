use thiserror::Error;

#[derive(Debug, Error)]
/// Lists error commonly caused by I/O
pub enum Error {
    #[error("Could not read FASTA record")]
    /// Could not read an entry in a fasta file
    FastaRecordError(#[from] std::io::Error),
    #[error("Could not convert bytes as it is invalid UTF-8")]
    /// Data is not in UTF-8 format
    NotUTF8(#[from] std::string::FromUtf8Error),
    #[error("Could not process the BAM file")]
    /// Bam Reading Error
    BamError(#[from] rust_htslib::bam::Error),
    #[error("Count too large for statistical testing")]
    /// Fisher's exact error
    FishersError(#[from] fishers_exact::TooLargeValueError),
}
