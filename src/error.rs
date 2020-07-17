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
    #[error("Could not create a reader for the FASTA file")]
    /// Read FASTA error
    FastaError(#[from] niffler::Error),
    #[error("Could not spawn threads")]
    /// Create thread pools erorr
    ThreadError,
    #[error("Only A, G, C and T nucleotides are allowed but got {0}")]
    /// Incorect nucleotide supplied
    IncorrectNuc(String),
    #[error("Could not read BED entry at line {0}")]
    /// Incorect nucleotide supplied
    BedRecordError(usize),
    #[error("Incorrect interval in BED entry on line {0}, end {1} is smaller than start {2} ")]
    /// Incorect nucleotide supplied
    IncorrectInterval(usize, u64, u64),
    #[error("Could not read VCF/BCF file")]
    /// Read VCF/BCF Error
    CouldNotReadVcf(#[from] rust_htslib::bcf::Error),
}
