use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
pub(crate) enum RustyNuc {
    Reference {
        #[structopt(
            short,
            long,
            help = "Forward reads",
            required = true,
            parse(from_os_str)
        )]
        forward: Vec<PathBuf>,
        #[structopt(short, long, help = "Reverse reads", parse(from_os_str))]
        reverse: Vec<PathBuf>,
        #[structopt(help = "Reference", parse(from_os_str))]
        reference: Vec<PathBuf>,
    },
    Variant {
        #[structopt(
            short,
            long,
            help = "VCF/BCF containing calls based on the supplied bam file",
            parse(from_os_str)
        )]
        vcf: Option<PathBuf>,
        #[structopt(
            help = "BAM files to process for possible location of oxo-g effects",
            parse(from_os_str)
        )]
        bam: Vec<PathBuf>,
    },
}
