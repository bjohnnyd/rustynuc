use crate::alignment::{Orientation, OxoPileup};
use crate::Result;
use crate::NUCLEOTIDES;
use log::debug;
use rust_htslib::{bcf, bcf::header::Id};
use std::collections::HashMap;

/// Header information added to VCF
pub const HEADER_RECORDS: [&[u8]; 11] = [
    br#"##FILTER=<ID=OxoG,Description="OxoG two-sided p-value < 0.05">"#,
    br#"##FILTER=<ID=InsufficientCount,Description="Insufficient number of reads aligning in the FF or FR orientation">"#,
    br#"##INFO=<ID=OXO_DEPTH,Number=1,Type=Integer,Description="OxoG Pileup Depth">"#,
    br#"##INFO=<ID=ADENINE_FF_FR,Number=2,Type=Integer,Description="Adenine counts in FF and FR orientations">"#,
    br#"##INFO=<ID=CYTOSINE_FF_FR,Number=2,Type=Integer,Description="Cytosine counts in FF and FR orientations">"#,
    br#"##INFO=<ID=GUANINE_FF_FR,Number=2,Type=Integer,Description="Guanine counts in FF and FR orientations">"#,
    br#"##INFO=<ID=THYMINE_FF_FR,Number=2,Type=Integer,Description="Thymine counts in FF and FR orientations">"#,
    br#"##INFO=<ID=AC_PVAL,Number=1,Type=Float,Description="A/C two-sided p-value">"#,
    br#"##INFO=<ID=GT_PVAL,Number=1,Type=Float,Description="G/T two-sided p-value">"#,
    br#"##INFO=<ID=FF_FR_AF,Number=A,Type=Float,Description="Alternate frequency calculations on the FF and FR">"#,
    br#"##INFO=<ID=OXO_CONTEXT,Number=1,Type=String,Description="3mer reference sequence context">"#,
];

/// Name for significant filter
pub const OXO_FILTER: &[u8] = b"OxoG";
/// Name for insufficient counts filter
pub const INSUFFICIENT_FILTER: &[u8] = b"InsufficientCount";

pub fn update_vcf_record(
    record: &mut bcf::Record,
    oxo: &OxoPileup,
    filter_id: (Id, Id),
    sig_threshold: f64,
    ff_fr_threshold: f32,
    oxo_ceiling: f32,
) -> Result<()> {
    let counts = oxo
        .nuc_counts(Orientation::FF)
        .iter()
        .zip(oxo.nuc_counts(Orientation::FR))
        .enumerate()
        .map(|(i, (ff_count, fr_count))| {
            let counts = [
                (ff_count - oxo.pseudocount) as i32,
                (fr_count - oxo.pseudocount) as i32,
            ];
            (NUCLEOTIDES[i], counts)
        })
        .collect::<HashMap<u8, [i32; 2]>>();

    record.push_info_integer(b"OXO_DEPTH", &[oxo.depth as i32])?;
    record.push_info_integer(
        b"ADENINE_FF_FR",
        counts.get(&b'A').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"CYTOSINE_FF_FR",
        counts.get(&b'C').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"GUANINE_FF_FR",
        counts.get(&b'G').unwrap_or_else(|| &[0, 0]),
    )?;
    record.push_info_integer(
        b"THYMINE_FF_FR",
        counts.get(&b'T').unwrap_or_else(|| &[0, 0]),
    )?;

    let ac_pval = oxo.get_nuc_pval(b'C')?;
    let gt_pval = oxo.get_nuc_pval(b'G')?;
    let is_oxog = ac_pval < sig_threshold || gt_pval < sig_threshold;
    record.push_info_float(b"AC_PVAL", &[ac_pval as f32])?;
    record.push_info_float(b"GT_PVAL", &[gt_pval as f32])?;

    if let Some(ref seq) = oxo.context {
        record.push_info_string(b"OXO_CONTEXT", &[seq])?;
    }

    let mut insufficient_count = false;
    let mut above_ff_fr_threshold = false;

    match record.alleles()[0] {
        ref_allele if ref_allele == b"G" || ref_allele == b"C" => {
            let oxo_af =
                record
                    .alleles()
                    .into_iter()
                    .skip(1)
                    .fold(Vec::new(), |mut oxo_af, alt_allele| {
                        let alt_counts = counts.get(&alt_allele[0]).unwrap_or_else(|| &[0, 0]);
                        let ref_counts = counts.get(&ref_allele[0]).unwrap_or_else(|| &[1, 1]);
                        let ff_af =
                            alt_counts[0] as f32 / (ref_counts[0] as f32 + alt_counts[0] as f32);
                        let fr_af =
                            alt_counts[1] as f32 / (ref_counts[1] as f32 + alt_counts[1] as f32);
                        above_ff_fr_threshold = (ff_af > oxo_ceiling || fr_af > oxo_ceiling)
                            || (ff_af > ff_fr_threshold && fr_af > ff_fr_threshold);
                        match (ref_allele, alt_allele) {
                            (b"G", b"T") | (b"C", b"A") if !oxo.occurence_sufficient() => {
                                insufficient_count = true
                            }
                            _ => {}
                        }
                        oxo_af.push(ff_af);
                        oxo_af.push(fr_af);
                        oxo_af
                    });

            if is_oxog && !above_ff_fr_threshold {
                debug!(
                    "Record {} is labeled as OxoG and is below the ff_fr_threshold of {}",
                    record.desc(),
                    ff_fr_threshold
                );
                record.push_filter(filter_id.0)
            }

            if insufficient_count {
                debug!("Record {} has insufficient count", record.desc());
                record.push_filter(filter_id.1)
            };
            record.push_info_float(b"FF_FR_AF", oxo_af.as_slice())?;
        }
        _ => {}
    }
    Ok(())
}
