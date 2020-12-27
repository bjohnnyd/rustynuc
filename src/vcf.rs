use crate::NUCLEOTIDES;
use crate::{
    alignment::{Orientation, OxoPileup},
    error::Error,
};
use crate::{genomic::is_any_g2t_or_c2a, Result};
use bcf::Record;
use log::debug;
use rust_htslib::bcf;
use std::collections::HashMap;

/// Header information updated in VCF
pub const HEADER_RECORDS: [&[u8]; 12] = [
    br#"##FILTER=<ID=FishersOxoG,Description="OxoG Fisher's exact p-value < 0.05">"#,
    br#"##FILTER=<ID=InsufficientCount,Description="Insufficient number of reads aligning in the FF or FR orientation for calculations">"#,
    br#"##FILTER=<ID=AfTooLow,Description="AF is below 0.04 on either FF or FR orientation">"#,
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

/// Name for Fisher's exact filter
pub const FISHERS_OXO_FILTER: &[u8] = b"FishersOxoG";
/// Name for insufficient counts filter
pub const INSUFFICIENT_FILTER: &[u8] = b"InsufficientCount";
/// Name for AF too low on FF/FR filter
pub const AF_LOW_FILTER: &[u8] = b"AfTooLow";

/// Minimum AF on both FF and FR orientation in order for a variant to pass AfTooLow
pub const AF_MIN: f32 = 0.02;

/// Returns the maximum frequency across FR and FF and returns it. Updates the supplied VCF record
/// with all calculated metadata from the provided oxo pileup.
pub fn update_vcf_record(record: &mut bcf::Record, oxo: &OxoPileup) -> Result<Option<f32>> {
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

    record.push_info_float(b"AC_PVAL", &[ac_pval as f32])?;
    record.push_info_float(b"GT_PVAL", &[gt_pval as f32])?;

    if let Some(ref seq) = oxo.context {
        record.push_info_string(b"OXO_CONTEXT", &[seq])?;
    }

    let mut g2t_or_c2a_max_freq: Option<f32> = None;

    match record.alleles()[0] {
        ref_allele if ref_allele == b"G" || ref_allele == b"C" => {
            let oxo_af =
                record
                    .alleles()
                    .into_iter()
                    .skip(1)
                    .fold(Vec::new(), |mut oxo_af, alt_allele| {
                        let alt_counts = counts.get(&alt_allele[0]).unwrap_or_else(|| &[0, 0]);
                        let ref_counts = counts.get(&ref_allele[0]).unwrap_or_else(|| &[0, 0]);
                        let ff_af =
                            alt_counts[0] as f32 / (ref_counts[0] as f32 + alt_counts[0] as f32);
                        let fr_af =
                            alt_counts[1] as f32 / (ref_counts[1] as f32 + alt_counts[1] as f32);
                        oxo_af.push(ff_af);
                        oxo_af.push(fr_af);

                        match (ref_allele, alt_allele) {
                            (b"G", b"T") | (b"C", b"A") => {
                                let max_orientation_af = ff_af.max(fr_af);
                                g2t_or_c2a_max_freq =
                                    Some(g2t_or_c2a_max_freq.map_or(max_orientation_af, |af| {
                                        af.max(max_orientation_af)
                                    }));
                            }
                            _ => {}
                        }

                        oxo_af
                    });

            record.push_info_float(b"FF_FR_AF", oxo_af.as_slice())?;
        }
        _ => {}
    }
    Ok(g2t_or_c2a_max_freq)
}

/// Applies the `FishersOxoG` filter to the VCF record based on the provided significance threshold and
/// frequency above which both ff and fr will be labelled as PASS
pub fn apply_fishers_filter(
    record: &mut bcf::Record,
    sig_threshold: f32,
    af_both_pass: f32,
) -> Result<()> {
    let ac_pval = record.get_float_value("AC_PVAL")?[0].to_owned();
    let gt_pval = record.get_float_value("GT_PVAL")?[0].to_owned();
    let above_both_pass = record
        .get_float_value("FF_FR_AF")?
        .chunks(2)
        .any(|freqs| freqs.iter().all(|freq| freq > &af_both_pass));

    if (ac_pval < sig_threshold || gt_pval < sig_threshold) & !above_both_pass {
        debug!("Record {} is labeled as OxoG", record.desc());
        let hview = record.header();
        let id = hview.name_to_id(FISHERS_OXO_FILTER)?;
        record.push_filter(id);
    } else if above_both_pass {
        debug!(
            "Record {} passes the Fisher's OxoG filter as it is above the AF threshold set at {}",
            record.desc(),
            af_both_pass
        );
    } else {
        debug!(
            "Record {} passes the Fisher's OxoG filter as the pValue is above the significance threshold set at {}",
            record.desc(),
            sig_threshold
        );
    }
    Ok(())
}

/// Applies the `AfTooLow` filter to VCF record by checking if any FF or FR frequencies across the
/// ALT alleles is below the provided `min_af`
pub fn apply_af_too_low_filter(record: &mut bcf::Record, min_af: f32) -> Result<()> {
    debug!(
        "Checking if record {} passes the AfTooLow Filter...",
        record.desc()
    );
    let ref_allele = record.alleles()[0].to_owned();
    let ff_fr_freqs = record.get_float_value("FF_FR_AF")?.to_vec();

    let af_too_low = record
        .alleles()
        .iter()
        .skip(1)
        .zip(ff_fr_freqs.chunks(2))
        // TODO: NEED TO DEAL WITH DINUCLEOTIDE CHANGES (not alt_allele[0])
        .any(|(alt_allele, freqs)| {
            if is_any_g2t_or_c2a(&ref_allele[..], alt_allele) {
                freqs.into_iter().any(|freq| freq < &min_af)
            } else {
                false
            }
        });
    if af_too_low {
        debug!("Record {} has an AF below the min of {} in at least one read orientation and will be labeled as AfTooLow", record.desc(), min_af);
        let hview = record.header();
        let id = hview.name_to_id(AF_LOW_FILTER)?;
        record.push_filter(id);
    } else {
        debug!("Record {} passes the AfTooLow filter", record.desc());
    }

    Ok(())
}

trait InfoExt {
    fn get_float_value<'a>(&'a mut self, tag: &'a str) -> Result<&'a [f32]>;
}

impl InfoExt for Record {
    fn get_float_value<'a>(&'a mut self, tag: &'a str) -> Result<&'a [f32]> {
        self.info(tag.as_bytes())
            .float()?
            .ok_or_else(|| Error::NoInfoTag(tag.to_string()))
    }
}
#[cfg(test)]
#[allow(unused)]
mod tests {
    use super::InfoExt;
    use rust_htslib::{bcf, bcf::Read};

    pub fn get_record_at(path: &str, idx: usize) -> bcf::Record {
        let mut vcf = bcf::Reader::from_path(path).unwrap();
        vcf.records().skip(idx).next().unwrap().unwrap()
    }
    fn test_get_frequencies_for_g2t_or_c2a() {
        let mut record = get_record_at(
            "/home/johnny/projects/rust_dev/rustynuc/tests/input/high_af/high_af.oxog.vcf.gz",
            0,
        );
        let freqs = record.get_float_value("FF_FR_AF").unwrap();
    }
}
