use assert_cmd::prelude::*;
use lazy_static::lazy_static;
use predicates::str::{contains, is_empty, is_match, PredicateStrExt};
use regex::Regex;
use std::process::Command;

lazy_static! {
    static ref FF_FR_AF_REGEX: Regex = Regex::new("FF_FR_AF=(.+?),(.+?)\t").unwrap();
}

fn extract_filter(output: &[u8]) -> String {
    String::from_utf8(output.to_vec())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("#"))
        .flat_map(|l| l.split("\t").skip(6).next())
        .collect()
}

#[test]
fn cli_no_args() {
    Command::cargo_bin("rustynuc").unwrap().assert().failure();
}
#[test]
fn cli_no_such_file() {
    Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&["tests/no_such_file.bam"])
        .assert()
        .failure()
        .stderr(contains("FileNotFound").trim());
}

#[test]
fn cli_normal_bam() {
    Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&["-r", "tests/input/ref.fa.gz", "tests/input/normal.bam"])
        .assert()
        .stdout(is_empty());
}

#[test]
fn cli_oxog_bam() {
    Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&["-r", "tests/input/ref.fa.gz", "tests/input/oxog.bam"])
        .assert()
        .stdout(contains("NC_009617.1_C_3963_3964"));
}

#[test]
fn cli_oxog_bcf() {
    Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "-r",
            "tests/input/ref.fa.gz",
            "-b",
            "tests/input/oxog.vcf.gz",
            "tests/input/oxog.bam",
        ])
        .assert()
        .stdout(contains("OxoG").count(5))
        .stdout(contains("PASS"));
}

#[test]
fn cli_print_all() {
    Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "--all",
            "-r",
            "tests/input/ref.fa.gz",
            "tests/input/oxog.bam",
        ])
        .assert()
        .stdout(is_match(r#"NC_009617.1\s1284\s1285"#).unwrap());
}

#[test]
fn cli_af_not_above_one() {
    let output = Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "--all",
            "-r",
            "tests/input/high_af/high_af.oxog.ref.fa.gz",
            "-b",
            "tests/input/high_af/high_af.oxog.vcf.gz",
            "tests/input/high_af/high_af.oxog.bam",
        ])
        .unwrap()
        .stdout;
    if let Some(captures) = FF_FR_AF_REGEX.captures(&String::from_utf8(output).unwrap()) {
        assert!(captures.len() == 3);
        assert!(captures.get(1).unwrap().as_str().parse::<f32>().unwrap() < 1.0);
        assert!(captures.get(2).unwrap().as_str().parse::<f32>().unwrap() < 1.0);
    } else {
        panic!("No FF_FR_AF in output")
    }
}

#[test]
fn cli_not_oxog_as_above_ff_fr_threshold() {
    let output = Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "-r",
            "tests/input/high_af/high_af.oxog.ref.fa.gz",
            "-b",
            "tests/input/high_af/high_af.oxog.vcf.gz",
            "tests/input/high_af/high_af.oxog.bam",
        ])
        .unwrap()
        .stdout;

    let filter_result = extract_filter(&output);
    assert_eq!(filter_result, "PASS".to_string())
}

#[test]
fn cli_oxog_false_positive_below_ff_fr_threshold_but_above_ceiling() {
    let output = Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "--fishers-af",
            "0.6",
            "-r",
            "tests/input/high_af/high_af.oxog.ref.fa.gz",
            "-b",
            "tests/input/high_af/high_af.oxog.vcf.gz",
            "tests/input/high_af/high_af.oxog.bam",
        ])
        .unwrap()
        .stdout;

    let filter_result = extract_filter(&output);

    assert_eq!(filter_result, "PASS".to_string());
}

#[test]
fn cli_oxog_false_positive_below_ff_fr_threshold_and_below_ceiling() {
    let output = Command::cargo_bin("rustynuc")
        .unwrap()
        .args(&[
            "--fishers-af",
            "0.6",
            "--oxo-af-ceiling",
            "1",
            "-r",
            "tests/input/high_af/high_af.oxog.ref.fa.gz",
            "-b",
            "tests/input/high_af/high_af.oxog.vcf.gz",
            "tests/input/high_af/high_af.oxog.bam",
        ])
        .unwrap()
        .stdout;

    let filter_result = extract_filter(&output);

    assert_eq!(filter_result, "OxoG".to_string());
}
