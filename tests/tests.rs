use assert_cmd::prelude::*;
use predicates::str::{contains, is_empty, is_match, PredicateStrExt};
use std::process::Command;
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
