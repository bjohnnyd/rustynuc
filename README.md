# rustynuc

[![Release][ico-version]][link-version]
[![Build Status][ico-travis]][link-travis]
[![Software License][ico-license]](./LICENSE.md)


Tool to calculate the likelihood of 8-oxoG damage based on alignment characteristics.

## Install

### Binary

Precompiled binaries are provided below:

| ![picture](static/64px-Tux.png) | ![picture](static/64px-MacOS_logo.png)  |
| :-----------------------------: | :-------------------------------------: |
| [TAR][linux-tar] | [TAR][osx-tar]  |
| [ZIP][linux-zip] | [ZIP][osx-zip]  |

### Cargo 

If you have cargo installed or have installed [RUSTUP](https://rustup.rs/), you can install directly from github
``` bash
$ cargo install --git https://github.com/bjohnnyd/rustynuc
```

### Build 

To compile from source rustup is required and can be obtained [HERE](https://rustup.rs/).  After installing rustup download the release archive file and build:

```bash
$ git clone https://github.com/bjohnnyd/rustynuc.git && cd rustynuc && cargo build --release 
```

All releases and associated binaries and archives are accessible under [Releases](https://github.com/bjohnnyd/rustynuc/releases).

## Usage

``` bash
$ ./rustynuc -h
```

```
rustynuc 0.2.0

USAGE:
    rustynuc [FLAGS] [OPTIONS] <bam>

FLAGS:
    -a, --all                Whether to just print results for all positions
    -h, --help               Prints help information
    -n, --no-qval            Skip calculating qvalue
    -p, --pseudocount        Whether to use pseudocounts (adds +1 to all counts) when calculating statistics
    -V, --version            Prints version information
    -v, --verbosity          Determines verbosity of the processing, can be specified multiple times -vvv
    -w, --with-track-line    Include track line (for correct visualization with IGV)

OPTIONS:
        --alpha <alpha>              FDR threshold [default: 0.2]
    -b, --bcf <bcf>                  BCF/VCF of variants called on the supplied alignment file
        --bed <bed>                  A BED file to restrict analysis to specific regions
        --fisher-sig <fisher-sig>    Significance threshold for Fisher's test [default: 0.05]
        --max-depth <max-depth>      Maximum pileup depth to use [default: 1000]
    -m, --min-reads <min-reads>      Minimum number of aligned reads in ff or fr orientation for a position to be
                                     considered [default: 4]
    -q, --quality <quality>          Minimum base quality to consider [default: 20]
    -r, --reference <reference>      Optional reference which will be used to determine sequence context and type of
                                     change
    -t, --threads <threads>          Number of threads [default: 4]

ARGS:
    <bam>    Alignments to correct for possible 8-oxoG damage

```

### Output

The default output (if not `--bcf/-b` is provided) is a BED file with the following info:

```
1. Chromosome
2. Start
3. End
4. Name (format is `<chromosome>_<start>_<end>` or if reference is provided `<chromosome>_<base>_<start>_<end>`
5. -log10 of p-value (p-value is the smallest of the A/C and G/T )
6. Strand
7. Depth
8. Adenine FF:FR counts
9. Cytosine FF:FR counts
10. Guanine FF:FR counts
11. Thymine FF:FR counts
12. A/C two-sided p-value Fisher's Exact Test
13. G/T two-sided p-value Fisher's Exact Test
(14). Sequnce Context (if reference provided)
14/15. adj. pvalue 
15/16. Significant at set FDR value (1 if yes, 0 if not)
```

To get only positions with p-value below 0.05: 

```bash
$ rustynuc -r tests/input/ref.fa.gz tests/alignments/oxog.bam | awk '$12 < 0.05 || $13 < 0.05'  | gzip > sig.bed.gz
```

Alternatively, if a VCF/BCF is provided the output is a VCF.  Multiple summaries are provided in the VCF file:

```
##FILTER=<ID=OxoG,Description="OxoG two-sided p-value < 0.05">,
##FILTER=<ID=OccurenceInsufficient,Description="There is not a sufficient number of reads aligning in the FF and FR orientation">,
##INFO=<ID=OXO_DEPTH,Number=1,Type=Integer,Description="OxoG Pileup Depth">,
##INFO=<ID=ADENINE_FF_FR,Number=2,Type=Integer,Description="Adenine counts at FF and FR">,
##INFO=<ID=CYTOSINE_FF_FR,Number=2,Type=Integer,Description="Cytosine counts at FF and FR">,
##INFO=<ID=GUANINE_FF_FR,Number=2,Type=Integer,Description="Guanine counts at FF and FR">,
##INFO=<ID=THYMINE_FF_FR,Number=2,Type=Integer,Description="Thymine counts at FF and FR">,
##INFO=<ID=A_C_PVAL,Number=1,Type=Float,Description="A/C two-sided p-value">,
##INFO=<ID=G_T_PVAL,Number=1,Type=Float,Description="G/T two-sided p-value">,
##INFO=<ID=AF_FF_FR,Number=A,Type=Float,Description="Alternate frequency calculations on the FF and FR">,
##INFO=<ID=OXO_CONTEXT,Number=1,Type=String,Description="3mer Context of reference">,
```

`AF_FF_FR` can be used to filter based on AF on the `FF` or `FR` orientations.  For each alternate allele, there are two AF provided so for example to filter the first alternate positions `AF_FF_FR[0]` and `AF_FF_FR[1]` can be used.  The command below will filter using the AF on FF/FR and also `FILTER=="PASS"` ensures only position with `p-val < 0.05` are returned.

```bash
$ FILTERCMD='TYPE =="snp" && AF > 0.04 && FILTER=="PASS" && (AF_FF_FR=="." | (AF_FF_FR[0] >= 0.04 && AF_FF_FR[1] >= 0.04))'
$ rustynuc --pseudocounts -r tests/input/ref.fa.gz --b tests/input/oxog.vcf.gz tests/alignments/oxog.bam | bcftools filter -Oz -i "$FILTERCMD" > nonoxog.vcf.gz
```

## Authors

- [Johnny Debebe][link-author]

## License

The MIT License (MIT). Please see [License File](LICENSE.md) for more information.


## Notes

- all test results are heavily dependent on depth so in cases where depth is not high the `AF_FF_FR` filter might be a better option
- fisher's exact is affected heavily by 0 counts so `pseudocounts` are going to produce less conservative results
- fisher's exact test is fairly conservative so might underestimate the true numbers
- FDR will be heavily dependent on %GC of the genome, size of the genome, whether a reference was provided, a VCF is provided or the test was restricted to specific regions.  All these factors will affect how many tests will need to be performed.  
- for all the filtering methods `DEPTH` will heavily determine power

[ico-version]: https://img.shields.io/github/v/release/bjohnnyd/rustynuc?include_prereleases&style=flat-square
[ico-license]: https://img.shields.io/github/license/bjohnnyd/rustynuc?color=purple&style=flat-square
[ico-travis]: https://img.shields.io/travis/com/bjohnnyd/rustynuc?style=flat-square
[ico-downloads]: https://img.shields.io/packagist/dt/:vendor/rustynuc.svg?style=flat-square

[link-version]: https://github.com/bjohnnyd/rustynuc/releases/latest
[link-travis]: https://travis-ci.com/bjohnnyd/rustynuc
[link-downloads]: https://packagist.org/packages/bjohnnyd/rustynuc
[link-author]: https://github.com/bjohnnyd

[linux-tar]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-unknown-linux-gnu.tar.gz
[linux-zip]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-unknown-linux-gnu.zip
[osx-tar]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-apple-darwin.tar.gz
[osx-zip]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-apple-darwin.zip
[windows-tar]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-pc-windows-gnu.tar.gz
[windows-zip]: https://github.com/bjohnnyd/rustynuc/releases/latest/download/x86_64-pc-windows-gnu.zip
