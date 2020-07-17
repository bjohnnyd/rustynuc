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
    -b, --bed <bed>                  A BED file to restrict analysis to specific regions
        --fisher-sig <fisher-sig>    Significance threshold for Fisher's test [default: 0.05]
        --max-depth <max-depth>      Maximum pileup depth to use [default: 400]
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

Output is a BED file with the following info:

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


## Authors and Citation

- [Johnny Debebe][link-author]

## License

The MIT License (MIT). Please see [License File](LICENSE.md) for more information.

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
