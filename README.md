# rustynuc

[![Release][ico-version]][link-version]
[![Build Status][ico-travis]][link-travis]
[![Software License][ico-license]](LICENSE.md)


Tool to calculate the likelihood of 8-oxoG damage based on alignment characteristics.

## Install

### Binary

The simplest way to install is using the precompiled binaries provided.

### Build 

To compile from source rustup is required and can be obtained [HERE](https://rustup.rs/).  After installing rustup download the release archive file and build:

```bash
$ git clone https://github.com/bjohnnyd/rustynuc.git && cd rustynuc && cargo build --release 
```

The compiled binary can then be ran using:

``` bash
$ ./target/release/rustynuc -h
```

All releases and associated binaries and archives are accessible under [Releases](https://github.com/bjohnnyd/rustynuc/releases).

## Usage

``` bash
$ ./rustynuc -h
```

```
rustynuc 0.1.0

USAGE:
    rustynuc [FLAGS] [OPTIONS] <bam>

FLAGS:
    -a, --all                Whether to just print results for all positions
    -h, --help               Prints help information
    -p, --pseudocount        Whether to use pseudocounts (adds +1 to all counts) when calculating statistics
    -V, --version            Prints version information
    -v, --verbosity          Determines verbosity of the processing, can be specified multiple times -vvv
    -w, --with-track-line    Include track line (for correct visualization with IGV)

OPTIONS:
        --alpha <alpha>              FDR threshold [default: 0.2]
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
4. Name (format is `chromosome_start_end` or if reference is provided `chromosome_base_start_end`
5. -log10 of p-value (p-value is the smalles of the A/C and G/T )
6. Strand
7. Depth
8. Adenine FF:FR counts
9. Cytosine FF:FR counts
10. Guanine FF:FR counts
11. Thymine FF:FR counts
12. A/C two-sided p-value Fisher's Exact Test
13. G/T two-sided p-value Fisher's Exact Test
14. adj. pvalue 
15. Significant at set FDR value (1 if yes, 0 if not)
```

## Authors and Citation

- [Bisrat Johnathan Debebe][link-author]

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

