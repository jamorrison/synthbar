# synthbar

## tl;dr

`synthbar` adds a synthetic cell barcode to each read in a provided FASTQ.

## Install and Build

```
git clone git@github.com:jamorrison/synthbar.git
cd synthbar
make
```

## Overview

Most single cell RNA-seq protocols include a cell barcode at the beginning of each read to distinguish which cell the
read came from. This barcode is followed by a unique molecular index (UMI), which allows for tracking individual
RNA fragments and correct deduplication of RNA sequencing data. Finally, a short linker sequence is required to connect
the cell barcode and UMI to the cDNA insert that will be sequenced.

While the majority of single cell protocols include both the cell barcode and the UMI, some protocols (particularly
plate based single cell protocols and even some bulk protocols) utilize a UMI without a cell barcode (see Note 1 below).
For historical reasons, most tools for analyzing RNA-seq data with UMIs require a cell barcode before the UMI and can't
handle reads without barcodes.

To overcome this hurdle, `synthbar` prepends a 7-base cell barcode (`CATATAC`) to the sequence string of each read in
the FASTQ, as well as a matching 7-base quality score (`IIIIIII`) in the quality string. If the user would rather use a
different barcode, the `--barcode` option allows for a user-defined barcode. In this case, a string of `I`'s of length
equal to the requested barcode will be prepended to the quality string. If desired, the linker sequence can be removed
as well, leaving just the barcode, UMI, and cDNA sequence.

_Note 1: With respect to the plate based single cell protocols, individual cells are sorted into each well in a plate
and then an index is added to distinguish cells for demultiplexing after sequencing (similar to distinguishing samples
on your standard Illumina sequencer). In the case of bulk protocols, a "cell barcode" by itself doesn't make sense as
there are many cells within the sample. However, even in bulk protocols, UMIs can serve a scientific purpose. In order
to properly account for the UMIs in these protocols, they should be processed in a similar manner to the various single
cell protocols with UMIs._

## Usage

```
Usage: synthbar [options] <FASTQ with UMIs>

Output options:
    -o, --output STR           name of output file [stdout]
Processing Options:
    -b, --barcode STR          barcode to prepend to each read [CATATAC]
    -U, --umi-first            add barcode to read after the UMI [off]
    -r, --remove-linker        remove linker from read [not removed]
    -l, --linker-length INT    length of linker to remove [6]
    -u, --umi-length INT       length of UMI before linker [8]
    -h, --help                 print usage and exit
        --version              print version and exit

Note 1: Input FASTQ can be gzip compressed or uncompressed
```

|       Option        |     Input      | Description                                                               |
|:--------------------|:---------------|:--------------------------------------------------------------------------|
| -o, --output        | string         | name of output file (defaults to stdout), does not write gzip'd files     |
| -b, --barcode       | string         | barcode to add instead of CATATAC (does not check if composed of ATCG's)  |
| -U, --umi-first     | -              | place the barcode after the UMI in the new read                           |
| -r, --remove-linker | -              | remove linker sequence from read (not removed by default)                 |
| -l, --linker-length | integer (>= 0) | length of linker to remove (default is 6), not used if `-r` not provided  |
| -u, --umi-length    | integer (>= 0) | length of UMI before linker (default is 8), not used if `-r` not provided |
| -h, --help          | -              | print usage and exit                                                      |
| --version           | -              | print version and exit                                                    |

Note, for protocols with no linking sequence, it is suggested to ignore the linker-related options, as this will ensure
everything is written after the UMI and eliminate the potential for inadvertently removing cDNA sequence.

## Read Structure

| In / Out | Linker? | UMI First? | Remove Linker? | Structure                                       |
|:--------:|:-------:|:----------:|:--------------:|:------------------------------------------------|
| Input    | Yes     | NA         | NA             | `( UMI ) + ( LINKER ) + ( cDNA )`               |
| Input    | No      | NA         | NA             | `( UMI ) + ( cDNA )`                            |
| Output   | Yes     | No         | No             | `( BARCODE ) + ( UMI ) + ( LINKER ) + ( cDNA )` |
| Output   | Yes     | No         | Yes            | `( BARCODE ) + ( UMI ) + ( cDNA )`              |
| Output   | Yes     | Yes        | No             | `( UMI ) + ( BARCODE ) + ( LINKER ) + ( cDNA )` |
| Output   | Yes     | Yes        | Yes            | `( UMI ) + ( BARCODE ) + ( cDNA )`              |
| Output   | No      | No         | No             | `( BARCODE ) + ( UMI ) + ( cDNA )`              |
| Output   | No      | Yes        | No             | `( UMI ) + ( BARCODE ) + ( cDNA )`              |

## Acknowledgments

  - `synthbar` uses two utilities from `klib`, `kseq` (for handling FASTQs) and `kstring` (for writing updated reads).

## Citation

At this time, there is no plan to publish `synthbar`. However, if you would kindly cite the GitHub page whenever you
use `synthbar` in your analysis, the reference would be greatly appreciated.
