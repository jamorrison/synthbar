# synthbar

Add synthetic cell barcode to FASTQ reads

## Usage

```
Usage: synthbar [options] <FASTQ with UMIs>

Output options:
    -o, --output STR           name of output file [stdout]
Processing Options:
    -r, --remove-linker        remove linker from read [not removed]
    -l, --linker-length INT    length of linker to remove [6]
    -u, --umi-length INT       length of UMI before linker [8]
    -h, --help                 print usage and exit
        --version              print version and exit

Note 1: Input FASTQ can be gzip compressed or uncompressed
```

|       Option        |     Input      | Description |
|:--------------------|:---------------|:------------|
| -o, --output        | string         | name of output file (defaults to stdout), does not write gzip'd files     |
| -r, --remove-linker | -              | remove linker sequence from read (not removed by default)                 |
| -l, --linker-length | integer (>= 0) | length of linker to remove (default is 6), not used if `-r` not provided  |
| -u, --umi-length    | integer (>= 0) | length of UMI before linker (default is 8), not used if `-r` not provided |
| -h, --help          | -              | print usage and exit                                                      |
| --version           | -              | print version and exit                                                    |
