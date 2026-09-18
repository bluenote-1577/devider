v0.1.0 - 9-18-2026

- Fixed a bug with overwriting output files when specifying many contigs to phase (e.g. with `--bed-file`), which could drop haplotypes for all but the last contig/region
- Fixed a crash (divide by zero) that could occur on some datasets when constructing haplotypes
- Fixed incorrect results when using k=1
- Fixed incorrect reporting of the number of SNPs in a contig
- Fixed a missing newline in the `ids.txt` output file


