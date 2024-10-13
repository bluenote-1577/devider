# devider - long-read haplotypes from mixtures of "small" sequences

**devider** is a method that separates long reads (Nanopore or PacBio) of a mixture of sequences into groups with similar alleles. This is called "phasing" or "haplotyping". 

devider is a "local haplotyping" method, so it works best when the sequence-of-interest is approximately the size of the reads. For bacterial genome-scale haplotyping, consider another tool such as [floria](https://github.com/bluenote-1577/floria).

### Example use cases:

* mixed viral long-read samples (e.g. co-infections or quasispecies)
* amplicon/enriched sequencing of specific genes
* haplotyping small sections of multi-strain bacterial communities

<p align="center">
  <img width="600" height="200" src="https://github.com/user-attachments/assets/c0a82bb5-7feb-4d13-ab59-04da2bce52b3", caption="asdf">
</p>
<p align="center">
  <i>
High-depth, heterogeneous sequencing that spans a 1kb gene.
  </i>
</p>

<p align="center">
  <img width="600" height="200" src="https://github.com/user-attachments/assets/34cb8bcf-8f23-47e4-b2f6-8515a21d3cf4", caption="asdf">
</p>
<p align="center">
  <i>
Separated groups ("haplotypes") after running devider.
  </i>
</p>

### Why devider?

Similar tools exist for detection of similar haplotypes in mixtures. devider was developed to fill the following gaps:

* **Speed and low-memory** - devider scales approximately linearly with sequencing depth and # of SNPs. > 30,000x coverage genes can be haplotyped in minutes. 
* **High heterogeneity and coverage** - devider uses a de Bruijn Graph approach, which works with very diverse samples (> 10 haplotypes)
* **Ease-of-use + interpretable outputs** - conda installable, engineered in rust, simple command line. Outputs are easy to interpret (haplotagged BAM or MSA). 

## Install

> [!NOTE]
> As of 2024-10-12, conda install is not ready yet. Will be available in the next few days. 

```sh
mamba install -c bioconda devider
devider -h 
```

See the [installation instructions on the wiki](https://github.com/bluenote-1577/devider/wiki/Installation) if you want to compile directly or want a static binary.

## Quick Start after install 

### Option 1 (more flexible): Running devider with VCF + BAM
```sh
git clone https://github.com/bluenote-1577/devider
cd devider
devider -b hiv_test/3000_95_3.bam  -v hiv_test/3000_95_3.vcf.gz  -r hiv_test/OR483991.1.fasta

# results folder
ls devider_output
```
### Option 2 (easier): Running devider with reads 
```sh
run_devider_pipeline -i reads.fq.gz -r reference.fa -o pipeline_output

# results folder
ls pipeline_output

# intermediate files (bam + vcf files)
ls pipeline_output/pipeline_files
```

> [!NOTE]
>  If you **did not** install via conda, do the following instead. 
>```sh
>mamba install -c bioconda tabix samtools lofreq minimap2
>git clone https://github.com/bluenote-1577/devider
>./devider/scripts/run_devider_pipeline -i reads.fq.gz -r reference.fa -o pipeline_output
>```

## Manuals, tutorials, and cookbook

### How to use devider

* [Cookbook](https://github.com/bluenote-1577/devider/wiki/Cookbook) - see here for usage examples.
* [Advanced usage manual](https://github.com/bluenote-1577/devider/wiki/Advanced-usage-manual) - see here for more detailed information about parameters and usage.
* [Output format](https://github.com/bluenote-1577/devider/wiki/Output-format) - for more information on how to interpret outputs.

### Tutorials

* [Tutorial 1 - getting started with devider](https://github.com/bluenote-1577/devider/wiki/Tutorial-1:-getting-started-with-devider)

## Citation

Forthcoming.
