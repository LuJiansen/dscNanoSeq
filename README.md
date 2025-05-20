# dscNanoSeq
Data preprocessing and analysis pipeline for dscNanoATAC &amp; dscNanoCUT

## Introduction
This repository contains scripts for preprocessing and analyzing data from dscNanoATAC and dscNanoCUT, droplet-based methods for single-cell epigenome profiling using Nanopore long-read seqeuencing.

## Install
This workflow is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://anaconda.org/) environment.

if you don't have mamba, install it first:
```
    conda update conda
    conda install mamba
```
create environmet for scNanoSeq
```
    mamba env create -f environment.yaml
```
Download UCSC Utils 
```
    rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedGraphToBigWig /path/to/envs/scNanoSeq/bin/
```

Then, install [ArchR](https://www.archrproject.com/) according to its tutorial. 
> [!CAUTION]
> ArchR (v1.0.1) should be used to create arrow files for each sequencing library due to its compatibility to long fragment length of dscNanoATAC/dscNanoCUT data (In the newer version, fragments longer than 1000 bp will be discarded)  

It takes ~ 30 min to finish the installation.

## Prepare references
1. Download the reference genome and chrom size file into the corresponding folder under [database](database). e.g.
```
    database/hg38/hg38.fa
    database/hg38/hg38.chrom.size
```
2. Build mimimap2 index:
```
    minimap2 -d ${ref}.mmi ${ref}.fa
```
3. Download the snp files into database/snp (OPTIONAL)
```
    # For GM12878
    wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz
    https://s3.eu-central-1.amazonaws.com/platinum-genomes/wget 2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz.tbi
```

## Usage
### 1. demultiplex of sequencing library
We use [Nanoplexer](https://github.com/hanyue36/nanoplexer) for library demultiplexing:
- inputs:
    - pass.fastq.gz: raw sequencing data;
    - barcode.fa: fasta files containing barcodes for demultiplex, see example in [barcode.fa](database/barcode/barcode.fa)
- output:
    - raw_data: containing the demultiplexed fastq files
- usage:
    ```
    sh scripts/nanoplexer.sh
    ```

### 2. preprocessing

- inputs:
    - raw_data（outputs from [Step.1](#1-demultiplex-of-sequencing-library)）
- outputs:
    - trim/{sample}_all_bc_merged.fastq.gz, fastq files containing valid reads 
    - trim/{sample}_fastq_stats.txt, statistics of intermediate fastq files
    - trim/{sample}_raw_barcode_stats.txt, statistics for each barcode
    - mapping/{sample}_duprm.bam, deduplicated bam files
    - mapping/{sample}_bc.bam.stats, statistics for raw bam files
    - mapping/{sample}_bc_q30.bam.stats, statistics for filtered bam files
    - mapping/{sample}_dedup.bam.stats, statistics for deduplicated bam files
    - fragment/{sample}_fragment.bed.gz, fragment files
    - arrow/{sample}.arrow, arrow files for ArchR analyses
- usage:
    ```
    snakemake -s pipeline/dscNanoSeq_preprocessing.smk -j 100 -k --profile slurm
    ```

### 3. downstream analyses
[ArchR](https://www.archrproject.com/index.html) is recommeded for downstream analyses.  


For allele-specific peaks analysie, please refer to [dscNanoSeq_haplotype_phasing.smk](pipeline/dscNanoSeq_haplotype_phasing.smk)

For co-accessibility analysis, please refer to [dscNanoSeq_coaccessibility.smk](pipeline/dscNanoSeq_coaccessibility.smk)

