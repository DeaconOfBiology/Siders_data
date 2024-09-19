# SIP-metagenomic analysis of cell and viral enrichment from Siders pond

## All the data and code needed to analysis the metagenomic data from DNA-SIP incubations of Siders pond sample

# Summary
 Viruses impact aquatic biogeochemistry via host cell lysis and rewiring host metabolic processes. However, their functional role in aquatic dark carbon cycling is mostly unexplored. Here, we established a method to identify active viruses targeting primary producers (chemoautotrophs) using <sup>13</sup>C-DNA stable isotopic probing combined with metagenomics. 

 # Overview of directories
```bash
.
├── code
│   ├── notebooks
│   └── scr
│       ├── python
│       └── R
├── data
│   ├── processed
│   │   └── smk_out
│   │       ├── Assemblies
│   │       │   ├── day0-DO-0-env-cell-control
│   │       │   │   ├── anvio
│   │       │   │   ├── bam
│   │       │   │   ├── megahit
│   │       │   │   └── virsorter_first
│   │       │   ├── day0-DO-0-env-virus-control
│   │       │   │   ├── anvio
│   │       │   │   ├── checkv
│   │       │   │   ├── megahit
│   │       │   │   ├── virsorter_final
│   │       │   │   └── virsorter_first
│   │       │   ├── day7-DO-0-12C-cell-enriched
│   │       │   │   ├── anvio
│   │       │   │   ├── bam
│   │       │   │   ├── megahit
│   │       │   │   └── virsorter_first
│   │       │   ├── day7-DO-0-12C-virus-enriched
│   │       │   │   ├── anivo
│   │       │   │   ├── checkv
│   │       │   │   ├── megahit
│   │       │   │   ├── virsorter_final
│   │       │   │   └── virsorter_first
│   │       │   ├── day7-DO-0-13C-cell-enriched
│   │       │   │   ├── anvio
│   │       │   │   ├── bam
│   │       │   │   ├── megahit
│   │       │   │   └── virsorter_first
│   │       │   ├── day7-DO-0-13C-virus-enriched
│   │       │       ├── anvio
│   │       │       ├── checkv
│   │       │       ├── megahit
│   │       │       ├── virsorter_final
│   │       │       └── virsorter_first
│   │       ├── bins
│   │       │   └── work_files
│   │       └── reads
│   │           ├── clean_reads
│   │           ├── merged_reads
│   │           ├── vcontig_cell_enriched_reads
│   │           └── vcontig_env_cell_enriched_reads
│   ├── raw
│   └── temp_files
├── documentation
└── smk_workflow
    ├── config
    ├── reports
    └── rules
```
# How to regenerate this repository

## Dependencies and locations
* R packages:
    * data.table
    * tidyverse
    * KEGGREST
    * ggplot2
    * pathview
    * ape
* Python packages:
    * pandas
    * sklearn
    * re
    * time
* Tools needed for snakemake pipeline.
    * snakemake 7.31.1
    * BBMap 38.97
    * Bowtie 2 2.5.1
    * samtools 1.18
    * megahit 1.2.9
    * multiqc 1.19
    * fastqc 0.12.1
    * quast 5.2.0
    * anvio 7
    * virsorter2 2.2.4
    * checkv 1.0.1

## Running the ananlysis
```
# Step 1: Run snakemake workflow:
cd ./code/scr/python/
sbatch pipe_line.slurm

# Step 2: Bin assemblies, dereplicate bins (final MAGs), map cell enriched reads to assemblies, 
# run GTDB on MAGs, and run kraken2 on assemblies:
sbatch driver.slurm

# Step 3: Dereplicate viral contigs to identify final vOTUs and map reads from the viral and cell enrichment to vOTU's
sbatch virus_driver.slurm

# Step 4: Annotate vOTUs and cell enrichement contigs:
sbatch metacerberus.slurm

# Step 5: Calculate EAF for cell enrichment contigs, vOTUs in viral enrichment, and vOTUs in cell enrichment"
cd ~/Siders_data/code/notebooks
jupyter nbconvert --to script --execute --stdout ExcessAtomFraction_calculations.ipynb | python

# Step 6: Write the paper:


```
