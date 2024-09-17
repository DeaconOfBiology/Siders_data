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
    * BBMap
    * Bowtie 2
    * metaSPAdes

## Running the ananlysis
```
```
