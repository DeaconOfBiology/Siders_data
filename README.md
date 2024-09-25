# SIP-metagenomic analysis of cell and viral enrichment from Siders pond

## All the data and code needed to analysis the metagenomic data from DNA-SIP incubations of Siders pond sample

# Summary
 Viruses impact aquatic biogeochemistry via host cell lysis and rewiring host metabolic processes. However, their functional role in aquatic dark carbon cycling is mostly unexplored. Here, we established a method to identify active viruses targeting primary producers (chemoautotrophs) using <sup>13</sup>C-DNA stable isotopic probing combined with metagenomics. 

 # Overview of directories
```bash
Siders_data
├── README  # the top level description of the content (this doc)
├── LICENSE # the license for this project
├── code
│   ├── notebooks
│   │   └── ExcessAtomFraction_calculations.ipynb # jupyter notebook to calculate relative abundance and excess atom fraction of cell and viral contigs
│   └── scr
│       ├── python
│       │   ├── anicalc.py                      # calculate the ani of the viral contigs
│       │   ├── aniclust.py                     # cluster viral contigs based on ani
│       │   ├── binning_individual.sh           # bin assemblies
│       │   ├── compress_reads.sh               # guzip '.fq' files
│       │   ├── driver.slurm                    # driver script for the cell contig analysis
│       │   ├── expand_reads.sh                 # decompress '.fq.gz' files
│       │   ├── file_paths.sh                   # file path identification
│       │   ├── gtdbtk.sh                       # taxonomically ID MAGs
│       │   ├── gtdb_to_ncbi_majority_vote.py   # translate gtdb taxonomy to ncbi taxonomy
│       │   ├── kaiju_contig_taxonomy.sh        # cell taxonomic identification with kaiju
│       │   ├── kraken2.sh                      # cell taxonomic identification with kraken2
│       │   ├── mapping_drep_mags.sh            # cell contig read mapping
│       │   ├── metacerberus.slurm              # metabolic predictions for all contigs
│       │   ├── pipe_line.slurm                 # script to run the snakemake pipeline
│       │   ├── vcontig_dereplication.sh        # dereplicate viral contigs
│       │   ├── viral_host_matching.sh          # viral host analysis with IPHoP
│       │   ├── virus_driver.slurm              # driver script for the viral contig analysis
│       │   ├── votu_genomad_taxa.sh            # taxonomic identification with genomad
│       │   └── votu_read_mapping.sh            # viral contig read mapping
│       └── R
│           └── EAF_plots.r # script to create EAF vs DNA density figures
├── data
│   ├── processed   # cleaned data, will not be altered once created
│   ├── raw         # raw data, will not be altered
│   └── temp_files  # temporary data created and deleted during workflow
├── smk_workflow
│   ├── config  # configuration file(s) for snakemake pipeline
│   ├── reports # report output of snakemake pipeline status
│   └── rules   
│       ├── preprocess_reads.smk                     # trim QC raw reads
│       ├── create_assemblies.smk                    # assemble cell enrichment (ce)
│       ├── cell_enriched_vcontig_identification.smk # idenitifying viral contigs in ce
│       ├── cell_enriched_vcontig_reads.smk          # ID reads mapping to viral contigs in ce 
│       └── virus_enriched_votu_identification.smk   # create and identify viral contigs
├── results 
│   ├── figures # graphs, likely designated for manuscript figures
│   └── tables  # tables out puts to be used in R for figures and other analyses
└── submission
    ├── manuscript.Rmd 
    ├── manuscript.pdf
    ├── manuscript.tex
    ├── references.bib
    ├── figures         
    │   ├── text         # final in text figures for submission
    │   └── supplemental # final supplemental figures for submission
    └── tables 
        ├── text         # final in text tables for submission
        └── supplemental # final supplemental tables for submission
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
* Tools needed for snakemake pipeline:
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
* Additional tools needed:
    * gtdbtk 2.3.2
    * kaiju 1.10.1
    * kraken2
    * metawrap
    * drep
    * genomad
    * metacerberus

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
