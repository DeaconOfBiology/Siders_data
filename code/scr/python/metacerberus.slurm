#!/bin/bash
#SBATCH --partition=Nebula
#SBATCH --job-name=MetaCerberus_vContigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=200gb
#SBATCH --time=48:00:00

#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out
cd /projects/luo_lab

## This code purges any module that may have been loaded, then loads the required module(s):
module purge
module load mambaforge/4.14
### This code calls in file paths:
work="/projects/luo_lab/Siders_data"
assemblies="$work/data/processed/Assemblies"
metacerberus_contigs="$work/data/processed/metabolic_predictions/contigs"
metacerberus_vContigs="$work/data/processed/metabolic_predictions/vContigs"

### Activate metacerberus mamba env:
mamba activate metacerberus

### This cunk of code uses metacerberus to annotate the vContigs (Need to test this when the cell contigs are finished annotating):
metacerberus.py --prodigal "$assemblies"/representative_votus.fa --hmm ALL --dir_out "$metacerberus_vContigs"

#When this is done, filter to only include contigs found in clean_dreped_contigs_removed_cell.fa. 
#metacerberus.py --prodigal "$assemblies"/day7-DO-0-13C-cell-enriched/anvio/clean_fasta_file/contigs_clean_headers.fa --hmm ALL --dir_out "$metacerberus_contigs"/day7-DO-0-13C-cell-enriched
#metacerberus.py --prodigal "$assemblies"/day7-DO-0-12C-cell-enriched/anvio/clean_fasta_file/contigs_clean_headers.fa --hmm ALL --dir_out "$metacerberus_contigs"/day7-DO-0-12C-cell-enriched
#metacerberus.py --prodigal "$assemblies"/day0-DO-0-env-virus-control/anvio/clean_fasta_file/contigs_clean_headers.fa --hmm ALL --dir_out "$metacerberus_contigs"/day0-DO-0-env-virus-control
