#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=cell_enrichment_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=200gb
#SBATCH --time=05-00:00:00
#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out

# Run your Bash script
#bash expand_reads.sh
#bash binning_individual.sh
#bash compress_reads.sh
#bash mapping_drep_mags.sh
#bash gtdbtk.sh
#bash kraken2.sh
bash kaiju_contig_taxonomy.sh
