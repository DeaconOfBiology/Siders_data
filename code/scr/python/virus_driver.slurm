#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=Derep_vcontigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240gb
#SBATCH --time=30-00:00:00
#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out

# Run your Bash script
#bash vcontig_dereplication.sh
#bash votu_read_mappinsg.sh
#bash votu_genomad_taxa.sh
bash viral_host_matching.sh
