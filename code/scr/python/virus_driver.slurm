#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=Binning_merged_cell_assemblies
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240gb
#SBATCH --time=30-00:00:00
#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out

# Run your Bash script
bash vcontig_dereplication.sh
#bash votu_read_mapping.sh

