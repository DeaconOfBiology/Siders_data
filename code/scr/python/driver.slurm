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
bash expand_reads.sh
bash binning_individual.sh
bash compress_reads.sh
bash mapping_drep_mags.sh
bash gtdbtk.sh
bash kraken2.sh
