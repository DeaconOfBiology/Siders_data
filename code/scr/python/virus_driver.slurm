#!/bin/bash
#SBATCH --partition=Nebula
#SBATCH --job-name=Derep_vcontigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240gb
#SBATCH --time=02-00:00:00
#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out

# Run your Bash script
#bash vcontig_dereplication.sh
bash votu_read_mapping.sh

