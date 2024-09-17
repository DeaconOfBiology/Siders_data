#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=gunzip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=240gb
#SBATCH --time=30-00:00:00
#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out
work="/projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX"
reads="data/processed/reads/clean_reads"
cd "$work"
for i in "$reads"/*_1.fastq; do mv -- "$i" "${i%_1.fastq}_r1.fq"; done
for i in "$reads"/*_2.fastq; do mv -- "$i" "${i%_2.fastq}_r2.fq"; done
gzip "$reads"/*.fq
