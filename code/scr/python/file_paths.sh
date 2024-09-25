# file_paths.sh
#!/bin/bash

#List of file paths
work="/projects/luo_lab/Siders_data"
assemblies="$work/data/processed/Assemblies"
reads="$work/data/processed/reads/clean_reads"
contig_db="$work/data/processed/anvio/contig_db"
profile_db="$work/data/processed/anvio/profile_db"
bowtie2="$work/data/processed/anvio/bowtie2"
merged_profile="$work/data/processed/anvio/merged_profile_db"
bins="$work/data/processed/bins"
final_bins="$work/data/processed/MAGs"
taxonomy="$work/data/processed/taxonomy"
kraken="$work/data/processed/taxonomy/contig/kraken"
kaiju="$work/data/processed/taxonomy/contig/kaiju"
viral_derep_blast_dbs="$work/data/processed/viral_derep_blast_dbs"
metacerberus_contigs="$work/data/processed/metabolic_predictions/contigs"
genomad_out="$work/data/processed/taxonomy/genomad/"
kaiju_db="/projects/luo_lab/Databases/kaijudb"

#File path to n needed python scripts
python_scripts="/projects/luo_lab/Siders_data/code/scr/python"