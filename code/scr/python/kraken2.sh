# kracken2.sh
#!/bin/bash
# List of file paths
source ./file_paths.sh 

#load modules
module purge
module load mambaforge/4.14

mamba activate metawrap-env

cd "$work"

kraken2 --use-names --db /projects/luo_lab/Databases/MY_KRAKEN2_DB --threads 10 "$assemblies"/clean_dreped_contigs_removed_cell.fa --confidence 0.5 > "$kraken"/kraken2_output.txt