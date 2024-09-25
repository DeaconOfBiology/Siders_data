#!/bin/bash
module load mambaforge/4.14

mamba activate kaiju
# List of file paths
source ./file_paths.sh 

##################################################################
# This chunk of code runs kaiju to taxonomically ID cell contigs #
##################################################################
kaiju -t "$kaiju_db"/nodes.dmp -f "$kaiju_db"/refseq/kaiju_db_refseq.fmi -i "$assemblies"/clean_dreped_contigs_removed_cell.fa -o "$kaiju"/kaiju.out.krona
