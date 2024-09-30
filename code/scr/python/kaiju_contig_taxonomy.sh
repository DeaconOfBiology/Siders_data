#!/bin/bash
module load mambaforge/4.14

mamba activate kaiju
# List of file paths
source ./file_paths.sh 

##################################################################
# This chunk of code runs kaiju to taxonomically ID cell contigs #
##################################################################
# Step 1. Taxonically identify using kaiju. Use the '-v' in the kaiju 
# step to keep the classification score.

#kaiju:
kaiju -t "$kaiju_db"/nodes.dmp \
    -f "$kaiju_db"/refseq/kaiju_db_refseq.fmi \
    -i "$assemblies"/clean_dreped_contigs_removed_cell.fa \
    -o "$kaiju"/kaiju.out -v

# Step 2. Add the names that correspond to the taxonomic ID:
kaiju-addTaxonNames -t "$kaiju_db"/nodes.dmp \
    -n "$kaiju_db"/names.dmp \
    -i "$kaiju"/kaiju.out  \
    -o "$kaiju"/kaiju_names.out \
    -r superkingdom,phylum,class,order,family,genus,species