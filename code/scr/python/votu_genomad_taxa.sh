#!/bin/bash
#Viral taxonomic identification

module load mambaforge/4.14
mamba activate genomad

# List of file paths
source ./file_paths.sh 


#Running genomad for taxonomic identification
genomad end-to-end --cleanup "$assemblies"/representative_votus.fa "$genomad_out" /users/troger50/genomad_db
