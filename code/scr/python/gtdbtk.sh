# gtdbtk.sh
#!/bin/bash

# List of file paths
source ./file_paths.sh 

#move to working directory
cd "$work"

########################################################
# Load mambaforge so we can use our mamba environments #
########################################################
module purge
module load mambaforge/4.14

#####################################
# Activate gtdbtk mamba environment #
#####################################
mamba activate gtdbtk-2.1.1 #This is actually gtdbtk-2.3.2

#######################################################
# This chunk of code runs the "classify_wf" of gtdbtk #
#######################################################
gtdbtk classify_wf --genome_dir "$final_bins"/drep/dereplicated_genomes \
    --out_dir "$taxonomy"/MAG -x fa \
    --cpus 32 \
    --pplacer_cpus 32 \
    --skip_ani_screen


################################################################
# This chunck of code runs the 'gtdb_to_ncbi_majority_vote.py' # 
# script to convert gtdbtk taxonomy to ncbi taxonomy           #
################################################################
python3 "$python_scripts"/gtdb_to_ncbi_majority_vote.py \
    --gtdbtk_output_dir "$taxonomy"/MAG \
    --bac120_metadata_file "$taxonomy"/MAG/metadata/bac120_metadata_r214.tar.gz \
    --ar53_metadata_file "$taxonomy"/MAG/metadata/ar53_metadata_r214.tar.gz \
    --output_file "$taxonomy"/MAG/gtdbtk_to_ncbi_taxonomy.csv \
    --gtdbtk_prefix gtdbtk

#########################################################
# This chunk runs programs of gtdbtk that were not ran  #
# in "classify_wf". The outputs from these programs are #
# needed for IPHoP.                                     #
#########################################################

gzip -d "$taxonomy"/MAG/align/*.msa.fasta.gz

gtdbtk infer --msa_file "$taxonomy"/MAG/align/gtdbtk.ar53.msa.fasta \
    --prefix gtdbtk.ar53 --out_dir "$taxonomy"/MAG/

gtdbtk root --input_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.ar53.unrooted.tree \
    --outgroup_taxon p__Altiarchaeota \
    --output_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.ar53.rooted.tree

gtdbtk decorate --input_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.ar53.rooted.tree \
    --output_tree "$taxonomy"/MAG/infer/gtdbtk.ar53.decorated.tree

gtdbtk infer --msa_file "$taxonomy"/MAG/align/gtdbtk.bac120.msa.fasta \
    --prefix gtdbtk.bac120 --out_dir "$taxonomy"/MAG/
    
gtdbtk root --input_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.bac120.unrooted.tree \
    --outgroup_taxon p__Patescibacteria \
    --output_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.bac120.rooted.tree

gtdbtk decorate --input_tree "$taxonomy"/MAG/infer/intermediate_results/gtdbtk.bac120.rooted.tree \
    --output_tree "$taxonomy"/MAG/infer/gtdbtk.bac120.decorated.tree