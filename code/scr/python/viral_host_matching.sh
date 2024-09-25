#!/bin/bash
#Viral host matching

module load mambaforge/4.14
mamba activate iphop_env

# List of file paths
source ./file_paths.sh 

#Running iPHoPs
#Run on iPHop database:
iphop predict --fa_file "$assemblies"/representative_votus.fa --db_dir /projects/luo_lab/Databases/iphop_db/Aug_2023_pub_rw/ --out_dir "$work"/data/processed/vh_matching/iphop/iphop_db -t 32

#Make personal database:
iphop add_to_db --fna_dir "$final_bins"/drep/dereplicated_genomes/ --gtdb_dir "$taxonomy"/MAG/ --out_dir /projects/luo_lab/Databases/Siders_iphop_inhouse_db --db_dir /projects/luo_lab/Databases/iphop_db/Aug_2023_pub_rw/

#Run on personal data based:
iphop predict --fa_file "$assemblies"/representative_votus.fa --db_dir /projects/luo_lab/Databases/Siders_iphop_inhouse_db/ --out_dir "$work"/data/processed/vh_matching/iphop/iphop_inhouse_db -t 32

