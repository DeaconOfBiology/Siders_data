# vcontig_dereplication.sh
#!/bin/bash

###########################################################################
# This script uses blastx, anicalc.py, and aniclust.py to dereplicate the #
# viral contigs. This will result in the identification of vOTUs          #
# A working installation of NCBI BLAST+ as well as checkv is required to  #
# to run this script.                                                     #
###########################################################################

### This chunck of code loads that required modules form hpc and calls in the file_path.sh script

module purge
module load blast/2.9.0+

#Call in file paths
source ./file_paths.sh 


### This chunk checks for whether the blasts directory exists
# and creates it if it doesn't. This directory is necessary for the
# steps in the code below.

if [ -d "$viral_derep_blast_dbs" ]
then
    	echo "Blasts folder already exists, continuing..."
        echo
else
    	echo "Blasts folder doesn't exist, creating and continuing..."
        echo
	mkdir "$viral_derep_blast_dbs"
fi


#### The below code chunks do the following:
# Step 1: Concatinates all vcontigs from the different assemblies into one.
# Step 2: Uses makeblastdb from the BLAST+ software to build a blast database
# from the combined_virus_contigs_clean_headers.fa created in the step 1.
# Step 3: Blasts combined_virus_contigs_clean_headers.fa against the data base
# created in step 1.
# Step 4: Calculates pairwise ANI by combining local aligments between sequence
# pairs.
# Step 5: Performs CD-HIT-like clustering using the MIUVIG recommended-parameters 
# (95% ANI + 85% AF).

#Step 1
find "$assemblies"/day?-DO-0-*-virus-*/anvio/clean_fasta_file -type f -name "contigs_clean_headers.fa" -exec cat {} + > "$assemblies"/combined_virus_contigs_clean_headers.fa
#Step 2
makeblastdb -in "$assemblies"/combined_virus_contigs_clean_headers.fa -dbtype nucl -out "$blast_dbs"/viral_blast_db
#Step 3
blastn -query "$assemblies"/combined_virus_contigs_clean_headers.fa -db "$blast_dbs"/viral_blast_db -outfmt '6 std qlen slen' -max_target_seqs 10000 -o "$blast_dbs"/viral_blast.tsv -num_threads 32 
#Step 4
anicalc.py -i "$blast_dbs"/viral_blast.tsv -o "$blast_dbs"/viral_blast_ani.tsv
#Step 5
aniclust.py --fna "$assemblies"/combined_virus_contigs_clean_headers.fa --ani "$blast_dbs"/viral_blast_ani.tsv --out "$blast_dbs"/viral_blast_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0 