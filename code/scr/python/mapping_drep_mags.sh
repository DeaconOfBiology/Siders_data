# mapping_drep_mags.sh
#!/bin/bash

#Call in file paths
source ./file_paths.sh 

#CD to working directory
cd "$work"
#Load required modules
module load anvio/7
module load bowtie2/2.5.1
module load samtools/1.9

#################################################
# Merge cell_enriched reads into one fasta file #
#################################################
find "$assemblies"/day?-DO-0-*-cell-*/anvio/clean_fasta_file -type f -name "contigs_clean_headers.fa" -exec cat {} + > "$assemblies"/combined_cell_contigs_clean_headers.fa

##########################################################################
# Remove contigs that belong to MAGs that are not the represenative MAGs #
##########################################################################
    #1. Create a list of MAGs that are not the representive MAGs
comm -23 <(ls "$final_bins"/all_mags | sort) <(ls "$final_bins"/drep/dereplicated_genomes | sort) > "$final_bins"/removed_MAGs.txt
    
    #2. Clean up removed_MAGs.txt fasta names
sed -e 's/\.fa$//' -e 's/bin\./bin_/' "$final_bins"/removed_MAGs.txt > "$final_bins"/modified_removed_MAGs.txt

    #3. Create a list of contigs in non-representive MAGs
        # Extract fasta filenames from removed_MAGs.txt and sort them
fasta_files=$(awk '{print $1}' "$final_bins"/modified_removed_MAGs.txt | sort)
        # Extract contigs associated with these fasta files from all_binning_results.txt
awk -v files="$fasta_files" '
BEGIN {
    split(files, file_array, "\n");
    for (i in file_array) file_set[file_array[i]] = 1;
}
{
    if ($2 in file_set) print $1;
}' "$final_bins"/all_binning_results.txt > "$final_bins"/contigs_in_removed_MAGs.txt

    #4. Remove contigs from "$assemblies"/combined_cell_contigs_clean_headers.fa found in the "$final_bins"/contigs_in_removed_MAGs.txt file
        # Extract contig identifiers into an array
contigs_file="$final_bins/contigs_in_removed_MAGs.txt"
        # Use awk to filter out contigs
awk '
# Read the contig file and build the set of contig headers to remove
BEGIN {
    while ((getline line < "'"$contigs_file"'") > 0) {
        contig_set[">" line] = 1;
    }
    close("'"$contigs_file"'");
}
# When processing the fasta file, track the current contig
{
    if (substr($0, 1, 1) == ">") {
        current_contig = $0;
    }
    if (current_contig in contig_set) {
        current_contig = "";
    } else if (current_contig != "") {
        print $0;
    }
}' "$assemblies"/combined_cell_contigs_clean_headers.fa > "$assemblies"/clean_dreped_contigs_removed_cell.fa

##################################
# Run bowtie2 and anvio profiles #
##################################
bowtie2-build "$assemblies"/clean_dreped_contigs_removed_cell.fa "$assemblies"/clean_dreped_contigs_removed_cell --threads 32
anvi-gen-contigs-database -f  "$assemblies"/clean_dreped_contigs_removed_cell.fa -o "$contig_db"/clean_dreped_contigs_removed_cell.db -n 'drep contigs removed cell' -T 32
anvi-run-hmms -c "$contig_db"/clean_dreped_contigs_removed_cell.db -T 32 --just-do-it

mkdir "$bowtie2"/drep_contigs_removed

for reads in "$reads"/clean_day*cell*_r1.fq.gz
    do
    rdname=$(dirname ${reads})
    rname=$(basename ${reads} _r1.fq.gz)
    bowtie2 --threads 32 -x "$assemblies"/clean_dreped_contigs_removed_cell -1 ${rdname}/${rname}_r1.fq.gz -2 ${rdname}/${rname}_r2.fq.gz -S "$bowtie2"/drep_contigs_removed/"$rname".sam
    samtools view -F 2316 -bS "$bowtie2"/drep_contigs_removed/"$rname".sam -o "$bowtie2"/drep_contigs_removed/"$rname"-RAW.bam
    samtools sort "$bowtie2"/drep_contigs_removed/"$rname"-RAW.bam -o "$bowtie2"/drep_contigs_removed/"$rname".bam
    samtools index "$bowtie2"/drep_contigs_removed/"$rname".bam
    anvi-profile -i "$bowtie2"/drep_contigs_removed/"$rname".bam -c "$contig_db"/clean_dreped_contigs_removed_cell.db -o "$profile_db"/drep_contigs_removed_"$rname"-profile

    #echo "Removing  merged/"$rname".sam and  merged/"$rname"-RAW.bam files..."
    rm -rf "$bowtie2"/drep_contigs_removed/"$rname".sam "$bowtie2"/drep_contigs_removed/"$rname"-RAW.bam
done


anvi-merge "$profile_db"/drep_contigs_removed*-profile/PROFILE.db -c "$contig_db"/clean_dreped_contigs_removed_cell.db -o "$merged_profile"/drep_contigs_removed

anvi-export-splits-and-coverages -p "$merged_profile"/drep_contigs_removed/PROFILE.db -c "$contig_db"/clean_dreped_contigs_removed_cell.db -o "$merged_profile"/drep_contigs_removed --use-Q2Q3-coverages --report-contigs
