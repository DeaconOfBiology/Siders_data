# votu_read_mapping.sh
#!/bin/bash

### This chunck of code loads that required modules form hpc and calls in the file_path.sh script

module purge
module load anvio/7
module load bowtie2/2.5.1
module load samtools/1.9

#Call in file paths
source ./file_paths.sh 


#Run bowtie2 and anvio profile
# bowtie2-build "$virsorter2"/derep_final-viral-combined.fa "$virsorter2"/derep_final-viral-combined --threads 32
# anvi-gen-contigs-database -f  "$virsorter2"/derep_final-viral-combined.fa -o "$contig_db"/derep_final-viral-combined.db -n 'dereplicated viral contig database' -T 32
# anvi-run-hmms -c "$contig_db"/derep_final-viral-combined.db -T 32 --just-do-it


# mkdir "$bowtie2"/derep_viral

# for reads in "$reads"/clean_day*virus*_1.fastq
#     do
#     rdname=$(dirname ${reads})
#     rname=$(basename ${reads} _1.fastq)
#     bowtie2 --threads 32 -x "$virsorter2"/derep_final-viral-combined -1 ${rdname}/${rname}_1.fastq -2 ${rdname}/${rname}_2.fastq -S "$bowtie2"/derep_viral/"$rname".sam
#     samtools view -F 2316 -bS "$bowtie2"/derep_viral/"$rname".sam -o "$bowtie2"/derep_viral/"$rname"-RAW.bam
#     samtools sort "$bowtie2"/derep_viral/"$rname"-RAW.bam -o "$bowtie2"/derep_viral/"$rname".bam
#     samtools index "$bowtie2"/derep_viral/"$rname".bam
#     anvi-profile -i "$bowtie2"/derep_viral/"$rname".bam -c "$contig_db"/derep_final-viral-combined.db -o "$profile_db"/derep_"$rname"-profile

#     #echo "Removing  merged/"$rname".sam and  merged/"$rname"-RAW.bam files..."
#     rm -rf "$bowtie2"/derep_viral/"$rname".sam "$bowtie2"/derep_viral/"$rname"-RAW.bam
# done


anvi-merge "$profile_db"/derep_*virus*-profile/PROFILE.db -c "$contig_db"/derep_final-viral-combined.db -o "$merged_profile"/derep_viral_ve

anvi-export-splits-and-coverages -p  "$merged_profile"/derep_viral_ve/PROFILE.db -c "$contig_db"/derep_final-viral-combined.db -o "$merged_profile"/derep_viral_ve --use-Q2Q3-coverages --report-contigs
