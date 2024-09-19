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
# bowtie2-build "$assemblies"/representative_votus_3kbp.fa "$assemblies"/representative_votus_3kbp --threads 32
# anvi-gen-contigs-database -f  "$assemblies"/representative_votus_3kbp.fa -o "$contig_db"/representative_votus_3kbp.db -n 'dereplicated viral contig database' -T 32
# anvi-run-hmms -c "$contig_db"/representative_votus_3kbp.db -T 32 --just-do-it


### This code chunk mappes the viral enrichment reads to the representative viral contigs:
# mkdir "$bowtie2"/derep_viral_ve

# for reads in "$reads"/clean_day*virus*_r?.fq.gz
#     do
#     rdname=$(dirname ${reads})
#     rname=$(basename ${reads} _r1.fq.gz)
#     bowtie2 --threads 32 -x "$assemblies"/representative_votus_3kbp -1 ${rdname}/${rname}_r1.fq.gz -2 ${rdname}/${rname}_r2.fq.gz -S "$bowtie2"/derep_viral_ve/"$rname".sam
#     samtools view -F 2316 -bS "$bowtie2"/derep_viral_ve/"$rname".sam -o "$bowtie2"/derep_viral_ve/"$rname"-RAW.bam
#     samtools sort "$bowtie2"/derep_viral_ve/"$rname"-RAW.bam -o "$bowtie2"/derep_viral_ve/"$rname".bam
#     samtools index "$bowtie2"/derep_viral_ve/"$rname".bam
#     anvi-profile -i "$bowtie2"/derep_viral_ve/"$rname".bam -c "$contig_db"/representative_votus_3kbp.db -o "$profile_db"/derep_viral_ve/derep_"$rname"-profile

#     #echo "Removing  merged/"$rname".sam and  merged/"$rname"-RAW.bam files..."
#     rm -rf "$bowtie2"/derep_viral_ve/"$rname".sam "$bowtie2"/derep_viral_ve/"$rname"-RAW.bam
# done

# anvi-merge "$profile_db"/derep_viral_ve/derep_*virus*-profile/PROFILE.db -c "$contig_db"/representative_votus_3kbp.db -o "$merged_profile"/derep_viral_ve

# anvi-export-splits-and-coverages -p  "$merged_profile"/derep_viral_ve/PROFILE.db -c "$contig_db"/representative_votus_3kbp.db -o "$merged_profile"/derep_viral_ve --use-Q2Q3-coverages --report-contigs



### This code chunk mappes the viral enrichment reads to the representative viral contigs: 
# mkdir "$bowtie2"/derep_viral_ce

# for reads in "$reads"/clean_day*cell*_r1.fq.gz
#     do
#     rdname=$(dirname ${reads})
#     rname=$(basename ${reads} _r1.fq.gz)
#     bowtie2 --threads 32 -x "$assemblies"/representative_votus_3kbp -1 ${rdname}/${rname}_r1.fq.gz -2 ${rdname}/${rname}_r2.fq.gz -S "$bowtie2"/derep_viral_ce/"$rname".sam
#     samtools view -F 2316 -bS "$bowtie2"/derep_viral_ce/"$rname".sam -o "$bowtie2"/derep_viral_ce/"$rname"-RAW.bam
#     samtools sort "$bowtie2"/derep_viral_ce/"$rname"-RAW.bam -o "$bowtie2"/derep_viral_ce/"$rname".bam
#     samtools index "$bowtie2"/derep_viral_ce/"$rname".bam
#     anvi-profile -i "$bowtie2"/derep_viral_ce/"$rname".bam -c "$contig_db"/representative_votus_3kbp.db -o "$profile_db"/derep_viral_ce/derep_"$rname"-profile

#     #echo "Removing  merged/"$rname".sam and  merged/"$rname"-RAW.bam files..."
#     rm -rf "$bowtie2"/derep_viral_ce/"$rname".sam "$bowtie2"/derep_viral_ce/"$rname"-RAW.bam
# done

anvi-merge "$profile_db"/derep_viral_ce/derep_*cell*-profile/PROFILE.db -c "$contig_db"/representative_votus_3kbp.db -o "$merged_profile"/derep_viral_ce

anvi-export-splits-and-coverages -p  "$merged_profile"/derep_viral_ce/PROFILE.db -c "$contig_db"/representative_votus_3kbp.db -o "$merged_profile"/derep_viral_ce --use-Q2Q3-coverages --report-contigs