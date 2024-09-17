# expand_reads.sh
#!/bin/bash

#Call in file paths
source ./file_paths.sh 

#Unzip fastq files and reformate for metawrap
if [ -e "$reads"/*_r2.fq ]; then
    echo "Decompressing read files..."
    echo
    gunzip "$reads"/*.fq.gz
    
    #Change extention on forward reads
    echo "Reads have been decompressed. Now, changing the file extentions to be used in metawrap..."
    echo
    for i in "$reads"/*_r1.fq; do
        mv -- "$i" "${i%_r1.fq}_1.fastq"
    done

    #Change extention on reverse reads
    for i in "$reads"/*_r2.fq; do
        mv -- "$i" "${i%_r2.fq}_2.fastq" 
    done
else
    echo "reads are already decompressed. Moving on to next step..."
    echo
fi
