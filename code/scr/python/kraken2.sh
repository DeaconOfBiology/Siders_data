# kracken2.sh
#!/bin/bash
# List of file paths
source ./file_paths.sh 

#load modules
module purge
module load mambaforge/4.14

mamba activate metawrap-env

cd "$work"

#kraken2 --threads 10 --db /projects/luo_lab/Databases/MY_KRAKEN2_DB /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/Assemblies/cell_retentate/scaffolds-clean_headers.fasta --output /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/kraken2_out/cells --confidence 0.5
#kraken2 --db /projects/luo_lab/Databases/MY_KRAKEN2_DB /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/Assemblies/cell_retentate/scaffolds-clean_headers.fasta  --output /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/kraken2_out --report
#metawrap kraken2 --no-preload -o /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/kraken2_out/cell -t 32  /projects/luo_lab/Rogers_SidersViralAnalysis_XXXX_20XX/data/processed/Assemblies/cell_retentate/scaffolds-clean_headers.fasta 
kraken2 --use-names --db /projects/luo_lab/Databases/MY_KRAKEN2_DB --threads 10 "$assemblies"/clean_dreped_contigs_removed_cell.fa --confidence 0.5 > "$kraken"/kraken2_output.txt