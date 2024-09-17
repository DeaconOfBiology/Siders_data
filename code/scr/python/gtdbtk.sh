# gtdbtk.sh
#!/bin/bash

# List of file paths
source ./file_paths.sh 

#Load mambaforge so we can use our mamba environments
module purge
module load mambaforge/4.14

#Activate gtdbtk mamba environment
mamba activate gtdbtk-2.1.1 #This is actually gtdbtk-2.3.2

#move to working directory
cd "$work"
#Run gtdbtk on MAGs
#gtdbtk classify_wf --genome_dir data/processed/MAGs/metawrap_50_10_bins --out_dir data/processed/Taxonomy/MAG -x fa --cpus 32 --pplacer_cpus 32 --skip_ani_screen
gtdbtk classify_wf --genome_dir "$final_bins"/drep/dereplicated_genomes --out_dir "$taxonomy"/MAG/drep -x fa --cpus 32 --pplacer_cpus 32 --skip_ani_screen
