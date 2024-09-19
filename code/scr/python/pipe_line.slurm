#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=Siders_Pipe_line
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=300gb
#SBATCH --time=30-00:00:00

#SBATCH --error=./log_reports/%x_%j.err
#SBATCH --output=./log_reports/%x_%j.out

#Load mambaforge so we can use our mamba environments
module purge
module load mambaforge/4.14
# Run Snakemake inside the SLURM job
mamba activate snakemake_pipe

# Navigate to the directory containing your Snakefile and scripts
cd /projects/luo_lab/Siders_data/smk_workflow/

# If you recieve a notice that directories or files are locked by snakemake, run the below script:
    #snakemake --unlock --cores all

#Create a map of the analysis
snakemake --dag all_run | dot -Tsvg > smk_workflow/reports/dag.svg 

#Run the pipeline
snakemake all_run --cores all --use-conda --resources mem_mb=250000 --scheduler greedy --rerun-incomplete

# If you need to test the pipeline, run the below in the smk_workflow/ directory:
    #snakemake all_run -n -r
