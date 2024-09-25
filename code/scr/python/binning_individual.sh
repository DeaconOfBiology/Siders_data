# binning_individual.sh
#!/bin/bash


# List of file paths
source ./file_paths.sh 

# cd into working directory
cd "$work"

# Load required modules
module purge
module load mambaforge/4.14
# mamba activate metawrap-env

#Create list to assembly files
assembly_files=$(find "$assemblies" -type f -name "contigs_clean_headers.fa") 

# Iterate over each assembly file
for assembly in $assembly_files; do
    # Extract the base name of the assembly file for creating a unique output directory
    parent_dir=$(dirname "$(dirname "$(dirname "$assembly")")") 
    assembly_name=$(basename "$parent_dir")
    output_dir="$bins"/"$assembly_name"

    # Run metawrap binning
    metawrap binning -o "$output_dir" -t 32 -a "$assembly" --metabat2 --maxbin2 --concoct --universal "$reads"/clean_"$assembly_name"*.fastq
    metawrap bin_refinement -o "$final_bins"/"$assembly_name" -t 32 -A "$output_dir"/metabat2_bins/ -B "$output_dir"/maxbin2_bins/ -C "$output_dir"/concoct_bins/ -c 50 -x 10
done
mamba deactivate

#Dereplicate MAGs across assemblies:
#-Checkm output
# -Add assembly name to bin files as a prefix
for f in "$final_bins"/day*/metawrap_50_10_bins/*.fa; do
  main_path=$(dirname "$f") 
  prefix_path=$(dirname "$main_path")
  prefix=$(basename "$prefix_path")
  filename=$(basename "$f")
  new_filename="$prefix"_"$filename"
  mv "$f" "$main_path/$new_filename"
done

# -Add the same prefix to the checkm output
for dir in "$final_bins"/day*; do
  # Extract the prefix from the directory name
  name=$(basename "$dir")
  # Update the checkm output file with the prefix
  awk -v prefix="$name" -F'\t' 'BEGIN {OFS="\t"} {if (NR == 1) {print $0} else {$1 = prefix "_" $1; print $0}}' "$dir"/metawrap_50_10_bins.stats >  "$dir"/prefix_metawrap_50_10_bins.stats
done

# -Set up checkm output to work with drep. 
mkdir "$final_bins"/drep
#  -Find the first file and extract header
header="genome,completeness,contamination"

find "$final_bins"/day* -name "prefix_metawrap_50_10_bins.stats" -exec sh -c 'tail -n +2 "$0" | cut -f1-3' {} \; \
  | sed 's/^\([^\t]*\)/\1.fa/' \
  | tr '\t' ',' > "$final_bins"/drep/dRep_temp_data

#  -Combine header and data
echo "$header" | tr '\t' ',' > "$final_bins"/drep/dRep.genomeInfo
cat "$final_bins"/drep/dRep_temp_data >> "$final_bins"/drep/dRep.genomeInfo


#  -Clean up
rm "$final_bins"/drep/dRep_temp_data 

# -Copy all bins to one directory
mkdir "$final_bins"/all_mags
find "$final_bins"/day* -type f -name "day*.fa" -exec cp {} "$final_bins"/all_mags \;

# -Run dRep
module purge
module load mambaforge/4.14
mamba activate drep

dRep dereplicate --genomeInfo "$final_bins"/drep/dRep.genomeInfo \
    --S_algorithm fastANI \
    -comp 50 \
    --SkipMash \
    -g "$final_bins"/all_mags/*.fa \
    -p 32 \
    "$final_bins"/drep/




