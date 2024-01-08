#!/bin/bash
#shellcheck disable=SC2016
# Check if the required arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <tfbs_motif_prefix_list> <fpscore_matrix_root_dir> <output_bed_dir>"
    exit 1
fi

# Read the input arguments
input_prefix=$1
root_dir=$2
output_dir=$3

# Loop through each prefix in the input file
while IFS= read -r prefix; do
    # Use the find command to search for files with the pattern
    file=$(find "$output_dir" -name "$prefix*.bed" -type f -print -quit)
    
    # check if the output file exists
    if [ -n "$file" ]; then
        echo "Output file ${file} already exists. Skipping..."
        continue
    else
        echo "Processing $prefix..."
        # Use find to search for files with the given prefix and .txt suffix
        matrix=$(find "$root_dir" -name "${prefix}_*matrix.txt" -type f)
        # check if the matrix variable is empty
        if [ -z "$matrix" ]; then
            echo "Matrix file for $prefix not found. Skipping..."
            continue
        else
            echo "Matrix file found: $matrix"
            # process the matrix file; grab the first 4 cols, add # to the header line for bedtools, then use sed to sort the bed file by chr and start position, ignoring the header line (1q)
            awk -v OFS='\t' '{print $1, $2, $3, $4}' "${matrix}" | sed '1s/^/#/' | (sed -u 1q; sort -k1,1V -k2,2n) > "$output_dir"/"${prefix}"_diffmode_TCGA-BRCA_fpscore_regions.bed
        fi
    fi
done < "$input_prefix"