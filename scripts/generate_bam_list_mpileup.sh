#!/usr/bin/env bash

# This script generates a list of bam files for use with mpileup

# Usage: generate_bam_list_mpileup.sh <path to root directory of bam files> <data id list>

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: generate_bam_list_mpileup.sh <path to root directory of bam files> <data id list>"
    exit 1
fi

# Set variables
BAM_DIR=$1
DATA_LIST=$2

# Traverse through data id list
readarray -t DATA_ID < "$DATA_LIST"
echo "List of data ids:" "${DATA_ID[@]}"

# Generate bam list
for i in "${DATA_ID[@]}"; do 
    # check subdirectories for bam files
    if [ -d "${BAM_DIR}/$i" ]; then
        echo "Generating bam list for $i..."
        find "${BAM_DIR}/$i" -name "*.bam" >> /home/users/ntu/suffiazi/scratch/inputs/mpileup_list/"$i".bam-list.txt
    else
        echo "No bam files found for $i"
    fi
done