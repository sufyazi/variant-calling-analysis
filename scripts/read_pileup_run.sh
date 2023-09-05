#!/usr/bin/env bash
# shellcheck disable=SC1091

# check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list_dir/> <tfbs_candidate.txt> <out_dir/>"
    exit 1
fi

# set up variables
bam_list_dir=$1
outfile=$3

# find all list.txt files in the directory
readarray -t lists < <(find "$bam_list_dir" -name "*.txt" -type f)
echo "Bam file lists loaded: " "${lists[@]}"

# load up tfbs list
readarray -t tfbs_list < "$2"
echo "List of TF motifs of interest: " "${tfbs_list[@]}"

# check if list.txt files exist
if [ -n "${lists[*]}" ] && [ -n "${tfbs_list[*]}" ]; then
    echo "Found list.txt files and candidate TF motifs of interest. Proceeding..."
    for tfbs in "${tfbs_list[@]}"; do
        echo "TF motif of interest: $tfbs"
        tfbs_loc=$(find "/scratch/users/ntu/suffiazi/outputs/brca-subtype-basal-UP" -name "${tfbs}*.bed" -type f)
        echo "tfbs bed file location: ${tfbs_loc}"
        for bams in "${lists[@]}"; do
            echo "Bam file list: $bams"
            id_name=$(basename "${bams%.bam-list.txt}")
            # submit job
            qsub -v TF_FILE="${tfbs_loc}",BAM_INP="${bams}",OUT_DIR="${outfile}",TF_NAME="${tfbs}",ID="${id_name}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/bcftools_mpileup_submit.pbs
        done
    done
else
    echo "No list.txt files found in $bam_list_dir"
    exit 1
fi