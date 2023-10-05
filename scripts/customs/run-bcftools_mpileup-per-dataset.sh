#!/usr/bin/env bash
# shellcheck disable=SC1091

# check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <bam_dir/> <dataset_id_list> <tfbs_bedfile_names.txt> <out_dir/>"
    exit 1
fi

# set up variables
bam_dir=$1
dataset_ids=$2
outfile=$4

# loop through the dataset IDs
readarray -t data_ids < "$dataset_ids"

for data_id in "${data_ids[@]}"; do
    echo "Dataset ID: $data_id"
    # find all bam files in the dataset directory
    readarray -t bams < <(find "$bam_dir"/"$data_id" \( -name "*nodup.no_chrM_MT.bam" -o -name "*nodup.rep-merged.bam" \) -type f)
    echo "Number of bam files: " "${#bams[@]}"
    # check if bam file variable is not empty
    if [ -n "${bams[*]}" ]; then
        echo "Found bam files. Proceeding..."
        # create a file to store bam list
        if [ ! -f /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-dataset/"${data_id}"_bam_list.txt ]; then
            mkdir -p /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-dataset
            for element in "${bams[@]}"; do
                echo "$element"
            done > /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-dataset/"${data_id}"_bam_list.txt
        fi
        bam_files="/home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-dataset/${data_id}_bam_list.txt"
        echo "Submitting PBS jobs sequentially..."
        # submit job
        qsub -v TF_LIST="$3",BAM_INP="${bam_files}",OUT_DIR="${outfile}",ID="${data_id}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/customs/run-bcftools_mpileup_parallel.pbs
    else
        echo "No bam files found in $bam_dir/$data_id"
        exit 1
    fi
done

