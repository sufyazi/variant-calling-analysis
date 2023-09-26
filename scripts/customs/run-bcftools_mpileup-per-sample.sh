#!/usr/bin/env bash
# shellcheck disable=SC1091

# check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir/> <dataset_id> <tfbs_bedfile_names.txt> <out_dir/>"
    exit 1
fi

# set up variables
bam_dir=$1
data_id=$2
outfile=$4

# get dataset ID and sample ID
readarray -t bams < <(find "$bam_dir"/"$data_id" -mindepth 1 -maxdepth 1 -name "sample*" -type d)

# set count
count=1

# loop through each sample bam directory
for bam_samp in "${bams[@]}"; do
    # test count
    if [ "$count" -eq 1 ]; then
        echo "Count: $count"
        echo "Sample bam: $bam_samp"
        # get sample ID
        samp_id=$(basename "$bam_samp")
        echo "Sample ID: $samp_id"
        # get bam list
        readarray -t bam_list < <(find "$bam_samp" \( -name "*nodup.no_chrM_MT.bam" -o -name "*rep-merged.bam" \) -type f)
        echo "Bam list: " "${bam_list[@]}"
        # check if bam list exists
        if [ -n "${bam_list[*]}" ]; then
            echo "Found bam list. Proceeding..."
            # create a file to store bam list
            if [ ! -f /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-sample/"${data_id}"/"${samp_id}"_bam_list.txt ]; then
                mkdir -p /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-sample/"${data_id}"
                echo "${bam_list[@]}" > /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-sample/"${data_id}"/"${samp_id}"_bam_list.txt
            fi
            bam_file="/home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/input_files/mpileup_lists/per-sample/${data_id}/${samp_id}_bam_list.txt"
            echo "Submitting variant calling job..."
            # submit job
            echo qsub -v TF_LIST="$3",BAM_INP="${bam_file}",OUT_DIR="${outfile}",ID="${data_id}",SAMP="${samp_id}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/customs/run-bcftools_mpileup_submit.pbs
        fi
        count=$((count + 1))
    else
        echo "Reached count limit. Exiting..."
    fi
done
