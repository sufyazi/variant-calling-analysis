#!/usr/bin/env bash
# shellcheck disable=SC1091

# check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_list_dir/> <tfbs_bedfile_names.txt> <out_dir/>"
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
    count=0
    for tfbs in "${tfbs_list[@]}"; do
        while [ "$count" -lt 10 ]; do
            echo "Count: $count"
            echo "TF motif of interest: $tfbs"
            # find tfbs bed file
            bedfile=$(find /home/users/ntu/suffiazi/scratch/inputs/brca_tfbs_bedfiles -name "${tfbs}" -type f)
            echo "tfbs bed file location: ${bedfile}"
            motif_id="${tfbs%_tfbs_merged_matrix-brca_brca.bed}"
            echo "Motif ID: $motif_id"
            for bams in "${lists[@]}"; do
                echo "Bam file list: $bams"
                id_name=$(basename "${bams%.bam-list.txt}")
                echo "ID name: $id_name"
                # submit job
                echo qsub -v TF_FILE="${bedfile}",BAM_INP="${bams}",OUT_DIR="${outfile}",TF_ID="${motif_id}",ID="${id_name}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/customs/run-bcftools_mpileup_submit.pbs
            done
            count=$((count + 1))

        done
    done
else
    echo "No list.txt files found in $bam_list_dir"
    exit 1
fi