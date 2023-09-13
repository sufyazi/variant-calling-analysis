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
        tfbs_loc=$(find "/scratch/users/ntu/suffiazi/outputs/BRCA-diff-footprinting/processed-tfbs" -name "${tfbs}*4col.bed" -type f)
        readarray -t tfbs_loc <<< "$tfbs_loc"
        echo "tfbs bed file location: " "${tfbs_loc[@]}"
        for bams in "${lists[@]}"; do
            echo "Bam file list: $bams"
            id_name=$(basename "${bams%.bam-list.txt}")
            echo "ID name: $id_name"
            # check if tfbs_loc is more than 1 in length
            echo "Length of tfbs_loc array: ${#tfbs_loc[@]}"
            if [ "${#tfbs_loc[@]}" -gt 1 ]; then
                echo "Motif binding site file array is more than 1 in length. Proceeding with for loop..."
                for txt in "${tfbs_loc[@]}"; do
                    echo "Motif BS file: $txt"
                    motif_id=$(basename "${txt%_binding_sites-basal-UP_4col.bed}")
                    echo "Motif ID: $motif_id"
                    # submit job
                    qsub -v TF_FILE="${txt}",BAM_INP="${bams}",OUT_DIR="${outfile}",TF_ID="${motif_id}",ID="${id_name}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/run-bcftools_mpileup_submit.pbs
                done
            elif [ "${#tfbs_loc[@]}" -eq 1 ]; then
                echo "There is precisely one motif TFBS file in the array. Proceeding with job submission..."
                motif_id=$(basename "${tfbs_loc[0]%_binding_sites-basal-UP_4col.bed}")
                echo "Motif ID: $motif_id"
                # submit job
                qsub -v TF_FILE="${tfbs_loc[0]}",BAM_INP="${bams}",OUT_DIR="${outfile}",TF_ID="${motif_id}",ID="${id_name}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/scripts/run-bcftools_mpileup_submit.pbs
            fi
        done
    done
else
    echo "No list.txt files found in $bam_list_dir"
    exit 1
fi