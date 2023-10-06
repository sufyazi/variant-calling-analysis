#!/usr/bin/env bash
# shellcheck disable=SC1091

# variables
bedfile=$1
bams=$2
outfile=$3
data_id=$4

echo "TF motif of interest: $bedfile"
# get motif ID
motif_id="${bedfile%_tfbs_merged_matrix-brca_brca.bed}"
echo "Motif ID: $motif_id"

# check if output vcf file already exists
if [ -f "${outfile}"/"${data_id}"_"${motif_id}"_qualgt10.var.flt.vcf ]; then
    echo "Output vcf file already exists. Skipping..."
else
    echo "Output vcf file does not exist. Proceeding..."
    # get bedfile path
    bed_path=$(find /scratch/users/ntu/suffiazi/inputs/brca_tfbs_bedfiles -name "${bedfile}" -type f)
    if [ -n "$bed_path" ]; then
        echo "Found bed file. Proceeding..."
        echo "Bed file path: $bed_path"
        # call and filter variants from bam files of BRCA datasets
        if bcftools mpileup -Ou -f /scratch/users/ntu/suffiazi/inputs/references/gatk4/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -T "${bed_path}" -b "${bams}" | bcftools call -Ou -mv | bcftools filter -i'QUAL>10' > "${outfile}"/"${data_id}"_"${motif_id}"_qualgt10.var.flt.vcf; then
                echo "Variant calling and filtering completed."
            else
                echo "Variant calling and filtering failed."
                exit 1
            fi
        else
            echo "Bed file ${bedfile} not found. Check again. Exiting..."
            exit 1
        fi
    fi

echo "Job completed."