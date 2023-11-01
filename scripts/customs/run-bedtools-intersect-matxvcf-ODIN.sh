#!/usr/bin/env bash
# shellcheck disable=SC1091

# set up variables
motif_prefixes=$1
motif_regions=$2
vcf_root_dir=$3
output_dir=$4

# load motif prefixes from the input text file
mapfile -t motif_prefixes < "$motif_prefixes"

# loop through each motif prefix
for motif_prefix in "${motif_prefixes[@]}"; do
    echo "Motif to process: $motif_prefix. Finding the associated vcf files..."
    # find specific motif vcf files
    readarray -t vcf_files < <( find "${vcf_root_dir}" -name "*_${motif_prefix}*.vcf" -type f)
    # construct an array of data_id from basename of vcf files
    data_ids=()
    for vcf in "${vcf_files[@]}"; do
        data_ids+=("$(basename "${vcf}" | cut -d'_' -f1)")
    done
    echo "Dataset ID in order of appearance: " "${data_ids[@]}"
    # find specific motif fps region file
    region_file=$(find "${motif_regions}" -name "*${motif_prefix}*stranded.txt" -type f)
    if [[ -z "${vcf_files[0]}" || -z "$region_file" ]]; then
        echo "VCF files or motif binding region file not found for ${motif_prefix}. Skipping..."
        continue
    fi
    # construct output filename
    output="${output_dir}/var-positions-per-subtype/${motif_prefix}_variant-overlapped-regions-with-position"
    # check if output file already exists
    if [ -f "${output}".txt ]; then
        echo "Output file ${output}.txt already exists. Checking if region-only file exists..."
        # construct region-only output filename
        output_reg="${output_dir}/collapsed-filt-regions/${motif_prefix}_variant-overlapped-regions-uniq"
        if [ -f "${output_reg}".txt ]; then
            echo "Region-only output file ${output_reg}.txt already exists. Skipping..."
            continue
        else
            # generate region-only output file
            cut -f1-4 "${output}".txt | sort -k1,1V -k2,2n -u > "${output_reg}".txt
            echo "Region-only output file ${output_reg}.txt did not exist, so it has been generated."
            continue
        fi
    fi
    # create a header for the output file
    printf "chr\tstart\tend\tstrand\tdataset\tchrom\tposition\n" > "${output}.txt"
    # intersect the consensus filtered regions with the vcf files of a particular motif
    if bedtools intersect -a "${region_file}" -b "${vcf_files[@]}" -wo -names "${data_ids[@]}" | cut -f1-7 >> "${output}.txt"; then
        echo "Successfully overlapped variant positions across all subtypes with the motif ID ${motif_prefix} region file and truncated the output file."
        # generate region-only output file
        echo "Generating region-only output file..."
        output_reg="${output_dir}/collapsed-filt-regions/${motif_prefix}_variant-overlapped-regions-uniq"
        cut -f1-4 "${output}".txt | sort -k1,1V -k2,2n -u > "${output_reg}".txt
    else
        echo "Failed to run bedtools intersect. Please investigate why."
    fi
done