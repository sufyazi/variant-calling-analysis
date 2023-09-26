#!/usr/bin/env bash
# shellcheck disable=SC1091

# set up variables
tfbs_mat_dir=$1
vcfs_raw_dir=$2
data_id_list=$3
output_dir=$4

# process each dataset ID
readarray -t data_ids < "$data_id_list"

# load all files in tfbs_beds
readarray -t tfbs_matrices < <(find "$tfbs_mat_dir" -name "*brca.txt" -type f)
for matrix in "${tfbs_matrices[@]}"; do
    echo "Processing $matrix"
    matrix_basename=$(basename "$matrix" .txt)
    # grab motif ID from matrix filename
    motif_id=$(echo "$matrix_basename" | cut -d "_" -f 1-3)

    # Build the find command dynamically
    find_command="find ${vcfs_raw_dir} \( -name"
    first_name=true
    # Loop through the search terms and add them to the find command
    for id in "${data_ids[@]}"; do
        if $first_name; then
            find_command+=" ${id}_${motif_id}*.vcf"
            first_name=false
        else
            find_command+=" -o -name ${id}_${motif_id}*.vcf"
        fi
    done

    # Add the closing parenthesis for the -o expressions
    find_command+=" \) -type f"
    # Execute the find command
    readarray -t vcfs < <(eval "$find_command")
    # if no VCFs found, skip to next matrix
    if [[ -z "${vcfs[*]}" ]]; then
        echo "No VCFs found for $motif_id"
        continue
    else
        echo "Found VCFs for $motif_id"
        echo "VCFs: " "${vcfs[@]}"
        # construct noheader filename
        noheader="${output_dir}/${matrix_basename}_noheader.tmp"
        # remove header and save to new file
        tail -n +2 "$matrix" > "${noheader}"
        
        for vcf in "${vcfs[@]}"; do
            # get basename without extension
            vcf_basename=$(basename "$vcf" _qualgt10.var.flt.vcf)
            echo "VCF basename: $vcf_basename"
            if [[ ! -f "${vcf_basename}_qualgt10.var.flt_noheadlines.tmp" ]]; then
                echo "No associated files found for $vcf_basename"
                # get rid of header lines in vcf file
                grep -v "^#" "$vcf" > /data/bank/booth/msazizan/tobias-analysis-outputs/tobias-tfbs-subset-matrices-datashop/reco/v0/prod/brca-called-variant-vcfs-nohead/"${vcf_basename}_qualgt10.var.flt_noheadlines.tmp" && \
                # get rid of ## header lines (keeping the #CHROM line)
                grep -v "^##" "$vcf" > /data/bank/booth/msazizan/tobias-analysis-outputs/tobias-tfbs-subset-matrices-datashop/reco/v0/prod/brca-called-variant-vcfs-nohead/"${vcf_basename}_qualgt10.var.flt_noheadlines.txt" && \
                # generate a bedtools-compatible bed file from the vcf
                grep -v "^#" "$vcf" | awk -v OFS="\t" '{print $1, $2, $2 + 1}' > /data/bank/booth/msazizan/tobias-analysis-outputs/tobias-tfbs-subset-matrices-datashop/reco/v0/prod/brca-called-variant-vcfs-nohead/"${vcf_basename}_qualgt10.var.flt_noheadlines.bed"
            else
                echo "Found associated files for $vcf_basename. Skipping..."
                continue
            fi  
        done
        # find the tmp vcf files
        readarray -t vcf_beds < <(find /data/bank/booth/msazizan/tobias-analysis-outputs/tobias-tfbs-subset-matrices-datashop/reco/v0/prod/brca-called-variant-vcfs-nohead -name "*${motif_id}_qualgt10.var.flt_noheadlines.bed" -type f)

        if [[ -z "${vcf_beds[*]}" ]]; then
            echo "No VCF beds found for $motif_id"
            continue
        else
            echo "Found VCF beds for $motif_id"
            echo "VCFs: " "${vcf_beds[@]}"
            # construct output filename and save header
            output="${output_dir}/${motif_id}_BRCA-vcfs-filtered-matrix.txt"
            head -n 1 "$matrix" > "$output"
            
            # check if output file exists
            if [[ -f "$output" && $(wc -l < "$output") -gt 1 ]]; then
                echo "Complete output file $output exists. Skipping..."
                continue
            else
                echo "Output file $output is incomplete. Proceeding..."
                # preprocess vcf then pipe to bedtools intersect to find the TFBSs that overlap with the VCF
                if bedtools intersect -a "${noheader}" -b "${vcf_beds[@]}" -wa >> "$output"; then
                    echo "Successfully intersected $matrix with all dataset VCFs"
                else
                    echo "Failed to intersect $matrix with ${#vcf_beds[@]} VCF bed files of $motif_id"
                fi
            fi
        fi
        rm "${noheader}"
    fi
done

