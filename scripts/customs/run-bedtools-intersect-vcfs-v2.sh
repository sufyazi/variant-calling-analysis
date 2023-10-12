#!/usr/bin/env bash
# shellcheck disable=SC1091

# check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <tfbs_matrix_dir> <vcfs_raw_dir> <data_id_list> <output_dir>"
    exit 1
fi

module load miniconda3/py38_4.8.3
conda activate /home/users/ntu/suffiazi/apps/mambaforge/envs/gatk-snv

# set up variables
tfbs_matrx_dir=$1
vcfs_raw_dir=$2
data_id_list=$3
output_dir=$4

###### HARCODED PATHS !! Change only if needed ##########
HEADERLESS="/scratch/users/ntu/suffiazi/outputs/brca-called-variants-diffmode-noheader"
BEDS="/scratch/users/ntu/suffiazi/inputs/BRCA-diffmode_tfbs_beds"
####################################################

# process each dataset ID
readarray -t data_ids < "$data_id_list"

# load motif ID from subdir name
readarray -t tfbs_motifs < <(find "$tfbs_matrx_dir" -mindepth 1 -type d)

#initialize count
count=0

for direc in "${tfbs_motifs[@]}"; do
    # increment count
    count=$((count + 1))
    motif_id=$(basename "$direc")
    matrix=$(find "$direc" -name "${motif_id}*matrix.txt" -type f)
    echo "Searching for $matrix..."
    
    # check for valid value in $matrix
    if [[ -z "$matrix" ]]; then
        echo "Matrix variable is empty. Please check your input file. Skipping..."
        continue
    else
        echo "Found a matrix file associated with $motif_id: $matrix"
        echo "Proceeding..."

        # Build the find command dynamically to search for associated VCF files across all dataset IDs
        find_command="find ${vcfs_raw_dir} \( -name"
        first_name=true
        # Loop through the search terms and add them to the find command
        for id in "${data_ids[@]}"; do
            if $first_name; then
                find_command+=" ${id}_${motif_id}*VAF.vcf"
                first_name=false
            else
                find_command+=" -o -name ${id}_${motif_id}*VAF.vcf"
            fi
        done
    
        # Add the closing parenthesis for the -o expressions
        find_command+=" \) -type f"
    
        # Execute the find command to find all associated VCF files
        readarray -t vcfs < <(eval "$find_command")

        # if no VCFs found, skip to next motif
        if [[ -z "${vcfs[*]}" ]]; then
            echo "No VCF files of all tested dataset IDs are found for $motif_id. Skipping..."
            continue
        else
            echo "Found VCFs for ${motif_id}." 
            echo "VCFs: " "${vcfs[@]}"
            echo "Motif count: $count"
        
            # generate headerless vcf file counterparts for the datasets
            for vcf in "${vcfs[@]}"; do
                # get basename without extension
                vcf_basename=$(basename "$vcf" .vcf)
                if [[ ! -f "${HEADERLESS}"/"${vcf_basename}_noheadlines.txt" ]]; then
                    echo "No headerless file of $vcf_basename was found. Creating headerless file..."
                    # get rid of ## header lines (keeping the #CHROM line)
                    if grep -v "^##" "$vcf" > "${HEADERLESS}/${vcf_basename}_noheadlines.txt"; then
                        echo "Created headerless file for $vcf_basename"
                    else
                        echo "Failed to create headerless file for $vcf_basename. Check logs"
                        exit 1
                    fi
                else
                    echo "Found associated files for $vcf_basename. Skipping..."
                    continue
                fi  
            done

            # construct masking file name
            mask="${output_dir}/${motif_id}_BRCA-subtype-vcf-filtered-mask"
            # construct output file name
            output="${output_dir}/${motif_id}_BRCA-subtype-vcf-filtered-matrix"
            
            echo "Checking output file existence..."

            # check if output file already exists
            if [[ -f "${output}.txt" ]]; then
                echo "Output filtered matrix file for $motif_id already exists. Skipping..."
                continue
            else
                echo "Output matrix for $motif_id not found. Proceeding with filtering..."
                
                # find bed file associated with the matrix file
                matrix_bed=$(find "${BEDS}" -name "${motif_id}_*.bed" -type f)

                if [ -n "$matrix_bed" ]; then
                    echo "Matrix bed file found: $matrix_bed"
                    # run bedtools intersect to find the TFBSs of the motif that overlap with the called variant positions in the VCFs
                    if bedtools intersect -a "${matrix_bed}" -b "${vcfs[@]}" -u -header > "${mask}.bed"; then
                        echo "Successfully intersected $matrix_bed with all dataset VCFs to create a masking bed file. Proceeding..."
                        # save the header line to a temporary file
                        head -n 1 "${matrix}" > "${output}.tmp"
                        # bedtools intersect the original matrix with the masking bed file
                        if tail -n +2 "${matrix}" | bedtools intersect -a - -b "${mask}.bed" -u >> "${output}.tmp"; then
                            echo "Successfully intersected $matrix with masking bed file tocreate a filtered matrix."
                            # rename the temporary file
                            mv "${output}.tmp" "${output}.txt"
                            # remove the masking bed file
                            rm "${mask}.bed"
                        else
                            echo "Failed to intersect $matrix with masking bed file. Check logs."
                            exit 1
                        fi
                    else
                        echo "Failed to intersect $matrix_bed with all dataset VCFs. Check logs."
                        exit 1
                    fi
                else
                    echo "Bed file for $motif_id is not found. Check again. Exiting..."
                    exit 1
                fi
            fi
        fi
    fi
done
