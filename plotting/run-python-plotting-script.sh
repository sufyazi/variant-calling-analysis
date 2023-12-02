#!/usr/bin/env bash
# shellcheck disable=SC1091

input_dir=$1
output_dir=$2

# traverse the input directory and create a list of all the input files ending with .tsv
# then, pass this list to the python script

# find and put into an array all the files ending with .tsv
readarray -t files < <(find "${input_dir}" -name "*.tsv" -type f | sort)
# Get the length of the array
array_length="${#files[@]}"
echo "${array_length} files have been found. Splitting into batches..."

# split the array into 10 batches of 136 files
batch=136

# batch counter
counter=0

# Use a loop to iterate over the array
for ((i=0; i<array_length; i+=batch)); do
    # Increment the counter
    ((counter++))
    echo "Currently on batch ${counter}"
    # Slice the array to get a part
    part=("${files[@]:i:batch}")
    # Print or use the part
    echo "${#part[@]}" " files in this batch"
    echo "${part[0]}" ": first file in this batch"
    echo "${part[-1]}" ": last file in this batch"

    # create a directory for this batch
    mkdir -p /home/users/ntu/suffiazi/scratch/outputs/tmp/input_batch-mean/batch-"${counter}"
    batch_dir="/home/users/ntu/suffiazi/scratch/outputs/tmp/input_batch-mean/batch-${counter}"

    # create output directory for this batch
    mkdir -p "${output_dir}"/batch-"${counter}"
    outbatch_dir="${output_dir}"/batch-"${counter}"
    
    # rsync the files in this batch to a directory
    rsync -avzHP "${part[@]}" "${batch_dir}"
    echo "Created directory ${batch_dir} and rsynced files to it"
    echo "Now submitting job for this batch..."

    # run the script
    if qsub -v INPUTDIR="${batch_dir}",OUTPUTDIR="${outbatch_dir}" /home/users/ntu/suffiazi/scripts/gatk-workflow-scripts/plotting/run-python-plotting-script.pbs; then
        echo "Submitted job for input file batch ${counter} via qsub for plotting."
        echo "--------------------"
    else
        echo "ERROR: Failed to submit job for input file batch ${counter} due to an error. Please see the log file for more details."
    fi
done




