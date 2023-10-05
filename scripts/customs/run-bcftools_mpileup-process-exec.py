#!/usr/bin/env python3

import concurrent.futures
import subprocess
import os
import sys

# set up variables
bams = sys.argv[1]
outfile = sys.argv[2]
data_id = sys.argv[3]
motif_id = sys.argv[4]

# List of file paths
file_paths = ['/path/to/file1.txt', '/path/to/file2.txt', '/path/to/file3.txt']

# Bash function to execute
def bash_exec(bed_path, bams, outfile, data_id, motif_id):
    try:
        # Define bash command
        # bash_command = f'bcftools mpileup -Ou -f /scratch/users/ntu/suffiazi/inputs/references/gatk4/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -T {bed_path} -b {bams} | bcftools call -Ou -mv | bcftools filter -i "QUAL>10" > {outfile}/{data_id}_{motif_id}_qualgt10.var.flt.vcf'
        
        bash_command = f'bcftools --help'

        # Execute the bash command using subprocess
        subprocess.run(bash_command, shell=True, check=True, text=True)
   
    except subprocess.CalledProcessError as e:
        return f"Error executing bash command: {e}\n"

# Create a ProcessPoolExecutor with the desired number of processes (adjust as needed)
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Submit the tasks for each file to the executor and collect the results
    results = list(executor.map(bash_exec, file_paths))

# Print or process the results as needed
for result in results:
    print(result)
