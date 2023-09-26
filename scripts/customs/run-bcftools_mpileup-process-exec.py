#!/usr/bin/env python3

import concurrent.futures
import subprocess
import os
import sys

# set up variables
input_files = sys.argv[1]
output_path = sys.argv[2]
process_count = sys.argv[3]

# load up the list of files to process by reading the file paths from variable
file_paths = [ for file in os.listdir(input_path) if file.endswith('brca.bed')]






# List of file paths
file_paths = ['/path/to/file1.txt', '/path/to/file2.txt', '/path/to/file3.txt']

# Function to execute the bash command on a file and return the output
def execute_bash_command(file_path):
    try:
        # Define your bash command
        bash_command = f'your_bash_command_here {file_path}'
        
        # Execute the bash command using subprocess
        result = subprocess.run(bash_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Return the command output as a string
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error executing command: {e}\n"

# Create a ProcessPoolExecutor with the desired number of processes (adjust as needed)
with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
    # Submit the tasks for each file to the executor and collect the results
    results = list(executor.map(execute_bash_command, file_paths))

# Print or process the results as needed
for result in results:
    print(result)
