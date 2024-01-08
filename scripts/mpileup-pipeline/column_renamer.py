#!/usr/bin/env python3

import sys
import polars as pl

# Usage: python3 column_renamer.py <path_to_vcf_nohead_file> <path_to_mapping_file> <output_path>
# Check for the correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python3 column_renamer.py <path_to_vcf_nohead_file> <path_to_mapping_file> <output_filename_with_path>")
    sys.exit(1)

# Get the path of the output noheader vcf file (.txt)
path = sys.argv[1]

# Get the path to the mapping file
mapping_file = sys.argv[2]

# Get the output path
out_path = sys.argv[3]

def rename_columns(path, mapping_dict):
    """
    Rename the column names based on the mapping dictionary
    """
    df = pl.read_csv(path, separator="\t")
    df.rename(mapping_dict)
    
    return df

if __name__ == "__main__":
    mapping_dict = {}
    with open(mapping_file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            mapping_dict[line[0]] = line[1]
    
    df = rename_columns(path, mapping_dict)
    df.write_csv(out_path, separator="\t")
    