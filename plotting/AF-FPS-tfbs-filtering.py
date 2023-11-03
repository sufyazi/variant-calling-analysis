#!/#!/usr/bin/env python3

# import libraries
import os
import sys
import fnmatch
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt

from functools import reduce
from sklearn.preprocessing import MinMaxScaler

# define a generator function for file paths
def path_generator(root_dir):
    for filename in os.listdir(root_dir):
        if filename.endswith('.txt'):
            yield os.path.join(root_dir, filename)
            
# create a vcf load function for the query vcfs
def load_vcf(vcf_path):
    # load up the vcf file with indels and multiallelic sites split into separate rows
    df_vcf = pd.read_csv(vcf_path, sep="\t")
    # rename columns in the dataframe
    df_vcf = df_vcf.rename(columns={"#[1]CHROM": "Chromosome", "[2]POS": "Start", "[3]REF": "ref_allele", "[4]ALT": "alt_allele", "[5]AF": "AF"})
    # add a column next to the "start" column called "end" with the same value as the "start" column
    df_vcf.insert(2, "End", df_vcf["Start"])
    return df_vcf

# create a function to find files with a specific pattern in their filenames
def find_files(path, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

#############
# load args #
#############

fps_path = sys.argv[1] # path to where the filtered matrices of footprinting scores are stored

af_path = sys.argv[2] # path to where the allelic frequency variant data are stored

# read in file containing dataset IDs with each line an element of a new list
dataset_ids = [line.rstrip('\n') for line in open(sys.argv[2])]
with open(sys.argv[3])) as file:
    dataset_ids = [line.rstrip('\n') for line in file]
print(f"Dataset IDs to be processed: {dataset_ids}")

output_path = sys.argv[4] # path to where the output files will be stored

#############
# load data #
############# 
      
for file in path_generator(fps_path):
    suffix = "_BRCA-subtype-vcf-filtered-matrix.txt"
    motif_id = os.path.basename(file).replace(suffix, '')
	print(f"Processing filtered TFBS matrix of {motif_id}...")
    # load the data
    df_fps = pd.read_csv(file, sep="\t")
    # drop the column "TFBS_strand" and "TFBS_score"
    df_fps = df_fps.drop(columns=["TFBS_strand", "TFBS_score"])
    # rename columns in the dataframe
    df_fps = df_fps.rename(columns={"TFBS_chr": "Chromosome", "TFBS_start": "Start", "TFBS_end": "End", "2GAMBDQ_Normal-like_score": "2GAMBDQ_Norm_fps"})
    # for all column names that end with the string 'score', replace the string with 'fps'
    df_fps = df_fps.rename(columns=lambda x: x.replace('score', 'fps') if x.endswith('score') else x)

    # load associated vcf files of the motif name for each dataset ID
    # first search for the associated vcf files based on the motif name
    vcf_paths = [os.path.join(af_path, dataset_id, f"{dataset_id}_{motif_name}_2.00_AF-per-site-with-indels.txt") for dataset_id in dataset_ids]

# create a list of tuples containing the dataset ID and the path to the corresponding vcf file

paths = [
    "/home/msazizan/hyperspace/gatk-workflow/plotting/raw_inputs/af/ANAB5F7/ANAB5F7_TFDP1_M08108_2.00_AF-per-site-with-indels.txt",
    "/home/msazizan/hyperspace/gatk-workflow/plotting/raw_inputs/af/98JKPD8/98JKPD8_TFDP1_M08108_2.00_AF-per-site-with-indels.txt",
    "/home/msazizan/hyperspace/gatk-workflow/plotting/raw_inputs/af/PU24GB8/PU24GB8_TFDP1_M08108_2.00_AF-per-site-with-indels.txt",
    "/home/msazizan/hyperspace/gatk-workflow/plotting/raw_inputs/af/S6R691V/S6R691V_TFDP1_M08108_2.00_AF-per-site-with-indels.txt",
    "/home/msazizan/hyperspace/gatk-workflow/plotting/raw_inputs/af/2GAMBDQ/2GAMBDQ_TFDP1_M08108_2.00_AF-per-site-with-indels.txt"
]

# create a list of IDs
ids = [ "ANAB5F7_basal", "98JKPD8_lumA", "PU24GB8_lumB", "S6R691V_her2", "2GAMBDQ_norm"]

# create a pair dictionary
path_id_dict = {id: load_vcf(path) for id in ids for path in paths if id.split("_")[0] in path}

print(path_id_dict)

gr_fpscore = pr.PyRanges(df_fpscore)
print(type(gr_fpscore))
gr_fpscore

count = 0
for key, val in grs.items():  
    if count == 0:
        overlap = gr_fpscore.join(val, how='left', suffix=f"_{key}_varsite_pos", preserve_order = True)
    else:
        overlap = filtered_gr.join(val, how='left', suffix=f"_{key}_varsite_pos", preserve_order = True)
    
    # drop the column "End" column
    overlap = overlap.drop([f"End_{key}_varsite_pos"])

    # cluster the pyRanges object by genomic range
    overlap = overlap.cluster(slack=-1)

    # cast back into a dataframe and filter by the AF column's max value (per cluster); this returns a filtered dataframe
    filtered_df = overlap.df.loc[overlap.df.groupby('Cluster')['AF'].idxmax()]

    # then rename the columns
    filtered_df = filtered_df.rename(columns={f"Start_{key}_varsite_pos": f"{key}_varsite_pos", "ref_allele": f"{key}_ref_allele", "alt_allele": f"{key}_alt_allele", "AF": f"{key}_AF"})
    
    # replace all the -1 values in column 'Start_varsites', 'ref_allele' and 'alt_allele', and AF with 0
    # Define a dictionary mapping column names to values to replace
    replace_dict = {f"{key}_varsite_pos": {-1: None}, f"{key}_ref_allele": {str(-1): None}, f"{key}_alt_allele": {str(-1): None}, f"{key}_AF": {-1: 0}}
    filtered_df = filtered_df.replace(replace_dict)

    # drop cluster column
    filtered_df = filtered_df.drop(columns=["Cluster"])

    # cast back into pyrange object
    filtered_gr = pr.PyRanges(filtered_df)

    # increment count
    count += 1

final_df = filtered_gr.df

# create a column called 'region_id'
final_df["region_id"] = final_df["Chromosome"].astype(str) + ":" + final_df["Start"].astype(str) + "-" + final_df["End"].astype(str)

# for all column name ending with the string '_fps', split the string, take the second element, change the first letter in the string to lowercase, and reconstruct the original string with the new first letter
final_df = final_df.rename(columns=lambda x: x.split('_')[0] + '_' + x.split('_')[1][0].lower() + x.split('_')[1][1:] + '_fps' if x.endswith('_fps') else x)

# save to file
final_df.to_csv("TFDP1_M08108_2.00_fpscore-af-varsites-combined-matrix-wide.tsv", sep="\t", index=False, na_rep='NULL')

# filtered_gr.df.shape
final_df.iloc[19:30]