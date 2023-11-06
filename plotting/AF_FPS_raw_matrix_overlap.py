#!/usr/bin/env python3

####################
# import libraries #
####################

import os
import sys
import fnmatch
import pandas as pd
import pyranges as pr
import concurrent.futures
import itertools

####################
# define functions #
####################

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
    for root, _, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

# function to overlap pyranges objects two at a time
def pyrange_obj_overlap(gr_fpscore, grs_vcf_dict):
    count = 0
    for key, val in grs_vcf_dict.items():  
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
      
    return filtered_gr

# define concurrent function to process multiple files at once
def process_file(file, af_path, dataset_ids, output_path):
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
    vcf_paths = find_files(af_path, f"*{motif_id}*.txt")
        
    if len(vcf_paths) != len(dataset_ids):
        print(f"ERROR: Number of vcf files ({len(vcf_paths)}) does not match the number of dataset IDs ({len(dataset_ids)})!")
        print("Printing both lists...")
        print(f"vcf_paths: {vcf_paths}")
        print(f"dataset_ids: {dataset_ids}")
        print("Exiting prematurely...")
        sys.exit(1)
        
    # create a dataset ID:af dataframe dictionary
    dataset_af_dict = {dataset: load_vcf(path) for dataset, path in zip(dataset_ids, vcf_paths)}
    print(dataset_af_dict)
        
    # create a pyranges object for the filtered TFBS footprint matrix
    gr_fpscore = pr.PyRanges(df_fps)
        
    # load up vcf dfs into pyranges 
    grs = {}
    for name, vcf in dataset_af_dict.items():
        gr_vcf = pr.PyRanges(vcf)
        grs[name] = gr_vcf

    target_gr = pyrange_obj_overlap(gr_fpscore, grs)
    target_df = target_gr.df

    # create a column called 'region_id'
    target_df["region_id"] = target_df["Chromosome"].astype(str) + ":" + target_df["Start"].astype(str) + "-" + target_df["End"].astype(str)

    # for all column name ending with the string '_fps', split the string, take the second element, change the first letter in the string to lowercase, and reconstruct the original string with the new first letter
    target_df = target_df.rename(columns=lambda x: x.split('_')[0] + '_' + x.split('_')[1][0].lower() + x.split('_')[1][1:] + '_fps' if x.endswith('_fps') else x)

    # construct output filename
    outfile = os.path.join(output_path, f"{motif_id}_fpscore-af-varsites-combined-matrix-wide.tsv")
        
    # save to file
    target_df.to_csv(outfile, sep="\t", index=False, na_rep='NULL')

    # print the dimensions of the dataframe
    print(f"Shape of the current motif ID ({motif_id}): {target_df.shape}")
    print(f"Output file for {motif_id} has been generated!")
    
##################
# load arguments #
##################

fps_path = sys.argv[1] # path to where the filtered matrices of footprinting scores are stored

af_path = sys.argv[2] # path to where the allelic frequency variant data are stored

# read in file containing dataset IDs with each line an element of a new list
with open(sys.argv[3]) as file:
    dataset_ids = [line.rstrip('\n') for line in file]
print(f"Dataset IDs to be processed: {dataset_ids}")

output_path = sys.argv[4] # path to where the output files will be stored

#############
# load data #
############# 

if __name__ == "__main__":
    # run concurrent processes
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(process_file, path_generator(fps_path), itertools.repeat(af_path), itertools.repeat(dataset_ids), itertools.repeat(output_path))  

    print ("All footprint matrices have been processed!")