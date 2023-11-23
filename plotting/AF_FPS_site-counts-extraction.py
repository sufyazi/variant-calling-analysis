#!/usr/bin/env python3

####################
# import libraries #
####################
import os
import sys
import pandas as pd
from pathlib import Path

##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 2:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_site-counts-extraction.py <root_dir> <output_path>')
	sys.exit(1)
else:
	root_dir = sys.argv[1]
	output_path = sys.argv[2]

##################
# define globals #
##################

if __name__ == '__main__':
	# Find all max_af_df files in the root directory
	target_dir = Path(root_dir)
	tsv_files = target_dir.glob('**/*maximum_af*count.tsv')
	print(f'Currently processing maximum_af files in target directory: {target_dir}')
	for i, tsv in enumerate(tsv_files):
		if i == 0:
			print(f'[File count:{i+1}] First file {os.path.basename(tsv)} has been loaded.')
			# get motif ID
			motif_id = tsv.stem.replace('_maximum_af_regions_count', '')
			print(motif_id)
			# load the first file
			df_af = pd.read_csv(tsv, sep='\t')
			# rename the column region_id to max_AF_region_count
			df_af = df_af.rename(columns={'region_id': 'max_AF_region_count'})
			# find the corresponding maximum fps file in the same target dir based on the motif ID
			max_fps_file = str(next(target_dir.glob(f'**/{motif_id}_maximum_fps-scaled_regions_count.tsv')))
			# load the maximum fps file
			df_fps = pd.read_csv(max_fps_file, sep='\t')
			# rename the column region_id to region_counts
			df_fps = df_fps.rename(columns={'region_id': 'max_FPS_region_count'})
			# merge the two dataframes by sample ID
			df_final = pd.merge(df_af, df_fps, on='sample_id')
			# insert motif_id column to the merged dataframe
			df_final.insert(0, 'motif_id', motif_id)
			print(df_final.tail())
			print('Moving on to the next file...')
		else:
			print(f'[File count:{i+1}] Subsequent file {os.path.basename(tsv)} has been loaded.')
			# get motif ID
			motif_id = tsv.stem.replace('_maximum_af_regions_count', '')
			print(motif_id)
			# load the max AF file
			df_af = pd.read_csv(tsv, sep='\t')
			# rename the column region_id to max_AF_region_count
			df_af = df_af.rename(columns={'region_id': 'max_AF_region_count'})
			# find the corresponding maximum fps file in the same target dir based on the motif ID
			max_fps_file = str(next(target_dir.glob(f'**/{motif_id}_maximum_fps-scaled_regions_count.tsv')))
			# load the maximum fps file
			df_fps = pd.read_csv(max_fps_file, sep='\t')
			# rename the column region_id to region_counts
			df_fps = df_fps.rename(columns={'region_id': 'max_FPS_region_count'})
			# merge the two dataframes by sample ID
			df_merged = pd.merge(df_af, df_fps, on='sample_id')
			# insert motif_id column to the merged dataframe
			df_merged.insert(0, 'motif_id', motif_id)
			print(df_merged.head())
			print('Concatenating this merged dataframe with the previous one...')
			print('Previous dataframe:')
			print(df_final.tail())
			# concatenate the two dataframes
			df_final = pd.concat([df_final, df_merged])
			print('New dataframe:')
			print(df_final.tail())
			print('Moving on to the next file...')
	
 	# sort the final dataframe
	df_final = df_final.sort_values(by='motif_id')
 	# save the final dataframe to file
	df_final.to_csv(f'{output_path}/1360_motifs_variable_site_counts-sorted.tsv', sep='\t', index=False)
	print(f'All files have been concatenated and saved to file. Done!')