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
	print(f'USAGE: python3 AF_FPS_common-max-region-extraction.py <root_dir> <output_path>')
	sys.exit(1)
else:
	root_dir = sys.argv[1] # where the motif ID directories are
	output_path = sys.argv[2]

##################
# define globals #
##################

if __name__ == '__main__':
	# Find all motif directories in the root directory
	target_dir = Path(root_dir)
	print(f'Currently processing motifs in target directory: {target_dir}')
	motif_dirs = target_dir.glob('*/')
	for i, m in enumerate(motif_dirs):
		motif_id = os.path.basename(m)
		print(f'Currently processing motif {motif_id} in target directory')
		# find and load the max_af_regions_count.tsv file
		max_af_file = str(next(m.glob(f'*_max_af_region-ids_unique.tsv')))
		print(f'[File count:{i+1}] Max AF region file of {motif_id} has been found: {os.path.basename(max_af_file)}')
		# load the first file
		df_af = pd.read_csv(max_af_file, sep='\t')
		# find the corresponding maximum fps file in the same target dir based on the motif ID
		max_fps_file = str(next(m.glob(f'*_max_fps-scaled_region-ids_unique.tsv')))
		print(f'[File count:{i+1}] Max FPS region file of {motif_id} has been found: {os.path.basename(max_fps_file)}')
	# 	# load the maximum fps file
		df_fps = pd.read_csv(max_fps_file, sep='\t')
		print('Successfully loaded both files. Merging...')
		# merge the two dataframes to get common rows
		df_final = df_af.merge(df_fps, on=['region_id', 'sample_id'], how='inner')
		# insert motif_id column to the merged dataframe
		df_final.insert(0, 'motif_id', motif_id)
		print(df_final.tail())
		# save to file
		df_final.to_csv(f'{output_path}/{motif_id}_common_max_AF-FPS_regions.tsv', sep='\t', index=False)
	
 	# concatenate all files generated
	print('Concatenating all merged files...')
	# find all files in the output path
	output_path = Path(output_path)
	files = output_path.glob('*_common_max_AF-FPS_regions.tsv')
	# loop through all files and concatenate
	df_concat = pd.concat([pd.read_csv(f, sep='\t') for f in files])
	print(f'All files have been concatenated and saved to file. Saving to file...')
	df_concat.to_csv(f'/data5/msazizan/plotting_data/IQR-region-counts/All-motifs_common_max_AF-FPS_regions.tsv', sep='\t', index=False)
	print('Done!')