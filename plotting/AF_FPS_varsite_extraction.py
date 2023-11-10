#!/usr/bin/env python3

####################
# import libraries #
####################
import sys
import pandas as pd
from pathlib import Path

##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 2:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_varsite_extraction.py <root_dir>')
	sys.exit(1)
else:
    root_dir = sys.argv[1]

##################
# define globals #
##################
output_path = '/home/msazizan/hyperspace/gatk-workflow/plotting'

if __name__ == '__main__':
	# Find all *.tsv files in root_dir
	target_dir = Path(root_dir)
	tsv_files = target_dir.glob('*counts.tsv')
    
	for i, tsv in enumerate(tsv_files):
		if i == 0:
			# create a dataframe from the first file
			df_first = pd.read_csv(tsv, sep='\t')
			# rename a column
			df_first = df_first.rename(columns={'region_id': 'unique_sites'})
			print(f'[File count:{i}] First file {tsv} has been loaded.')
		else:
			# load the next file
			df_current = pd.read_csv(tsv, sep='\t')
   			# rename a column
			df_current = df_current.rename(columns={'region_id': 'unique_sites'})
			# concatenate the dataframes
			df_first = pd.concat([df_first, df_current]).reset_index(drop=True)
			print(f'[File count:{i}] {tsv} has been concatenated.')
	
 	# sort the final dataframe
	df_final = df_first.sort_values(by='motif_id')
 	# save the concatenated dataframe to file
	df_final.to_csv(f'{output_path}/output-data/combined_varsite_counts/1360-motifs-combined-AF-filtered-varsite-counts.tsv', sep='\t', index=False)
	print(f'All files have been concatenated and saved to file. Done!')