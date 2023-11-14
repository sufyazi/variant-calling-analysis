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
	tsv_files = target_dir.glob('*-variance.tsv')
    
	for i, tsv in enumerate(tsv_files):
		if i == 0:
			# get motif ID
			motif_id = tsv.stem.replace('_afps-variance', '')
			# load the first file
			df_first = pd.read_csv(tsv, sep='\t')
			# find the number of unique sites
			counts = df_first['region_id'].nunique()
			# create an empty df
			final_df = pd.DataFrame(columns=['motif_id', 'filtered_variable_site_counts'])
			# append the motif ID and counts to the empty df
			final_df = final_df.append({'motif_id': motif_id, 'filtered_variable_site_counts': counts}, ignore_index=True)

			print(f'[File count:{i}] First file {tsv} has been loaded.')
		else:
			# load the next file
			df_current = pd.read_csv(tsv, sep='\t')
			# get motif ID
			motif_id = tsv.stem.replace('_afps-variance', '')
			# find the number of unique sites
			counts = df_current['region_id'].nunique()
			# append the motif ID and counts to the empty df
			final_df = final_df.append({'motif_id': motif_id, 'filtered_variable_site_counts': counts}, ignore_index=True)
			print(f'[File count:{i}] {tsv} has been loaded.')
	
 	# sort the final dataframe
	final_df = final_df.sort_values(by='filtered_variable_site_counts', ascending=False)
 	# save the final dataframe to file
	final_df.to_csv(f'{output_path}/output-data/combined_varsite_counts/1360-motifs-sorted_variable-TFBS-counts-filtered-AF.tsv', sep='\t', index=False)
	print(f'All files have been concatenated and saved to file. Done!')