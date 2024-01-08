#!/usr/bin/env python3

####################
# import libraries #
####################
import os
import sys
import pandas as pd
import scipy.stats as stats
from pathlib import Path

##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 2:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_extract-discrete-outlier-counts.py <root_dir> <output_path>')
	sys.exit(1)
else:
	root_dir = sys.argv[1] # where the motif ID table directories are
	output_path = sys.argv[2]

##################
# define globals #
##################

if __name__ == '__main__':
	# Find all motif directories in the root directory
	target_dir = Path(root_dir)
	print(f'Currently processing motifs in target directory: {target_dir}')
	motif_dirs = target_dir.glob('*/')
	
	# initialize a dict to store the correlation sets
	motif_corr_dict = {}

	for i, m in enumerate(motif_dirs):
		motif_id = os.path.basename(m)
		print(f'Currently processing motif {motif_id} in target directory')
		# find the filtered starting df
		hiaf_file = str(next(m.glob(f'*_high_AF_regs_abv_IQR_threshold_sorted_by_FPS_var_table.tsv')))
		print(f'[File count:{i+1}] HI-AF filtered region file of {motif_id} has been found: {os.path.basename(hiaf_file)}')
		# load the file
		hiaf_df = pd.read_csv(hiaf_file, sep='\t')
		
		# create a dictionary to store motif correlation values
		subtype_correlations = {}

		# calculate the Spearman correlation
		for subtype in hiaf_df['sample_id'].unique():
			corr_df = hiaf_df[hiaf_df['sample_id'] == subtype].copy()
			correlation, pvalue = stats.spearmanr(corr_df['FPS_scaled'], corr_df['AF'])
			print(subtype, correlation, pvalue)
			# create a tuple
			corr_tuple = (correlation, pvalue)
			# append to subtype dictionary
			subtype_correlations[subtype] = corr_tuple

		# append to the motif dictionary
		motif_corr_dict[motif_id] = subtype_correlations

	# initialize a dataframe to store the correlation values
	df_final = pd.DataFrame()
	print(df_final.head())
	print(f'Now printing the correlation dictionary: {motif_corr_dict}')
	# now for each motif, create a dataframe of the correlation values and pvalues
	for i, motif_id in enumerate(motif_corr_dict.keys()):
		print(f'Currently processing motif {motif_id} in the correlation dictionary: motif [{i+1}]')
		for j, (subtype, corr) in enumerate(motif_corr_dict[motif_id].items()):
			print(f'subtype {j+1}: {subtype}, {corr}')
			# create a dataframe from the tuple
			df = pd.DataFrame([corr], columns=['correlation', 'pvalue'])
			print(df.head())
			# add a column for the motif ID
			df['motif_id'] = motif_id
			# add a column for the subtype
			df['subtype'] = subtype
			if j == 0:
				# create a dataframe from the dictionary
				df_temp = df
				j += 1
			else:
				# concat the dataframes
				df_temp = pd.concat([df_temp, df], ignore_index=True)
				print(f'now printing df_temp when j > 0: {df_temp.head()}')
				j += 1
		# concat to the final dataframe
		df_final = pd.concat([df_final, df_temp], ignore_index=True)
		
	# clean up the dataframe
	df_final = df_final[['motif_id', 'subtype', 'correlation', 'pvalue']]
	# sort by ascending pvalue
	df_final = df_final.sort_values(by=['pvalue'])
	# save to file
	df_final.to_csv(f'{output_path}/AF-FPS_region_Spearman-corr-by-subtype-sorted-by-pval.tsv', sep='\t', index=False)