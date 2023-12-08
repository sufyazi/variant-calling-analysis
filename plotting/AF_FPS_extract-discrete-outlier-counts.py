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


# create a 2x2 contingency matrix
def create_contingency_table(df1, df2, df3, df4):
    # Get the unique counts
    hi_af_hi_fps = df1['region_id'].nunique()
    print(f'There are {hi_af_hi_fps} unique region_ids in the HI AF HI FPS table')
    hi_af_lo_fps = df2['region_id'].nunique()
    print(f'There are {hi_af_lo_fps} unique region_ids in the HI AF LO FPS table')
    lo_af_hi_fps = df3['region_id'].nunique()
    print(f'There are {lo_af_hi_fps} unique region_ids in the LO AF HI FPS table')
    lo_af_lo_fps = df4['region_id'].nunique()
    print(f'There are {lo_af_lo_fps} unique region_ids in the LO AF LO FPS table')

    # Define the data
    data = {'HI-FPS': [hi_af_hi_fps, lo_af_hi_fps], 'LO-FPS': [hi_af_lo_fps, lo_af_lo_fps]}

    # Create the DataFrame and set the index
    df = pd.DataFrame(data, index=['HI-AF', 'LO-AF'])

    return df

##################
# define globals #
##################

if __name__ == '__main__':
    
    # initialize a dictionary to store the significantly mutated motifs
	significant_motifs = {}

	# Find all motif directories in the root directory
	target_dir = Path(root_dir)
	print(f'Currently processing motifs in target directory: {target_dir}')
	motif_dirs = target_dir.glob('*/')
	for i, m in enumerate(motif_dirs):
		motif_id = os.path.basename(m)
		print(f'Currently processing motif {motif_id} in target directory')
		# find and load the HI-AF_HI-FPS df
		hiaf_hifps_file = str(next(m.glob(f'*_HI_AF_regs_GT_FPS-mean_sorted_by_FPS_var_table.tsv')))
		print(f'[File count:{i+1}] HI-AF_HI-FPS region file of {motif_id} has been found: {os.path.basename(hiaf_hifps_file)}')
		# load the file
		hiaf_hifps_df = pd.read_csv(hiaf_hifps_file, sep='\t')
		# find and load the HI-AF_LO-FPS df
		hiaf_lofps_file = str(next(m.glob(f'*_HI_AF_regs_LT_FPS-mean_sorted_by_FPS_var_table.tsv')))
		print(f'[File count:{i+1}] HI-AF_LO-FPS region file of {motif_id} has been found: {os.path.basename(hiaf_lofps_file)}')
	 	# load the file
		hiaf_lofps_df = pd.read_csv(hiaf_lofps_file, sep='\t')
		# find and load the LO-AF_HI-FPS df
		loaf_hifps_file = str(next(m.glob(f'*_LO_AF_regs_GT_FPS-mean_sorted_by_FPS_var_table.tsv')))
		print(f'[File count:{i+1}] LO-AF_HI-FPS region file of {motif_id} has been found: {os.path.basename(loaf_hifps_file)}')
	 	# load the file
		loaf_hifps_df = pd.read_csv(loaf_hifps_file, sep='\t')
		# find and load the LO-AF_LO-FPS df
		loaf_lofps_file = str(next(m.glob(f'*_LO_AF_regs_LT_FPS-mean_sorted_by_FPS_var_table.tsv')))
		print(f'[File count:{i+1}] LO-AF_LO-FPS region file of {motif_id} has been found: {os.path.basename(loaf_lofps_file)}')
	 	# load the file
		loaf_lofps_df = pd.read_csv(loaf_lofps_file, sep='\t')
  
		# create a contingency table
		df_final = create_contingency_table(hiaf_hifps_df, hiaf_lofps_df, loaf_hifps_df, loaf_lofps_df)
		# save to file
		df_final.to_csv(f'{output_path}/{motif_id}_AF-FPS_region_contingency_table.tsv', sep='\t', index=True)

		# calculate the fisher exact test
		oddsratio, pvalue = stats.fisher_exact(df_final)
		print(f'Fisher exact test for {motif_id}: {pvalue}')

		if pvalue < 0.05:
			print(f'Fisher exact test for {motif_id} is significant. Adding to the dictionary of significantly mutated motifs...')
			significant_motifs[motif_id] = pvalue
		else:
			print(f'Fisher exact test for {motif_id} is not significant. Skipping...')

	# save the dictionary of significant motifs to file
	print(f'Finished processing all motifs in target directory: {target_dir}')
	print(f'Saving the dictionary of significant motifs to file...')
	with open(f'{output_path}/AF-FPS_significant_motifs_dictionary.tsv', 'w') as f:
		for key in significant_motifs.keys():
			f.write("%s\t%s\n"%(key,significant_motifs[key]))
	print('Done!')