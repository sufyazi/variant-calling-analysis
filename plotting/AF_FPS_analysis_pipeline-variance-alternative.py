#!/usr/bin/env python3

####################
# import libraries #
####################

import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import concurrent.futures as cf

from pathlib import Path
from natsort import index_natsorted
from sklearn.preprocessing import MinMaxScaler

####################
# define functions #
####################
def plot_stacked_barplot(longdf, motif_id, output_path, rotate_xticks=False, xticks_fontsize=5, output_subfolder='unfilt-barplots', plot_suffix='unfilt'):
	try:
		# Get a list of unique 'sample_id' values
		sample_ids = longdf['sample_id'].unique()
		# plot the sorted stacked bar plot
		plt.figure(figsize=(12,6), dpi=300)

		# Initialize a zero array for the 'bottom' parameter of the bar plot
		bottom = np.zeros(len(longdf['region_id'].unique()))

		# For each 'sample_id'
		for i, sample_id in enumerate(sample_ids):
    		# Get the 'AF' values for this 'sample_id'
			data = longdf[longdf['sample_id'] == sample_id]

    		# Create a bar plot for this 'sample_id', stacked on top of the previous one
			sns.barplot(data=data, x='region_id', y='AF', bottom=bottom, color=sns.color_palette()[i])

    		# Update 'bottom' for the next 'sample_id'
			bottom += data['AF'].values

		# Create a patch for each 'sample_id'
		patches = [mpatches.Patch(color=sns.color_palette()[i], label=sample_id) for i, sample_id in enumerate(sample_ids)]

		# Add the legend to the plot
		plt.legend(handles=patches)
	
		if rotate_xticks:
			plt.xticks(rotation=90, fontsize=xticks_fontsize)
		else:
			plt.xticks([])
    
		plt.ylabel('Cumulative allelic frequency (AF)', fontsize=12)
		plt.xlabel(f'{motif_id} motif sites with called variants ({longdf['region_id'].nunique()})', fontsize=12)
	
		output = f'{output_path}/graphs/{output_subfolder}/{motif_id}_afps-barplot-{plot_suffix}.png'
		plt.savefig(output, dpi=300, bbox_inches='tight')
		plt.close()
		return True
	except Exception as e:
		print(f'ERROR: {e}')
		return False


def variance_calc_alt_df(filepath, output_path):
	# extract motif id from filename
	motif_id = os.path.basename(filepath).replace('_fpscore-af-varsites-combined-matrix-wide.tsv', '')
	# print message
	print(f'Processing {motif_id}...')

	# load data file
	matrix_df = pd.read_csv(filepath, sep='\t')
	# filter df
	fps_df = matrix_df.filter(regex='_fps$|_id$').copy()
	# calculate variance of fps values across samples per region_id and add to a new column called 'fps_var'
	fps_df = fps_df.set_index('region_id')
	fps_df['fps_var'] = fps_df.var(axis=1)
	# subset the dataframe to just the 'fps_var' column and reset the index
	fps_var_df = fps_df[['fps_var']].reset_index()
	# sort the df by fps_var in descending order
	fps_var_df = fps_var_df.sort_values(by='fps_var', ascending=False)
	fps_var_df = fps_var_df.set_index('region_id')
	
	#### process AF data ####
	af_df = matrix_df.filter(regex='_AF$|_id$').copy()
	af_df = af_df.set_index('region_id')
	# convert to long format
	af_df_long = af_df.melt(id_vars=["region_id"], var_name="variable", value_name="value")
	# split the variable column into sample_id and type columns using reverse split string method, which returns a dataframe of columns based on the number of splits (n=x); this can directly be assigned to new columns in the original dataframe
	af_df_long[['sample_id', 'type']] = af_df_long['variable'].str.rsplit('_', n=1, expand=True)
	# drop the redundant 'variable' column
	af_df_long = af_df_long.drop(columns=["variable"])
	# now pivot the dataframe to create new columns based on the type column
	af_df_lpv = af_df_long.pivot(index=['region_id', 'sample_id'], columns='type', values='value').reset_index()

	# remove the index name
	af_df_lpv = af_df_lpv.rename_axis(None, axis=1)

	# sort the dataframe by region_id naturally
	af_df_lpv = af_df_lpv.reindex(index=index_natsorted(af_df_lpv['region_id']))
	af_df_lpv = af_df_lpv.reset_index(drop=True)
	
	# merge af_df_sorted_lpv and fps_var_df on region_id index
	merged_df = pd.merge(fps_var_df, af_df_lpv, right_on='region_id', left_index=True, how='outer')
	# reposition region_id column to the first column, sample_id to the second column and fps_var to the third column
	cols = list(merged_df.columns)
	cols.insert(0, cols.pop(cols.index('region_id')))
	cols.insert(1, cols.pop(cols.index('sample_id')))
	cols.insert(2, cols.pop(cols.index('fps_var')))
	merged_df = merged_df.loc[:, cols]

	if os.path.isfile(f'{output_path}/graphs/variance-scatterplots/{motif_id}_afps-jointplot-variance-alt.png'):
		print(f'WARNING: {motif_id}_afps-jointplot-variance-alt.png already exists! Skipping...')
	else:
		print(f'Saving {motif_id}_afps-jointplot-variance-alt.png to file...')
		# plot scatterplot of AF variance vs FPS variance
		sns.set_context("poster", rc={"figure.dpi": 300})
		g = sns.jointplot(data=merged_df, x='AF', y='fps_var', hue='sample_id', palette='bright', height=16, alpha=0.8)

		plt.xlabel('Allelic frequency (AF)')
		plt.ylabel('Footprint score (FPS) variance')
		g.fig.suptitle(f'AF vs FPS variance of {motif_id} sites', y=1)
		plt.savefig(f'{output_path}/graphs/variance-scatterplots/{motif_id}_afps-jointplot-variance-alt.png', dpi=300, bbox_inches='tight')
		plt.close()

	return merged_df, motif_id
 
def filter_top5percent(regsorted_df, output_path, motif_id):
    # truncate the dataframe to the first 5% of region_id values of all unique region_id values
	regsorted_df_filt = regsorted_df.head(int(len(regsorted_df['region_id'].unique())*0.05)*5)
	regsorted_df_filt.to_csv(f'{output_path}/output-data/afps_region-sorted/{motif_id}_afps-regionsorted-top5pc.tsv', sep='\t', index=False)
	# return the path to the saved file
	file_output = f'{output_path}/output-data/afps_region-sorted/{motif_id}_afps-regionsorted-top5pc.tsv'
	return file_output

def scale_fps(regsorted_df):
	# scale the FPS values to a range of 0-1
	# Initialize a MinMaxScaler
	scaler = MinMaxScaler()

	# Fit the MinMaxScaler to the 'FPS' column and transform it
	regsorted_df['FPS_scaled'] = scaler.fit_transform(regsorted_df[['FPS']])
 
	return regsorted_df

def filter_df(regsorted_df, iqr_filter=False):
	if iqr_filter == True:
		# filter out outliers using IQR method, but the inverse, where we retain the outliers and remove the inliers
		# Calculate Q1, Q3 and IQR for the 'AF' column
		Q1 = regsorted_df['FPS_scaled'].quantile(0.25)
		Q3 = regsorted_df['FPS_scaled'].quantile(0.75)
		# then filter the inliers
		iqr_filtered_df = regsorted_df[(regsorted_df['AF'] > 0.5) & ((regsorted_df['FPS_scaled'] <= Q1) | (regsorted_df['FPS_scaled'] >= Q3))]
		return iqr_filtered_df
	else:
		# filter for just AF > 0.5
		af_filtered_df = regsorted_df[regsorted_df['AF'] > 0.5]
		return af_filtered_df
    
def plot_jointplot(df, output_path, motif_id, subfolder='scaled-jointplots', out_suffix='scaled'):
	try:
  		# Create a jointplot of 'AF' and 'FPS'
		sns.jointplot(data=df, x='AF', y='FPS_scaled', kind='scatter', hue='sample_id', height=12)
		#save the plot
		plt.savefig(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-jointplot-{out_suffix}.png', dpi=300, bbox_inches='tight')
		plt.close()
		return True
	except Exception as e:
		print(f'ERROR: {e}')
		return False
 
def plot_lmplot(df, output_path, motif_id, type='split', subfolder='scaled-lmplots'):
	try:
		if type != 'split':
 			# Create a lmplot of 'AF' and 'FPS'
			sns.set_context("talk", rc={"figure.dpi": 300})
			sns.lmplot(data=df, x='AF', y='FPS_scaled', hue='sample_id', height=12)
			#save the plot
			plt.savefig(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-lmplot-scaled.png', dpi=300, bbox_inches='tight')
			plt.close()
			return True
		else:
			sns.set_context("talk", rc={"figure.dpi": 300})
			sns.lmplot(data=df, x='AF', y='FPS_scaled', col='sample_id', height=8, facet_kws=dict(sharex=False, sharey=False))
			#save the plot
			plt.savefig(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-lmplot-scaled-split.png', dpi=300, bbox_inches='tight')
			plt.close()
			return True
	except Exception as e:
		print(f'ERROR: {e}')
		return False
 
def process_data_into_variance(target_file, output_path):
	# load and process input file
	print(f'Processing {target_file}...')
	afps_variance_alt_df, motif_id = variance_calc_alt_df(target_file, output_path)
	# save the merged variance dataframe to file
	afps_variance_alt_df.to_csv(f'{output_path}/output-data/variance_dataframes/{motif_id}_afps-variance-alt.tsv', sep='\t', index=False)
	return True


def process_input_tsv(root_dir):
    # Find all *.tsv files in root_dir
    target_dir = Path(root_dir)
    tsv_files = target_dir.glob('*.tsv')
    return tsv_files

##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 2:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_analysis_pipeline-variance-alt.py <root_dir>')
	sys.exit(1)
else:
    root_dir = sys.argv[1]

##################
# define globals #
##################
output_path = '/home/msazizan/hyperspace/gatk-workflow/plotting'

if __name__ == '__main__':
	inputs = process_input_tsv(root_dir)
	# for target_file in inputs:
	# 	process_data(target_file, output_path)
    # run concurrent processes
	with cf.ProcessPoolExecutor(max_workers=5) as executor:
		executor.map(process_data_into_variance, inputs, it.repeat(output_path))
	print ("Pipeline finished! All footprint matrices have been processed with an alternative approach.")