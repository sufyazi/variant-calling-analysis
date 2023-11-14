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


def regionsort_df(filepath, output_path):
	# extract motif id from filename
	motif_id = os.path.basename(filepath).replace('_fpscore-af-varsites-combined-matrix-wide.tsv', '')
	# print message
	print(f'Processing {motif_id}...')

	# load data file
	afps_df = pd.read_csv(filepath, sep='\t')
	# filter df
	afps_df = afps_df.filter(regex='_AF$|_fps$|_id$')
	# convert df to long format
	afps_dfl = afps_df.melt(id_vars=["region_id"], var_name="variable", value_name="value")

	# split the variable column into sample_id and type columns using reverse split string method, which returns a dataframe of columns based on the number of splits (n=x); this can directly be assigned to new columns in the original dataframe
	afps_dfl[['sample_id', 'type']] = afps_dfl['variable'].str.rsplit('_', n=1, expand=True)
	# drop the redundant 'variable' column
	afps_dfl = afps_dfl.drop(columns=["variable"])
	# now pivot the dataframe to create new columns based on the type column
	afps_df_lpv = afps_dfl.pivot(index=['region_id', 'sample_id'], columns='type', values='value').reset_index()
	# remove the index name and rename the columns to match the type values
	afps_df_lpv = afps_df_lpv.rename_axis(None, axis=1).rename(columns={'fps': 'FPS'})
	# sort the dataframe by region_id naturally
	afps_df_lpv = afps_df_lpv.reindex(index=index_natsorted(afps_df_lpv['region_id']))
	afps_df_lpv = afps_df_lpv.reset_index(drop=True)

	# generate a jointplot of the unfiltered data
	# check if the output file exists, if yes, skip plotting
	if os.path.isfile(f'{output_path}/graphs/unfilt-scatterplots/{motif_id}_afps-jointplot-unfilt.png'):
		print(f'Jointplot for {motif_id} already exists. Skipping...')
	else:
		print(f'Plotting jointplot for {motif_id}...')
		sns.jointplot(x="AF", y="FPS", data=afps_df_lpv, hue='sample_id', height=12)
		plt.savefig(f'{output_path}/graphs/unfilt-scatterplots/{motif_id}_afps-jointplot-unfilt.png', dpi=300, bbox_inches='tight')
		plt.close()
 
	# Calculate the cumulative 'AF' for each 'region_id'
	cumulative_af = afps_df_lpv.groupby('region_id', observed=True)['AF'].sum().reset_index().rename(columns={'AF': 'cumulative_AF'})
	# set the index to 'region_id' and then sort the dataframe by 'AF' in descending order
	cumulative_af = cumulative_af.set_index('region_id').sort_values(by='cumulative_AF', ascending=False)
	# now use the index of cumulative_af to sort afps_df_lpv by 'region_id' in descending order, and then sort the sample_id 
	# Create a categorical variable with ordered categories
	afps_df_lpv['region_id'] = pd.Categorical(afps_df_lpv['region_id'], categories=cumulative_af.index.unique(), ordered=True)
	# Sort by the categorical 'region_id'
	afps_df_lpv = afps_df_lpv.sort_values('region_id')
	# get unique sample_id values into a list to define a categorical order
	datasets = afps_df_lpv['sample_id'].unique().tolist()
	datasets = sorted(datasets)
	# Create a categorical variable with ordered categories
	afps_df_lpv['sample_id'] = pd.Categorical(afps_df_lpv['sample_id'], categories=datasets, ordered=True)
	# Sort 'sample_id' within each 'region_id'
	afps_df_regsorted = afps_df_lpv.groupby('region_id', sort=False, observed=False).apply(lambda x: x.sort_values('sample_id')).reset_index(drop=True)
 
	# save the sorted dataframe to file
	# check existence first
	if os.path.isfile(f'{output_path}/output-data/afps_region-sorted/{motif_id}_afps-regionsorted.tsv'):
		print(f'Sorted dataframe for {motif_id} already exists. Skipping...')
	else:
		print(f'Saving sorted dataframe for {motif_id}...')
		afps_df_regsorted.to_csv(f'{output_path}/output-data/afps_region-sorted/{motif_id}_afps-regionsorted.tsv', sep='\t', index=False)
	
	return afps_df_regsorted, motif_id
 
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
		# check if the output file exists, if yes, skip plotting
		if os.path.isfile(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-jointplot-{out_suffix}.png'):
			print(f'Jointplot for {motif_id} already exists. Skipping...')
			return True
		else:
			print(f'Plotting jointplot for {motif_id}...')
			sns.jointplot(data=df, x='AF', y='FPS_scaled', kind='scatter', hue='sample_id', height=12)
			#save the plot
			plt.ylim(0, 1)
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
			# check if the output file exists, if yes, skip plotting
			if os.path.isfile(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-lmplot-scaled.png'):
				print(f'lmplot for {motif_id} already exists. Skipping...')
				return True
			else:
				sns.set_context("talk", rc={"figure.dpi": 300})
				sns.lmplot(data=df, x='AF', y='FPS_scaled', hue='sample_id', height=12)
				#save the plot
				plt.savefig(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-lmplot-scaled.png', dpi=300, bbox_inches='tight')
				plt.close()
				return True
		else:
			# check if the output file exists, if yes, skip plotting
			if os.path.isfile(f'{output_path}/graphs/{subfolder}/{motif_id}_afps-lmplot-scaled-split.png'):
				print(f'lmplot for {motif_id} already exists. Skipping...')
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
 
def process_data(target_file, output_path, filter_fps=False):
	# load and process input file
	afps_regsorted_df, motif_id = regionsort_df(target_file, output_path)

	# plot the sorted stacked bar plot (unfiltered)
	# print(f'Plotting stacked barplot for all variant sites sorted by descending AF values...')
	# plot_stacked_barplot(afps_regsorted_df, motif_id, output_path)
	# print(f'Done!') 
 
	print('Filtering for top 5% variant sites...')
	# filter the top 5% of variant sites
	afps_top = filter_top5percent(afps_regsorted_df, output_path, motif_id)
	print(f'Now plotting stacked barplot for top 5% variant sites...')
	# load the filtered dataframe
	afps_df_top = pd.read_csv(afps_top, sep='\t')
	if plot_stacked_barplot(afps_df_top, motif_id, output_path, rotate_xticks=True, xticks_fontsize=5, output_subfolder='topsite-barplots', plot_suffix='top_5pc'):
		print(f'Done!')
	else:
		print(f'ERROR: Unable to plot stacked barplot for top 5% variant sites of {motif_id}')
	# scale and filter the FPS values
	print(f'Min-max scaling FPS values...')
	afps_regsorted_df = scale_fps(afps_regsorted_df)

	if filter_fps == True:
		iqr_filt_df = filter_df(afps_regsorted_df, iqr_filter=True)
		# plot the scaled data jointplot and the filtered data jointplot
		print(f'Plotting jointplots for scaled and filtered data...')
		if plot_jointplot(afps_regsorted_df, output_path, motif_id) and plot_jointplot(iqr_filt_df, output_path, motif_id, subfolder='scaled-filtered-jointplots', out_suffix='scaled-iqr-filtered'):
			print(f'Done!')
		else:
			print(f'ERROR: Unable to plot jointplots for either scaled or filtered data of {motif_id}')
		# plot the combined lmplot and separated lmplot
		print(f'Plotting lmplots for filtered data...')
		if plot_lmplot(iqr_filt_df, output_path, motif_id) and plot_lmplot(iqr_filt_df, output_path, motif_id, type='combined'):
			print(f'Done!')
		else:
			print(f'ERROR: Unable to plot lmplots for filtered data of {motif_id}')
   
		# finally extract the counts of filtered variant sites and save to file
		print(f'Extracting counts of filtered variant sites...')
		filtered_sites_df = iqr_filt_df.groupby('sample_id', observed=True)['region_id'].nunique().reset_index()
		filtered_sites_df['motif_id'] = motif_id
		try:
			filtered_sites_df.to_csv(f'{output_path}/output-data/iqr-filtered_varsites/{motif_id}_iqr-filt-varsite-counts.tsv', sep='\t', index=False)
			print(f'Done! IQR-filtered variant sites have been saved to file.')
			return True
		except Exception as e:
			print(f'ERROR: Unable to save IQR-filtered variant sites to file.')
			return False
	else:
		# filter df for AF > 0.5
		af_filtered_df = filter_df(afps_regsorted_df)
		# plot the scaled data jointplot
		print(f'Plotting jointplot for scaled FPS data...')
		if plot_jointplot(af_filtered_df, output_path, motif_id):
			print(f'Done!')
		else:	
			print(f'ERROR: Unable to plot jointplot for scaled FPS data of {motif_id}')
		# plot the combined lmplot and separated lmplot
		print(f'Plotting lmplots for scaled FPS data...')
		if plot_lmplot(af_filtered_df, output_path, motif_id) and plot_lmplot(af_filtered_df, output_path, motif_id, type='combined'):
			print(f'Done!')
		else:
			print(f'ERROR: Unable to plot lmplots for scaled FPS data of {motif_id}')

		# finally extract the counts of filtered variant sites and save to file
		print(f'Extracting counts of all variant sites when AF > 0.5...')
		filtered_sites_df = af_filtered_df.groupby('sample_id', observed=True)['region_id'].nunique().reset_index()
		filtered_sites_df['motif_id'] = motif_id
		try:
			filtered_sites_df.to_csv(f'{output_path}/output-data/af-filtered_varsites/{motif_id}_af-filt-varsite-counts.tsv', sep='\t', index=False)
			print(f'Done! AF-filtered variant sites have been saved to file.')
			return True
		except Exception as e:
			print(f'ERROR: Unable to save AF-filtered variant sites to file.')
			return False

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
	print(f'USAGE: python3 AF_FPS_analysis_pipeline.py <root_dir>')
	sys.exit(1)
else:
    root_dir = sys.argv[1]

##################
# define globals #
##################
output_path = '/home/msazizan/hyperspace/gatk-workflow/plotting/test-ground'

if __name__ == '__main__':
	inputs = process_input_tsv(root_dir)
	filter_fps = True
	# for target_file in inputs:
	# 	process_data(target_file, output_path)
    # run concurrent processes
	with cf.ProcessPoolExecutor(max_workers=5) as executor:
		executor.map(process_data, inputs, it.repeat(output_path), it.repeat(filter_fps))
	print ("Pipeline finished! All footprint matrices have been processed.")