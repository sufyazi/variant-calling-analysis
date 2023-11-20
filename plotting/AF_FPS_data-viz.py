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

def process_input_tsv(root_dir):
    # Find all *.tsv files in root_dir
    target_dir = Path(root_dir)
    tsv_files = target_dir.glob('*.tsv')
    return tsv_files

def process_data(target_file, output_path, plot=True):
	print(f'Processing {target_file}...')
	# import the data
	matrix_afps = pd.read_csv(target_file, sep='\t')
	# extract motif id from filename
	motif_id = os.path.basename(target_file).replace('_fpscore-af-varsites-combined-matrix-wide.tsv', '')
	print(f'{motif_id} data has been loaded.')
	afps_df = matrix_afps.filter(regex='_AF$|_fps$|_id$').copy()
	# convert to long format
	# convert to long format
	afps_df_long = afps_df.melt(id_vars=["region_id"], var_name="variable", value_name="value")

	# split the variable column into sample_id and type columns using reverse split string method, which returns a dataframe of columns based on the number of splits (n=x); this can directly be assigned to new columns in the original dataframe
	afps_df_long[['sample_id', 'type']] = afps_df_long['variable'].str.rsplit('_', n=1, expand=True)

	# drop the redundant 'variable' column
	afps_df_long = afps_df_long.drop(columns=["variable"])

	# now pivot the dataframe to create new columns based on the type column
	afps_df_lpv = afps_df_long.pivot(index=['region_id', 'sample_id'], columns='type', values='value').reset_index()

	# remove the index name and rename the columns to match the type values
	afps_df_lpv = afps_df_lpv.rename_axis(None, axis=1).rename(columns={'fps': 'FPS'})

	# sort the dataframe by region_id naturally
	afps_df_lpv = afps_df_lpv.reindex(index=index_natsorted(afps_df_lpv['region_id']))
	afps_df_lpv = afps_df_lpv.reset_index(drop=True)

	# use MinMaxScaler to scale the raw fps values to range between 0 and 1
	# scale the FPS values to a range of 0-1
	# Initialize a MinMaxScaler
	scaler = MinMaxScaler()
	# copy df
	fps_df_scaled = matrix_afps.filter(regex='_fps$|_id$').copy()
	# set the index to 'region_id'
	fps_df_scaled = fps_df_scaled.set_index('region_id')
	# Fit the MinMaxScaler to the 'FPS' column and transform it
	fps_df_scaled = pd.DataFrame(scaler.fit_transform(fps_df_scaled), columns=fps_df_scaled.columns, index=fps_df_scaled.index)
	# rename columns by adding '_scaled' to the column names
	fps_df_scaled = fps_df_scaled.add_suffix('_scaled')
	
	# reset index
	fps_df_scaled_long = fps_df_scaled.reset_index()
	# convert to long format
	fps_df_scaled_long = fps_df_scaled_long.melt(id_vars=["region_id"], var_name="variable", value_name="value")
	# split the variable column into sample_id and type columns using reverse split string method, which returns a dataframe of columns based on the number of splits (n=x); this can directly be assigned to new columns in the original dataframe
	# Split the 'variable' column into three parts
	fps_df_scaled_long[['part1', 'part2', 'part3']] = fps_df_scaled_long['variable'].str.rsplit('_', n=2, expand=True)
	# Assign part1 to 'sample_id' and concatenate the other parts to form 'type'
	fps_df_scaled_long['sample_id'] = fps_df_scaled_long['part1']
	fps_df_scaled_long['type'] = fps_df_scaled_long['part2'].str.upper() + '_' + fps_df_scaled_long['part3']
	# Drop the unnecessary columns
	fps_df_scaled_long = fps_df_scaled_long.drop(['variable', 'part1', 'part2', 'part3'], axis=1)
	# now pivot the dataframe to create new columns based on the type column
	fps_df_scaled_lpv = fps_df_scaled_long.pivot(index=['region_id', 'sample_id'], columns='type', values='value').reset_index()
	# remove the index name and rename the columns to match the type values
	fps_df_scaled_lpv = fps_df_scaled_lpv.rename_axis(None, axis=1)
	# sort the dataframe by region_id naturally
	fps_df_scaled_lpv = fps_df_scaled_lpv.reindex(index=index_natsorted(fps_df_scaled_lpv['region_id']))
	fps_df_scaled_lpv = fps_df_scaled_lpv.reset_index(drop=True)

	# merge the two dataframes
	afps_full_dfl = afps_df_lpv.merge(fps_df_scaled_lpv, on=['region_id', 'sample_id'])
	# save to file
	afps_full_dfl.to_csv(f'/home/msazizan/hyperspace/gatk-workflow/plotting/output-data/AF-FPS_regionsorted_full_longtable/{motif_id}_afps_full_scaled_longtable.tsv', sep='\t', index=False)
 
	# extract fps columns
	fps_df = matrix_afps.filter(regex='_fps$|_id$')
	# calculate variance of fps values across samples per region_id and add to a new column called 'fps_var'
	fps_df = fps_df.set_index('region_id')
	fps_df['FPS_var'] = fps_df.var(axis=1)
	# calculate the coefficient of variation (CV) of fps values across samples per region_id and add to a new column called 'fps_cv'
	fps_df['FPS_cv'] = fps_df.drop(columns=['FPS_var']).std(axis=1) / fps_df.drop(columns=['FPS_var']).mean(axis=1) * 100
	# calculate the quartile coefficient of dispersion (QCD) of fps values across samples per region_id and add to a new column called 'fps_qcd'
	fps_df['FPS_qcd'] = (fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.75, axis=1) - fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.25, axis=1)) / (fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.75, axis=1) + fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.25, axis=1))
 
	# calculate variance of fps_scaled values across samples per region_id and add to a new column called 'fps_scaled_var'
	fps_df_scaled['FPS_scaled_var'] = fps_df_scaled.var(axis=1)
	# calculate the coefficient of variation (CV) of fps values across samples per region_id and add to a new column called 'fps_scaled_cv'
	fps_df_scaled['FPS_scaled_cv'] = fps_df_scaled.drop(columns=['FPS_scaled_var']).std(axis=1) / fps_df_scaled.drop(columns=['FPS_scaled_var']).mean(axis=1) * 100
	# calculate the quartile coefficient of dispersion (QCD) of fps values across samples per region_id and add to a new column called 'fps_scaled_qcd'
	fps_df_scaled['FPS_scaled_qcd'] = (fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.75, axis=1) - fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.25, axis=1)) / (fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.75, axis=1) + fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.25, axis=1))
 
	# subset only the fps stats columns
	fps_stats_df = fps_df.filter(regex='_var$|_cv$|_qcd$|_id$').copy()
	fps_scaled_stats_df = fps_df_scaled.filter(regex='_var$|_cv$|_qcd$|_id$').copy()
	# merge on region_id index from both tables
	fps_stats_all_df = fps_stats_df.merge(fps_scaled_stats_df, left_index=True, right_index=True)

	afps_full_dfl_ind = afps_full_dfl.set_index('region_id')
	# merge afps_full_dfl_ind with stats df on region_id
	fps_merged_stats = afps_full_dfl_ind.merge(fps_stats_all_df, left_index=True, right_index=True, how='left')

	# sort naturally by region_id
	fps_merged_stats_sorted =fps_merged_stats.reset_index().reindex(index=index_natsorted(fps_merged_stats.index))
 
	# save file 
	fps_merged_stats_sorted.to_csv(f'/home/msazizan/hyperspace/gatk-workflow/plotting/output-data/AF-FPS_regionsorted_full_longtable/{motif_id}_afps_full_scaled_longtable_with_stats.tsv', sep='\t', index=False)
 
	# filter out unique region_id rows that have fps == 0
	# group by 'region_id' first 
	merged_filt = fps_merged_stats_sorted.groupby('region_id').filter(lambda x: x['FPS'].sum() > 0)
	
	# filter out unique region_id rows that have AF == 0 in all subtypes (sample_id grouping)
	# this has to be per group, so we need to groupby first
	merged_filt = merged_filt.groupby('region_id').filter(lambda x: (x['AF'] != 0).any())
	
	# for each unique region_id, find those that has AF == 0 in at least one subtype
	# this has to be per group, so we need to groupby first
	region_af_zero = merged_filt.groupby('region_id').filter(lambda x: (x['AF'] == 0).any())

	# then find the max AF value for each region_id and add to a new column called 'max_AF'
	df = region_af_zero.groupby('region_id').agg({'AF': 'max'}).rename(columns={'AF': 'max_AF'}).reset_index()
	max_af = df.reindex(index=index_natsorted(df['region_id']))
	# now return a boolean mask for regions that have max_AF <= 0.5
	mask = max_af['max_AF'] <= 0.5
	# this is done by subsetting the long dataframe with an expression that makes use of isin() method on the max_af dataframe masked by the boolean series, and then taking the inverse of the expression using ~
	merged_filt = merged_filt[~(merged_filt['region_id'].isin(max_af[mask]['region_id']))]
	
	# copy the filtered dataframe
	mf_df = merged_filt.copy()
	mf_df = mf_df.reset_index(drop=True)
	# calculate AF median per region_id
	mf_df['AF_median'] = mf_df.groupby('region_id')['AF'].transform('median')
 
	# calculate the variance of AF per region_id and add to a new column called 'AF_var'
	mf_df['AF_var'] = mf_df.groupby('region_id')['AF'].transform('var')
	mf_df['AF_cv'] = mf_df.groupby('region_id')['AF'].transform(lambda x: x.std() / x.mean() * 100)

	# set a small constant for QCD calculation to prevent NaN due to division by zero
	constant = 0.0001
	mf_df['AF_qcd'] = mf_df.groupby('region_id')['AF'].transform(lambda x: (x.quantile(q=0.75) - x.quantile(q=0.25)) / (x.quantile(q=0.75) + x.quantile(q=0.25) + constant))
	
	# to ensure that each unique region_id is retained as a group of subtype rows, we need to filter after grouping per region_id
	high_af = mf_df.groupby('region_id').filter(lambda x: (x['AF_median'] > 0.5).any())
	low_af = mf_df.groupby('region_id').filter(lambda x: (x['AF_median'] <= 0.5).any())

	if plot == True:
		# plot violin plot for sites with at least one subtype with AF == 0
		plt.figure(figsize=(10, 10), dpi=300)
		plt.subplot(2, 1, 1)
		sns.violinplot(x='region_id', y='AF', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), color='lightgray', inner='quartile', linecolor='black', linewidth=1.5)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='AF', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), hue='sample_id', size=4, jitter=True, palette='bright')

		# plot legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=1, fontsize=8)
		plt.xlabel('')
		plt.ylabel('Allele frequency (AF)', fontsize=8)

		# then plot the scaled FPS values for these sites
		plt.subplot(2, 1, 2)
		sns.violinplot(x='region_id', y='FPS_scaled', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), color='lightgray', inner='quartile', linecolor='black', linewidth=1.5)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='FPS_scaled', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), hue='sample_id', size=4, jitter=True, palette='bright', legend=False)

		# plot legend outside of the plot
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (at least one subtype with AF = 0)', fontsize=8)
		plt.ylabel(f'Footprint scores (FPS) (min-max scaled)', fontsize=8)
		plt.subplots_adjust(hspace=0.5)
		
  		# check existence of output directory
		if not os.path.exists(f'{output_dir}/output-data/plots/{motif_id}'):
			os.makedirs(f'{output_dir}/output-data/plots/{motif_id}')
		# save the plot
		plt.savefig(f'{output_dir}/output-data/plots/{motif_id}/{motif_id}_filtered_sites_with_at_least_subtype_with_AF0.pdf', dpi=300)
		# close the plot
		plt.close()

		# plot boxplot of AF distribution across high_af sites and the corresponding FPS_scaled variances
		plt.figure(figsize=(10, 10), dpi=300)
	# specify subplot
	plt.subplot(2, 1, 1)
	# plot count plot for high AF set
	sns.boxplot(x='region_id', y='AF', data=high_af, color='gold', showfliers=False)
	plt.xticks(rotation=90, fontsize=2)
	plt.xlabel(f'{motif_id} TF binding sites with allelic variants (median AF > 0.5)', fontsize=8)
	plt.ylabel('AF Distribution Per Site', fontsize=8)

	plt.subplot(2, 1, 2)
	# plot count plot for high AF set
	sns.barplot(x='region_id', y='FPS_scaled_var', data=high_af, color='darkmagenta')
	plt.xticks(rotation=90, fontsize=2)
	plt.xlabel(f'{motif_id} TF binding sites with allelic variants (median AF > 0.5)', fontsize=8)
	plt.ylabel('FPS (scaled) Variance', fontsize=8)
	
	# save the plot
	plt.savefig(f'{output_dir}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars.pdf', dpi=300)
	# close the plot
	plt.close()
 
##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 2:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_data-viz.py <root_dir>')
	sys.exit(1)
else:
    root_dir = sys.argv[1]
    output_dir = sys.argv[2]

##################
# define globals #
##################

if __name__ == '__main__':
	inputs = process_input_tsv(root_dir)
	# uncomment this to run serially
	# for target_file in inputs:
	# 	process_data(target_file, output_path)
 
	# uncomment this to run in parallel
	with cf.ProcessPoolExecutor(max_workers=8) as executor:
		executor.map(process_data, inputs, it.repeat(output_dir), it.repeat(True))
	print ("Pipeline finished! All footprint matrices have had the site variances calculated.")