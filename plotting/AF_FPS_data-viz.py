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

	os.makedirs(f'{output_path}/output-data/tables/{motif_id}', exist_ok=True)
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_full_scaled_longtable.tsv'):
		# save to file
		afps_full_dfl.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_full_scaled_longtable.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table already exists. Skipping...')
	# extract fps columns
	fps_df = matrix_afps.filter(regex='_fps$|_id$')
	# calculate variance of fps values across samples per region_id and add to a new column called 'fps_var'
	fps_df = fps_df.set_index('region_id')
	fps_df['FPS_var'] = fps_df.var(axis=1)
	# calculate the coefficient of variation (CV) of fps values across samples per region_id and add to a new column called 'fps_cv'
	# fps_df['FPS_cv'] = fps_df.drop(columns=['FPS_var']).std(axis=1) / fps_df.drop(columns=['FPS_var']).mean(axis=1) * 100
	# calculate the quartile coefficient of dispersion (QCD) of fps values across samples per region_id and add to a new column called 'fps_qcd'
	# fps_df['FPS_qcd'] = (fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.75, axis=1) - fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.25, axis=1)) / (fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.75, axis=1) + fps_df.drop(columns=['FPS_var', 'FPS_cv']).quantile(q=0.25, axis=1))

	# calculate variance of fps_scaled values across samples per region_id and add to a new column called 'fps_scaled_var'
	fps_df_scaled['FPS_scaled_var'] = fps_df_scaled.var(axis=1)
	# calculate the coefficient of variation (CV) of fps values across samples per region_id and add to a new column called 'fps_scaled_cv'
	# fps_df_scaled['FPS_scaled_cv'] = fps_df_scaled.drop(columns=['FPS_scaled_var']).std(axis=1) / fps_df_scaled.drop(columns=['FPS_scaled_var']).mean(axis=1) * 100
	# calculate the quartile coefficient of dispersion (QCD) of fps values across samples per region_id and add to a new column called 'fps_scaled_qcd'
	# fps_df_scaled['FPS_scaled_qcd'] = (fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.75, axis=1) - fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.25, axis=1)) / (fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.75, axis=1) + fps_df_scaled.drop(columns=['FPS_scaled_var', 'FPS_scaled_cv']).quantile(q=0.25, axis=1))

	# subset only the fps stats columns
	# fps_stats_df = fps_df.filter(regex='_var$|_cv$|_qcd$|_id$').copy()
	# fps_scaled_stats_df = fps_df_scaled.filter(regex='_var$|_cv$|_qcd$|_id$').copy()
	fps_stats_df = fps_df.filter(regex='_var$|_id$').copy()
	fps_scaled_stats_df = fps_df_scaled.filter(regex='_var$|_id$').copy()
	# merge on region_id index from both tables
	fps_stats_all_df = fps_stats_df.merge(fps_scaled_stats_df, left_index=True, right_index=True)

	afps_full_dfl_ind = afps_full_dfl.set_index('region_id')
	# merge afps_full_dfl_ind with stats df on region_id
	fps_merged_stats = afps_full_dfl_ind.merge(fps_stats_all_df, left_index=True, right_index=True, how='left')

	# sort naturally by region_id
	fps_merged_stats_sorted =fps_merged_stats.reset_index().reindex(index=index_natsorted(fps_merged_stats.index))

	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_full_scaled_longtable_with_stats-var.tsv'):
		# save file
		fps_merged_stats_sorted.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_full_scaled_longtable_with_stats-var.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table with stats already exists. Skipping...')
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
	# mf_df['AF_cv'] = mf_df.groupby('region_id')['AF'].transform(lambda x: x.std() / x.mean() * 100)

	# set a small constant for QCD calculation to prevent NaN due to division by zero
	# constant = 0.0001
	# mf_df['AF_qcd'] = mf_df.groupby('region_id')['AF'].transform(lambda x: (x.quantile(q=0.75) - x.quantile(q=0.25)) / (x.quantile(q=0.75) + x.quantile(q=0.25) + constant))

	# to ensure that each unique region_id is retained as a group of subtype rows, we need to filter after grouping per region_id
	high_af = mf_df.groupby('region_id').filter(lambda x: (x['AF_median'] > 0.5).any())
	# low_af = mf_df.groupby('region_id').filter(lambda x: (x['AF_median'] <= 0.5).any())

	# subset high_af for sites with FPS_var > 75th percentile
	high_af_fps_var = high_af[high_af['FPS_scaled_var'] > high_af['FPS_scaled_var'].quantile(q=0.75)]

	# sort the region_id in the filtered high_af dataframe based on descending order of AF_var
	high_af_uniq_reg = high_af_fps_var[['region_id', 'AF_var', 'FPS_scaled_var']].drop_duplicates()
	afvar_uniqsort = high_af_uniq_reg.sort_values(by='AF_var', ascending=False)
	fpsvar_uniqsort = high_af_uniq_reg.sort_values(by='FPS_scaled_var', ascending=False)

	# extract the region_id from the high_af_fps_var_df_uniqsort dataframe
	fps_sorted_region_ids = fpsvar_uniqsort['region_id']
	# Change 'region_id' to a categorical variable with the categories ordered by 'fps_sorted_region_ids'
	high_af['region_id'] = pd.Categorical(high_af['region_id'], categories=fps_sorted_region_ids, ordered=True)
	# Filter the DataFrame
	high_af_fpsvar_filt = high_af[high_af['region_id'].isin(fps_sorted_region_ids)]
	# Sort the DataFrame by 'region_id'
	high_af_fpsvar_filt = high_af_fpsvar_filt.sort_values('region_id')

	# Find the index of the max FPS_scaled value for each region_id
	idx_max = high_af_fpsvar_filt.groupby('region_id', observed=True)['FPS_scaled'].idxmax()
	idx_min = high_af_fpsvar_filt.groupby('region_id', observed=True)['FPS_scaled'].idxmin()
	# Select the corresponding rows
	max_fps_scaled = high_af_fpsvar_filt.loc[idx_max]
	min_fps_scaled = high_af_fpsvar_filt.loc[idx_min]

	# do the same for AF
	idx_max_af = high_af_fpsvar_filt.groupby('region_id', observed=True)['AF'].idxmax()
	idx_min_af = high_af_fpsvar_filt.groupby('region_id', observed=True)['AF'].idxmin()
	# Select the corresponding rows
	max_af_raw = high_af_fpsvar_filt.loc[idx_max_af]
	min_af_raw = high_af_fpsvar_filt.loc[idx_min_af]

	# create masks for the inverse of max and min values
	mask = ~high_af_fpsvar_filt.index.isin(idx_max)
	max_fps_scaled_inv = high_af_fpsvar_filt[mask]

	mask = ~high_af_fpsvar_filt.index.isin(idx_min)
	min_fps_scaled_inv = high_af_fpsvar_filt[mask]

	# do the same for AF
	mask = ~high_af_fpsvar_filt.index.isin(idx_max_af)
	max_af_raw_inv = high_af_fpsvar_filt[mask]

	mask = ~high_af_fpsvar_filt.index.isin(idx_min_af)
	min_af_raw_inv = high_af_fpsvar_filt[mask]

	if plot == True:
		# define color palettes
		dutchfield = ["#e60049", "#0bb4ff", "#87bc45", "#ef9b20", "#b33dc6"]
		springpastel = ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a"]
		# create color dictionary
		color_dict = {'S6R691V_her2': "#e60049", 'ANAB5F7_basal': "#0bb4ff", '98JKPD8_lumA': "#87bc45", 'PU24GB8_lumB': "#ef9b20", '2GAMBDQ_norm': "#b33dc6"}
		gray = 'lightgray'
		color_gray_dict = {'S6R691V_her2': gray, 'ANAB5F7_basal': gray, '98JKPD8_lumA': gray, 'PU24GB8_lumB': gray, '2GAMBDQ_norm': gray}
  		####################################### plot violin plot for sites with at least one subtype with AF == 0
		plt.figure(figsize=(10, 10), dpi=300)
		plt.subplot(2, 1, 1)
		sns.violinplot(x='region_id', y='AF', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), color='gainsboro', inner='quartile', linecolor='black', linewidth=1.5)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='AF', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), hue='sample_id', size=4, jitter=True, palette=dutchfield)

		# plot legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=1, fontsize=8)
		plt.xlabel('')
		plt.ylabel('Allele frequency (AF)', fontsize=8)

		# then plot the scaled FPS values for these sites
		plt.subplot(2, 1, 2)
		sns.violinplot(x='region_id', y='FPS_scaled', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), color='gainsboro', inner='quartile', linecolor='black', linewidth=1.5)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='FPS_scaled', data=mf_df.groupby('region_id').filter(lambda x: (x['AF'] == 0).any()), hue='sample_id', size=4, jitter=True, palette=dutchfield, legend=False)

		# plot legend outside of the plot
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (at least one subtype with AF = 0)', fontsize=8)
		plt.ylabel(f'Footprint scores (FPS) (min-max scaled)', fontsize=8)
		plt.subplots_adjust(hspace=0.5)

		# check existence of output directory
		os.makedirs(f'{output_path}/output-data/plots/{motif_id}', exist_ok=True)
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_filtered_sites_with_at_least_subtype_with_AF0.pdf'):
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_filtered_sites_with_at_least_subtype_with_AF0.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
		else:
			print(f'{motif_id} subset violin plot already exists. Skipping...')
			plt.close()

		###################################### plot boxplot of AF distribution on filtered sorted sites as well as the substype-hued stripplot on top
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=high_af_fpsvar_filt, hue='sample_id', palette=dutchfield, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		plt.ylabel('AF Distribution Per Site', fontsize=8)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)

		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_af_fpsvar_filt, color='palevioletred', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=8)

		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=high_af_fpsvar_filt, hue='sample_id', palette=dutchfield, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		plt.ylabel('FPS Scaled Distribution Per Site', fontsize=8)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_af_fpsvar_filt, color='mediumpurple', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (median AF > 0.5)', fontsize=8)
		plt.ylabel('FPS Variance', fontsize=8)
		plt.subplots_adjust(hspace=0.5)

		# check existence of output directory
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots.pdf'):
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
		else:
			print(f'{motif_id} subset sorted violin plot already exists. Skipping...')
			plt.close()

		###################################### plot boxplot of AF distribution on filtered sorted sites as well as the substype-hued stripplot on top (with maxima)
		import textwrap
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=max_af_raw_inv, hue='sample_id', palette=color_gray_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=max_af_raw, hue='sample_id', palette=color_dict, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		ylabel = textwrap.fill('AF Per Site (maxima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)
		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_af_fpsvar_filt, color='palevioletred', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=8)

		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps_scaled_inv, hue='sample_id', palette=color_gray_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps_scaled, hue='sample_id', palette=color_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS Per Site (maxima)', width=20)
		plt.ylabel(ylabel, fontsize=8)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_af_fpsvar_filt, color='mediumpurple', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (median AF > 0.5)', fontsize=8)
		plt.ylabel('FPS Variance', fontsize=8)
		plt.subplots_adjust(hspace=0.5)

		# save the plot
		# check existence of output directory
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots-maxima.pdf'):
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots-maxima.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
		else:
			print(f'{motif_id} subset sorted integrated plots of maxima already exist. Skipping...')
			plt.close()

		###################################### plot boxplot of AF distribution on filtered sorted sites as well as the substype-hued stripplot on top (with minima)
		import textwrap
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=min_af_raw_inv, hue='sample_id', palette=color_gray_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=min_af_raw, hue='sample_id', palette=color_dict, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		ylabel = textwrap.fill('AF Per Site (minima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)
		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_af_fpsvar_filt, color='palevioletred', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=8)

		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_af_fpsvar_filt, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps_scaled_inv, hue='sample_id', palette=color_gray_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps_scaled, hue='sample_id', palette=color_dict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS Per Site (minima)', width=20)
		plt.ylabel(ylabel, fontsize=8)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_af_fpsvar_filt, color='mediumpurple', edgecolor='black')
		plt.xticks(rotation=90, fontsize=2)
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (median AF > 0.5)', fontsize=8)
		plt.ylabel('FPS Variance', fontsize=8)
		plt.subplots_adjust(hspace=0.5)

		# save the plot
		# check existence of output directory
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots-minima.pdf'):
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_dist_boxplot_persite_with_FPS_scaled_vars-filtsorted-dots-minima.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
		else:
			print(f'{motif_id} subset sorted integrated plots of minima already exist. Skipping...')
			plt.close()

	# quantify the number of filtered sites per sample_id
	max_af_raw_subset = max_af_raw[['region_id', 'sample_id']]
	max_af_raw_df = max_af_raw_subset.groupby('sample_id')['region_id'].count().to_frame()
	max_fps_scaled_subset = max_fps_scaled[['region_id', 'sample_id']]
	max_fps_scaled_df = max_fps_scaled_subset.groupby('sample_id')['region_id'].count().to_frame()

	# save high_af_fpsvar_filt and max_af_raw and max_fps_scaled to csv
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_af_filtsorted_longtable.tsv'):
		high_af_fpsvar_filt.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_af_filtsorted_longtable.tsv', sep='\t', index=True)
		max_af_raw_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_af_regions-unique.tsv', sep='\t', index=True)
		max_fps_scaled_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_fps_scaled_regions-unique.tsv', sep='\t', index=True)
	else:
		print(f'{motif_id} high_af_filtsorted_longtable already exists. Skipping...')
		print(f'{motif_id} maximum_af_regions-unique already exists. Skipping...')
		print(f'{motif_id} maximum_fps_scaled_regions-unique already exists. Skipping...')

##################
# load arguments #
##################
# check for the required arguments
if len(sys.argv) < 3:
	print(f'ERROR: Missing required arguments!')
	print(f'USAGE: python3 AF_FPS_data-viz.py <root_dir> <output_dir>')
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
	for target_file in inputs:
		process_data(target_file, output_dir, True)

	# uncomment this to run in parallel
	#with cf.ProcessPoolExecutor(max_workers=8) as executor:
	#	executor.map(process_data, inputs, it.repeat(output_dir), it.repeat(True))

	print ("Pipeline finished! All footprint matrices have been processed.")
