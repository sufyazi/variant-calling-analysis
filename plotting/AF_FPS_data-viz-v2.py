#!/usr/bin/env python3

####################
# import libraries #
####################

import os
import sys
import textwrap
import matplotlib
matplotlib.use("Agg")
import pandas as pd
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

def load_data(tsv_filepath):
	# import the data
	matrix_afps = pd.read_csv(tsv_filepath, sep='\t')
	# extract motif id from filename
	motif_id = os.path.basename(tsv_filepath).replace('_fpscore-af-varsites-combined-matrix-wide.tsv', '')
	print(f'{motif_id} data has been loaded.')
	afps_df = matrix_afps.filter(regex='_AF$|_fps$|_id$').copy()
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
	print(f'{motif_id} matrix has been loaded and converted to long format.')
	return matrix_afps, motif_id, afps_df_lpv

def scale_data(matrix):
	# scale the FPS values to a range of 0-1
	# Initialize a MinMaxScaler
	scaler = MinMaxScaler()
	# copy df
	fps_df_scaled = matrix.filter(regex='_fps$|_id$').copy()
	# set the index to 'region_id'
	fps_df_scaled = fps_df_scaled.set_index('region_id')
	# Fit the MinMaxScaler to the 'FPS' column and transform it
	fps_df_scaled = pd.DataFrame(scaler.fit_transform(fps_df_scaled), columns=fps_df_scaled.columns, index=fps_df_scaled.index)
	# rename columns by adding '_scaled' to the column names
	fps_df_scaled = fps_df_scaled.add_suffix('_scaled')
	##### Now convert to long format #####
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
	return fps_df_scaled, fps_df_scaled_lpv

def variance_stats(matrix, motif_id, fps_df_scaled):
	# calculate the variance of AF and FPS values across samples per region_id
	# extract fps columns
	fps_df = matrix.filter(regex='_fps$|_id$').copy()
	print(f'Calculating {motif_id} FPS and AF variance statistics...')
	# calculate variance of fps values across samples per region_id and add to a new column called 'fps_var'
	fps_df = fps_df.set_index('region_id')
	fps_df['FPS_var'] = fps_df.var(axis=1)
	# now calculate variance of fps_scaled values across samples per region_id and add to a new column called 'fps_scaled_var'
	fps_df_scaled_cp = fps_df_scaled.copy()
	fps_df_scaled_cp['FPS_scaled_var'] = fps_df_scaled_cp.var(axis=1)
	# Do the same for AF values
	# extract af columns
	af_df = matrix.filter(regex='_AF$|_id$')
	af_df = af_df.set_index('region_id')
	# then calculate variance of af values across samples per region_id and add to a new column called 'af_var'
	af_df['AF_var'] = af_df.var(axis=1)
	# now merge the stats columns
	fps_stats_df = fps_df.filter(regex='_var$|_id$').copy()
	fps_scaled_stats_df = fps_df_scaled_cp.filter(regex='_var$|_id$').copy()
	af_stats_df = af_df.filter(regex='_var$|_id$').copy()
	# merge on region_id index from both fps tables
	fps_stats_df_temp = fps_stats_df.merge(fps_scaled_stats_df, left_index=True, right_index=True)
	# now merge with af_df
	afps_stats_df = fps_stats_df_temp.merge(af_stats_df, left_index=True, right_index=True)
	return afps_stats_df

def accessory_df(afps_stats_mergesorted):
	# subset a copy of the dataframe with only the AF_var and FPS_scaled_var columns and region_id
	unfiltered_df = afps_stats_mergesorted[['region_id', 'FPS_scaled_var', 'AF_var']].copy().drop_duplicates().reset_index(drop=True)
	# set region_id as index
	unfiltered_df = unfiltered_df.set_index('region_id')
	# filter rows using the condition FPS_scaled_var > 0.001 and AF_var > 0.001
	filtered_df = unfiltered_df[(unfiltered_df['FPS_scaled_var'] > 0.001) & (unfiltered_df['AF_var'] > 0.001)]
	return filtered_df

def accessory_plot(filtered_df, motif_id, output_path):
	################ PLOT ASIDE ################
	# plot scatter plot of AF_var vs FPS_scaled_var unfiltered
	# print(f'Plotting {motif_id} scatter plot of AF_var vs FPS_scaled_var (unfiltered)...')
	# g = sns.jointplot(data=unfiltered_df, x="FPS_scaled_var", y="AF_var", height=10, ratio=5, color='darkslateblue' )
	# g.fig.suptitle(f'AF vs FPS var for {motif_id} (unfiltered: {len(unfiltered_df)} regions)')
	# # save the plot
	# print(f'Saving {motif_id} scatter plot of AF and scaled FPS variance (unfiltered)...')
	# g.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-unfilt.pdf', dpi=300, bbox_inches="tight")
	# g.fig.clear()
	# now plot the filtered scatterplot
	print(f'Plotting {motif_id} scatter plot of AF_var vs FPS_scaled_var (filtered)...')
	g = sns.jointplot(data=filtered_df, x="FPS_scaled_var", y="AF_var", height=10, ratio=5, color='darkslateblue' )
	g.fig.suptitle(f'AF vs FPS var for {motif_id} (filtered: {len(filtered_df)} regions)')
	# save the plot
	print(f'Saving {motif_id} scatter plot of AF and scaled FPS variance (filtered)...')
	g.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-filt.pdf', dpi=300, bbox_inches="tight")
	# close the plot
	plt.close("all")
	del g
	print('Plot space closed and plots have been saved to file.')
	################ END PLOT ASIDE ################

def basic_filtering(afps_stats_mergesorted):
	# filter out unique region_id rows that have ALL fps == 0; group by 'region_id' first 
	merged_filt = afps_stats_mergesorted.groupby('region_id').filter(lambda x: x['FPS'].sum() > 0)
	# filter out unique region_id rows that have AF == 0 in all subtypes (sample_id)
	merged_filt_nozero = merged_filt.groupby('region_id').filter(lambda x: (x['AF'] != 0).any())
	# now for each unique region_id, find those that have AF == 0 in at least one subtype
	region_af_zero = merged_filt_nozero.groupby('region_id').filter(lambda x: (x['AF'] == 0).any())
	# then find the max AF value for each region_id and add to a new column called 'max_AF'
	df = region_af_zero.groupby('region_id').agg({'AF': 'max'}).rename(columns={'AF': 'max_AF'}).reset_index()
	max_af = df.reindex(index=index_natsorted(df['region_id']))
	# now return a boolean mask for regions that have max_AF <= 0.5
	mask = max_af['max_AF'] <= 0.5
	# filter the merged_filt_nozero such that rows that have `region_id` in the Boolean-masked `max_af` <= 0.5 dataframe will be discarded (so we retain only `region_id` rows that have `max_af` > 0.5)
	merged_filt_nolowaf = merged_filt_nozero[~(merged_filt_nozero['region_id'].isin(max_af[mask]['region_id']))]
	# additionally, filter out unique region_id rows that have AF_var == 0 in even one subtype (sample_id) regardless of max_AF
	merged_filt_nozeroaf = merged_filt_nozero[~(merged_filt_nozero['region_id'].isin(max_af['region_id']))]
	# also do the inverse filtering
	atleast_one_zero_af = merged_filt_nozero[merged_filt_nozero['region_id'].isin(max_af['region_id'])]
	return merged_filt_nolowaf, merged_filt_nozeroaf, atleast_one_zero_af

def thresholding_strat(high_af_df, threshold, low_af_df=None, unsplit_df=None):
	if threshold == 'iqr':
		# thresholding strategy 1: use IQR method to capture outliers and return only outlier regions
		# subset high_af for sites with FPS_var > 75th percentile + 1.5 * IQR
		iqr = high_af_df['FPS_scaled_var'].quantile(q=0.75) - high_af_df['FPS_scaled_var'].quantile(q=0.25)
		high_af_outliers = high_af_df[high_af_df['FPS_scaled_var'] > high_af_df['FPS_scaled_var'].quantile(q=0.75) + (1.5 * iqr)]
		print('Thresholding strategy 1: IQR method: Returning one dataframe only.')
		return high_af_outliers, None, None, None
	elif threshold == 'central':
		# thresholding strategy 2: compute mean of FPS_scaled data points from the unsplit (high+low af) dataframe and return only regions with FPS_scaled > mean
		global_mean = unsplit_df['FPS_scaled'].mean()
		high_af_abovemean = high_af_df.groupby('region_id').filter(lambda x: (x['FPS_scaled'] > global_mean).all())
		high_af_belowmean = high_af_df.groupby('region_id').filter(lambda x: (x['FPS_scaled'] <= global_mean).all())
		low_af_abovemean = low_af_df.groupby('region_id').filter(lambda x: (x['FPS_scaled'] > global_mean).all())
		low_af_belowmean = low_af_df.groupby('region_id').filter(lambda x: (x['FPS_scaled'] <= global_mean).all())
		print('Thresholding strategy 2: Central method: Returning two dataframes.')
		return high_af_abovemean, high_af_belowmean, low_af_abovemean, low_af_belowmean
	else:
		raise ValueError('Invalid thresholding strategy. Please choose either "iqr" or "central".')
	
def filtersort_df(input_df):
	# discard sites if it has AF_var == 0
	nzaf_df = input_df[input_df['AF_var'] != 0]
	# extract unique region_ids from the high_nzaf_outliers df
	af_filt_uniq_reg = nzaf_df[['region_id', 'AF_var', 'FPS_scaled_var']].drop_duplicates()
	# sort the region_id in the filtered high_af dataframe based on descending order of AF_var or FPS_scaled_var
	# afvar_uniqsort = high_af_uniq_reg.sort_values(by='AF_var', ascending=False)
	fpsvar_uniqsort = af_filt_uniq_reg.sort_values(by='FPS_scaled_var', ascending=False)
	# extract the region_id from the sorted dataframe
	sorted_region_ids = fpsvar_uniqsort['region_id']
	# Change 'region_id' to a categorical variable with the categories ordered by 'fps_sorted_region_ids'
	df_copy = nzaf_df.copy()
	df_copy['region_id'] = pd.Categorical(nzaf_df['region_id'], categories=sorted_region_ids, ordered=True)
	# Filter the DataFrame
	nzaf_df_filt = df_copy[df_copy['region_id'].isin(sorted_region_ids)]
	# Sort the DataFrame by 'region_id'
	nzaf_df_filt = nzaf_df_filt.sort_values('region_id')
	# get unique sample_id values into a list to define a categorical order
	datasets = nzaf_df_filt['sample_id'].unique().tolist()
	datasets = sorted(datasets)
	# Create a categorical variable with ordered categories
	dataset_copy = nzaf_df_filt.copy()
	dataset_copy['sample_id'] = pd.Categorical(dataset_copy['sample_id'], categories=datasets, ordered=True)
	# Sort 'sample_id' within each 'region_id'
	nzaf_df_filtsorted = dataset_copy.groupby('region_id', sort=False, observed=False).apply(lambda x: x.sort_values('sample_id')).reset_index(drop=True)
	return nzaf_df_filtsorted

def find_min_max(input_df):
	######### FIND MAXIMA AND MINIMA #########
	# Find the index of the max FPS_scaled value for each region_id
	idx_max = input_df.groupby('region_id', observed=True)['FPS_scaled'].idxmax()
	idx_min = input_df.groupby('region_id', observed=True)['FPS_scaled'].idxmin()
	# Select the corresponding rows
	max_fps_scaled = input_df.loc[idx_max]
	min_fps_scaled = input_df.loc[idx_min]
	# do the same for AF
	idx_max_af = input_df.groupby('region_id', observed=True)['AF'].idxmax()
	idx_min_af = input_df.groupby('region_id', observed=True)['AF'].idxmin()
	# Select the corresponding rows
	max_af_raw = input_df.loc[idx_max_af]
	min_af_raw = input_df.loc[idx_min_af]
	# create masks for the inverse of max and min values
	mask = ~input_df.index.isin(idx_max)
	max_fps_scaled_inv = input_df[mask]
	mask = ~input_df.index.isin(idx_min)
	min_fps_scaled_inv = input_df[mask]
	# do the same for AF
	mask = ~input_df.index.isin(idx_max_af)
	max_af_raw_inv = input_df[mask]
	mask = ~input_df.index.isin(idx_min_af)
	min_af_raw_inv = input_df[mask]
	return max_fps_scaled, min_fps_scaled, max_af_raw, min_af_raw, max_fps_scaled_inv, min_fps_scaled_inv, max_af_raw_inv, min_af_raw_inv

def boxplot_region_filtsorted(input_df, motif_id, output_path, palette_col, threshold, central_stat=None):
	################# START box plot of AF distr on filtered sorted sites as well as the substype-hued stripplot ################
	# check existence of output directory
	if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-IQR.pdf') or not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-gbl-mean.pdf'):
		print(f'Plotting {motif_id} boxplot of AF distribution on filtered sites sorted by FPS variance...')
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=input_df, hue='sample_id', palette=palette_col, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF per site', width=15)
		plt.ylabel(ylabel, fontsize=10)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)

		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF variance', width=15)
		plt.ylabel(ylabel, fontsize=10)
		
		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=input_df, hue='sample_id', palette=palette_col, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		# plot horizontal line at fps_scaled_global_mean
		if threshold == 'central':
			plt.axhline(y=central_stat, color='black', linestyle='--')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS per site', width=15)
		plt.ylabel(ylabel, fontsize=10)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=4)
		if threshold == 'iqr':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=10)
		elif threshold == 'central':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and scaled FPS above global mean)', fontsize=10)
		plt.ylabel('Scaled FPS variance', fontsize=10)
		plt.subplots_adjust(hspace=0.05)
		# save the plot
		if threshold == 'iqr':
			print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-IQR.pdf', dpi=300, bbox_inches="tight")
			print('Plot saved.')
			# close the plot
			plt.close('all')
			print('Plot space closed.')
		elif threshold == 'central':
			print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-gbl-mean.pdf', dpi=300, bbox_inches="tight")
			print('Plot saved.')
			# close the plot
			plt.close('all')
			print('Plot space closed.')
	else:
		print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
	################# END box plot of AF distr on filtered sorted sites as well as the substype-hued stripplot ################

def boxplot_maxima(input_df, max_af, max_af_inv, max_fps, max_fps_inv, motif_id, output_path, palette_col, gray_palette):
	###################################### START box plot of AF distr on filtered sorted sites (with maxima) ################
	if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima.pdf'):
		print(f'Plotting {motif_id} boxplot of AF distribution on filtered sorted sites, highlighting maxima...')
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=max_af_inv, hue='sample_id', palette=gray_palette, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=max_af, hue='sample_id', palette=palette_col, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF per site (maxima)', width=15)
		plt.ylabel(ylabel, fontsize=10)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)

		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF variance', width=15)
		plt.ylabel(ylabel, fontsize=10)

		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps_inv, hue='sample_id', palette=gray_palette, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps, hue='sample_id', palette=palette_col, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS per site (maxima)', width=15)
		plt.ylabel(ylabel, fontsize=10)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=4)
		plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=10)
		plt.ylabel('Scaled FPS variance', fontsize=10)
		plt.subplots_adjust(hspace=0.05)
		# save the plot
		print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
		plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima.pdf', dpi=300, bbox_inches="tight")
		print('Plot saved.')
		# close the plot
		plt.close('all')
		print('Plot space closed.')
	else:
		print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
		################### END box plot of AF distr on filtered sorted sites (with maxima) ################

def boxplot_minima(input_df, min_af, min_af_inv, min_fps, min_fps_inv, motif_id, output_path, palette_col, gray_palette):
	################### START box plot of AF distr on filtered sorted sites (with minima) ################
	if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-minima.pdf'):
		print(f'Plotting {motif_id} boxplot of AF distribution on filtered sorted sites, highlighting minima...')
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=min_af_inv, hue='sample_id', palette=gray_palette, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=min_af, hue='sample_id', palette=palette_col, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF per site (minima)', width=15)
		plt.ylabel(ylabel, fontsize=10)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)

		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('AF variance', width=15)
		plt.ylabel(ylabel, fontsize=10)

		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=input_df, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps_inv, hue='sample_id', palette=gray_palette, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps, hue='sample_id', palette=palette_col, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(ticks=plt.xticks()[0], labels=[])
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS per site (maxima)', width=15)
		plt.ylabel(ylabel, fontsize=10)

		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=input_df, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=4)
		plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=10)
		plt.ylabel('Scaled FPS variance', fontsize=10)
		plt.subplots_adjust(hspace=0.05)
		# save the plot
		print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
		plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima.pdf', dpi=300, bbox_inches="tight")
		print('Plot saved.')
		# close the plot
		plt.close('all')
		print('Plot space closed.')
	else:
		print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
	#################### END box plot of AF distr on filtered sorted sites (with minima) ################

def process_data(target_file, output_path, threshold, plot=True):
	################ define color palettes ################
	dutchfield_colordict = {'S6R691V_her2': "#e60049", 'ANAB5F7_basal': "#0bb4ff", '98JKPD8_lumA': "#87bc45", 'PU24GB8_lumB': "#ef9b20", '2GAMBDQ_norm': "#b33dc6"}
	gray = 'lightgray'
	gray_colordict = {'S6R691V_her2': gray, 'ANAB5F7_basal': gray, '98JKPD8_lumA': gray, 'PU24GB8_lumB': gray, '2GAMBDQ_norm': gray}
	################ START ################
	print(f'Processing {target_file}...')
	# load the data
	matrix_afps, motif_id, afps_df_lpv = load_data(target_file)
	################ Create output directories ################
	# check existence of output directory
	os.makedirs(f'{output_path}/output-data/plots/{motif_id}', exist_ok=True)
	os.makedirs(f'{output_path}/output-data/tables/{motif_id}', exist_ok=True)
	# use MinMaxScaler to scale the raw fps values to range between 0 and 1
	print(f'Scaling {motif_id} FPS matrix...')
	fps_df_scaled, fps_df_scaled_lpv = scale_data(matrix_afps)
	# merge the two dataframes
	afps_full_dfl = afps_df_lpv.merge(fps_df_scaled_lpv, on=['region_id', 'sample_id'])
	print(f'{motif_id} matrix has been scaled and processed.')
	################ SAVEPOINT ################
	###########################################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable.tsv'):
		print(f'Saving {motif_id} FPS-scaled data table...')
		# save to file
		afps_full_dfl.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table already exists. Skipping...')
	###########################################
	################ SAVEPOINT ################
	# calculate variance statistics for AF and FPS values
	afps_stats_df = variance_stats(matrix_afps, motif_id, fps_df_scaled)
	#### merge now with afps_full_dfl
	# set the index to 'region_id'
	afps_full_dfli = afps_full_dfl.set_index('region_id')
	# merge afps_full_dfli with stats df on region_id
	afps_merged_stats = afps_full_dfli.merge(afps_stats_df, left_index=True, right_index=True, how='left')
	# sort naturally by region_id
	afps_stats_mergesorted =afps_merged_stats.reset_index().reindex(index=index_natsorted(afps_merged_stats.index))
	print(f'Stats calculated for {motif_id} FPS matrix.')
	################ SAVEPOINT ################
	###########################################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable_with_stats.tsv'):
		print(f'Saving {motif_id} FPS-scaled data table with stats...')
		# save file
		afps_stats_mergesorted.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable_with_stats.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table with stats already exists. Skipping...')
	# generate accessory plots
	filtered_df = accessory_df(afps_stats_mergesorted)
	print(f'Saving {motif_id} FPS_scaled and AF variances filtered for values more than 0.001...')
	# save file
	filtered_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_var_filtered_unique_regions.tsv', sep='\t', index=True)
	# plot scatter plot of AF_var vs FPS_scaled_var
	print(f'Plotting {motif_id} scatter plot of AF_var vs FPS_scaled_var...')
	accessory_plot(filtered_df, motif_id, output_path)
	###########################################
	################ SAVEPOINT ################
	
	############# BASIC FILTERING #############
	# as the first returned output would not be used here, assign it to _
	_, merged_filt_nozeroaf, atleast_one_zero_af = basic_filtering(afps_stats_mergesorted)
	################ SAVEPOINT ################
	###########################################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afzero_atleast_one_regions_datatable.tsv'):
		print(f'Saving data table of {motif_id} regions with at least one zero AF value per region ID...')
		# save file
		atleast_one_zero_af.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afzero_atleast_one_regions_datatable.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table where at least one zero AF value per region ID already exists. Skipping...')
	###########################################
	################ SAVEPOINT ################
	# copy the filtered dataframe
	mf_df = merged_filt_nozeroaf.copy()
	mf_df = mf_df.reset_index(drop=True)
	# to ensure that each unique region_id is retained as a group of subtype rows, we need to filter after grouping per region_id
	print(f'Thresholding {motif_id} processed matrix...')
	high_af = mf_df.groupby('region_id').filter(lambda x: (x['AF'] > 0.5).all())
	low_af = mf_df.groupby('region_id').filter(lambda x: (x['AF'] <= 0.5).all())
	######## THRESHOLDING ########
	if threshold == 'iqr':
		high_af_fps_outliers, *_ = thresholding_strat(high_af, 'iqr')
		high_af_fps_outliers_filtsorted = filtersort_df(high_af_fps_outliers)
		################ SAVEPOINT ################
		###########################################
		if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regs_abv_IQR_threshold_sorted_by_FPS_var_table.tsv'):
			print(f'Saving data table of {motif_id} regions with high AF (>0.5) passing IQR threshold and sorted by FPS variance...')
			# save file
			high_af_fps_outliers_filtsorted.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regs_abv_IQR_threshold_sorted_by_FPS_var_table.tsv', sep='\t', index=False)
		else:
			print(f'{motif_id} data table of regions passing IQR threshold already exists. Skipping...')
		###########################################
		# find maxima and minima
		max_fps_scaled, min_fps_scaled, max_af_raw, min_af_raw, max_fps_scaled_inv, min_fps_scaled_inv, max_af_raw_inv, min_af_raw_inv = find_min_max(high_af_fps_outliers_filtsorted)
		print(f'{motif_id} processed matrix has been extensively filtered and sorted. Entering plotting phase...')

		if plot == True:
			boxplot_region_filtsorted(high_af_fps_outliers_filtsorted, motif_id, output_path, dutchfield_colordict, threshold)
			boxplot_maxima(high_af_fps_outliers_filtsorted, max_af_raw, max_af_raw_inv, max_fps_scaled, max_fps_scaled_inv, motif_id, output_path, dutchfield_colordict, gray_colordict)
			boxplot_minima(high_af_fps_outliers_filtsorted, min_af_raw, min_af_raw_inv, min_fps_scaled, min_fps_scaled_inv, motif_id, output_path, dutchfield_colordict, gray_colordict)
			# close all plots
			plt.close('all')
		else:
			print('Skipping plotting phase...')
		
		print('Now quantifying the number of sites where a subtype has maximum FPS_scaled value...')
		# quantify the number of filtered sites per sample_id
		max_af_raw_subset = max_af_raw[['region_id', 'sample_id']]
		max_af_raw_df = max_af_raw_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
		max_fps_scaled_subset = max_fps_scaled[['region_id', 'sample_id']]
		max_fps_scaled_df = max_fps_scaled_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
		max_af_raw_subset.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_max_af_region-ids_unique.tsv', sep='\t', index=False)
		max_af_raw_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_max_af_regions_count.tsv', sep='\t', index=True)
		max_fps_scaled_subset.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_max_fps-scaled_region-ids_unique.tsv', sep='\t', index=False)
		max_fps_scaled_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_max_fps-scaled_regions_count.tsv', sep='\t', index=True)
		print('Region counts have been performed and saved to files.')
		print('Data analysis complete. Exiting...')

	elif threshold == 'central':
		high_af_abovemean, high_af_belowmean, low_af_abovemean, low_af_belowmean = thresholding_strat(high_af, 'central', low_af, mf_df)
		high_af_abvmean_fs = filtersort_df(high_af_abovemean)
		high_af_blwmean_fs = filtersort_df(high_af_belowmean)
		low_af_abvmean_fs = filtersort_df(low_af_abovemean)
		low_af_blwmean_fs = filtersort_df(low_af_belowmean)
		###########################################
		################ SAVEPOINT ################
		if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_HI_AF_regs_abv_FPS-mean_sorted_by_FPS_var_table.tsv'):
			print(f'Saving data table of {motif_id} regions above or below mean scaled FPS and sorted by FPS variance...')
			# save files
			high_af_abvmean_fs.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_HI_AF_regs_GT_FPS-mean_sorted_by_FPS_var_table.tsv', sep='\t', index=False)
			# quantify the number of filtered sites per sample_id
			high_af_abvmean_fs_subset = high_af_abvmean_fs[['region_id', 'sample_id']]
			high_af_abvmean_fs_subset_df = high_af_abvmean_fs_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
			high_af_abvmean_fs_subset_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_HI-AF_GT-mean_regions_count.tsv', sep='\t', index=True)

			high_af_blwmean_fs.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_HI_AF_regs_LT_FPS-mean_sorted_by_FPS_var_table.tsv', sep='\t', index=False)
			# quantify the number of filtered sites per sample_id
			high_af_blwmean_fs_subset = high_af_blwmean_fs[['region_id', 'sample_id']]
			high_af_blwmean_fs_subset_df = high_af_blwmean_fs_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
			high_af_blwmean_fs_subset_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_HI-AF_LT-mean_regions_count.tsv', sep='\t', index=True)

			low_af_abvmean_fs.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_LO_AF_regs_GT_FPS-mean_sorted_by_FPS_var_table.tsv', sep='\t', index=False)
			# quantify the number of filtered sites per sample_id
			low_af_abvmean_fs_subset = low_af_abvmean_fs[['region_id', 'sample_id']]
			low_af_abvmean_fs_subset_df = low_af_abvmean_fs_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
			low_af_abvmean_fs_subset_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_LO-AF_GT-mean_regions_count.tsv', sep='\t', index=True)

			low_af_blwmean_fs.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_LO_AF_regs_LT_FPS-mean_sorted_by_FPS_var_table.tsv', sep='\t', index=False)
			# quantify the number of filtered sites per sample_id
			low_af_blwmean_fs_subset = low_af_blwmean_fs[['region_id', 'sample_id']]
			low_af_blwmean_fs_subset_df = low_af_blwmean_fs_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
			low_af_blwmean_fs_subset_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_LO-AF_LT-mean_regions_count.tsv', sep='\t', index=True)
			print(f'All filtered data tables have been saved to file.')
			print('Data analysis complete. Exiting...')
		else:
			print(f'{motif_id} data tables of regions above or below FPS mean already exist. Skipping...')


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
		process_data(target_file, output_dir, 'iqr', True)

	# uncomment this to run in parallel
	# with cf.ProcessPoolExecutor(max_workers=8) as executor:
	# 	executor.map(process_data, inputs, it.repeat(output_dir), it.repeat('central'), it.repeat(True))

	print ("Pipeline finished! All footprint matrices have been processed.")
