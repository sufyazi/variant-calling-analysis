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

def thresholding_strat(input_df, threshold, unsplit_df=None):
	if threshold == 'iqr':
		# thresholding strategy 1: use IQR method to capture outliers and return only outlier regions
		# subset high_af for sites with FPS_var > 75th percentile + 1.5 * IQR
		iqr = input_df['FPS_scaled_var'].quantile(q=0.75) - input_df['FPS_scaled_var'].quantile(q=0.25)
		high_af_outliers = input_df[input_df['FPS_scaled_var'] > input_df['FPS_scaled_var'].quantile(q=0.75) + (1.5 * iqr)]
		return high_af_outliers
	elif threshold == 'median':
		# thresholding strategy 2: compute median of FPS_scaled data points from the unsplit (high+low af) dataframe and return only regions with FPS_scaled > median
		global_median = unsplit_df['FPS_scaled'].median()
		high_af_abovemedian = input_df[input_df['FPS_scaled'] > global_median]
		return high_af_abovemedian
	else:
		print('Invalid thresholding strategy. Please choose either "iqr" or "median".')
		sys.exit(1)

def process_data(target_file, output_path, plot=True, threshold='iqr'):
	print(f'Processing {target_file}...')
	# import the data
	matrix_afps = pd.read_csv(target_file, sep='\t')
	# extract motif id from filename
	motif_id = os.path.basename(target_file).replace('_fpscore-af-varsites-combined-matrix-wide.tsv', '')
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

	# use MinMaxScaler to scale the raw fps values to range between 0 and 1
	# scale the FPS values to a range of 0-1
	# Initialize a MinMaxScaler
	print(f'Scaling {motif_id} FPS matrix...')
	scaler = MinMaxScaler()
	# copy df
	fps_df_scaled = matrix_afps.filter(regex='_fps$|_id$').copy()
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
	# merge the two dataframes
	afps_full_dfl = afps_df_lpv.merge(fps_df_scaled_lpv, on=['region_id', 'sample_id'])
	print(f'{motif_id} matrix has been scaled and processed.')
	################ SAVEPOINT ################
	os.makedirs(f'{output_path}/output-data/tables/{motif_id}', exist_ok=True)
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable.tsv'):
		print(f'Saving {motif_id} FPS-scaled data table...')
		# save to file
		afps_full_dfl.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table already exists. Skipping...')
	################ SAVEPOINT ################

	# calculate the variance of AF and FPS values across samples per region_id
	# extract fps columns
	fps_df = matrix_afps.filter(regex='_fps$|_id$').copy()
	print(f'Calculating {motif_id} FPS and AF variance statistics...')
	# calculate variance of fps values across samples per region_id and add to a new column called 'fps_var'
	fps_df = fps_df.set_index('region_id')
	fps_df['FPS_var'] = fps_df.var(axis=1)
	# now calculate variance of fps_scaled values across samples per region_id and add to a new column called 'fps_scaled_var'
	fps_df_scaled_cp = fps_df_scaled.copy()
	fps_df_scaled['FPS_scaled_var'] = fps_df_scaled_cp.var(axis=1)
	# Do the same for AF values
	# extract af columns
	af_df = matrix_afps.filter(regex='_AF$|_id$')
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
	
	#### merge now with afps_full_dfl
	# set the index to 'region_id'
	afps_full_dfli = afps_full_dfl.set_index('region_id')
	# merge afps_full_dfli with stats df on region_id
	afps_merged_stats = afps_full_dfli.merge(afps_stats_df, left_index=True, right_index=True, how='left')
	# sort naturally by region_id
	afps_stats_mergesorted =afps_merged_stats.reset_index().reindex(index=index_natsorted(afps_merged_stats.index))
	print(f'Stats calculated for {motif_id} FPS matrix.')
	################ SAVEPOINT ################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable_with_stats.tsv'):
		print(f'Saving {motif_id} FPS-scaled data table with stats...')
		# save file
		afps_stats_mergesorted.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_fullscaled_longtable_with_stats.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table with stats already exists. Skipping...')
	################ SAVEPOINT ################
	### ASIDE
	# subset a copy of the dataframe with only the AF_var and FPS_scaled_var columns and region_id
	afps_stats_ms_var = afps_stats_mergesorted[['region_id', 'FPS_scaled_var', 'AF_var']].copy().drop_duplicates().reset_index(drop=True)
	# set region_id as index
	afps_stats_ms_var = afps_stats_ms_var.set_index('region_id')
	# filter rows using the condition FPS_scaled_var > 0.001 and AF_var > 0.001
	test_df = afps_stats_ms_var[(afps_stats_ms_var['FPS_scaled_var'] > 0.001) & (afps_stats_ms_var['AF_var'] > 0.001)]
	### END ASIDE
	################ SAVEPOINT ################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_var_filtered_unique_regions.tsv'):
		print(f'Saving {motif_id} FPS_scaled and AF variances filtered for values more than 0.001...')
		# save file
		test_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afps_var_filtered_unique_regions.tsv', sep='\t', index=True)
	else:
		print(f'{motif_id} data table of AF and FPS scaled variances already exists. Skipping...')
	################ SAVEPOINT ################
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
	# also do the inverse filtering and save as a separate dataframe
	atleast_one_zero_af = merged_filt_nozero[merged_filt_nozero['region_id'].isin(max_af['region_id'])]
	################ SAVEPOINT ################
	if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afzero_atleast_one_regions_datatable.tsv'):
		print(f'Saving data table of {motif_id} regions with at least one zero AF value per region ID...')
		# save file
		atleast_one_zero_af.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_afzero_atleast_one_regions_datatable.tsv', sep='\t', index=False)
	else:
		print(f'{motif_id} data table where at least one zero AF value per region ID already exists. Skipping...')
	################ SAVEPOINT ################
	######## FILTERING ########
	# copy the filtered dataframe
	mf_df = merged_filt_nozeroaf.copy()
	mf_df = mf_df.reset_index(drop=True)
	# calculate AF median per region_id
	mf_df['AF_median'] = mf_df.groupby('region_id')['AF'].transform('median')
	# to ensure that each unique region_id is retained as a group of subtype rows, we need to filter after grouping per region_id
	print(f'Thresholding {motif_id} processed matrix...')
	high_af = mf_df.groupby('region_id').filter(lambda x: (x['AF'] > 0.5).all())
	######## THRESHOLDING ########
	if threshold == 'iqr':
		high_af_fps_outliers = thresholding_strat(high_af, 'iqr')
	elif threshold == 'median':
		high_af_fps_outliers = thresholding_strat(high_af, 'median', mf_df)
	# discard sites if it has AF_var == 0
	high_nzaf_outliers = high_af_fps_outliers[high_af_fps_outliers['AF_var'] != 0]
	# extract unique region_ids from the high_nzaf_outliers df
	high_af_uniq_reg = high_nzaf_outliers[['region_id', 'AF_var', 'FPS_scaled_var']].drop_duplicates()
	# sort the region_id in the filtered high_af dataframe based on descending order of AF_var or FPS_scaled_var
	# afvar_uniqsort = high_af_uniq_reg.sort_values(by='AF_var', ascending=False)
	fpsvar_uniqsort = high_af_uniq_reg.sort_values(by='FPS_scaled_var', ascending=False)
	# extract the region_id from the sorted dataframe
	sorted_region_ids = fpsvar_uniqsort['region_id']
	# Change 'region_id' to a categorical variable with the categories ordered by 'fps_sorted_region_ids'
	df_copy = high_nzaf_outliers.copy()
	df_copy['region_id'] = pd.Categorical(high_nzaf_outliers['region_id'], categories=sorted_region_ids, ordered=True)
	# Filter the DataFrame
	high_nzaf_outlie_filt = df_copy[df_copy['region_id'].isin(sorted_region_ids)]
	# Sort the DataFrame by 'region_id'
	high_nzaf_outlie_filt = high_nzaf_outlie_filt.sort_values('region_id')
	# get unique sample_id values into a list to define a categorical order
	datasets = high_nzaf_outlie_filt['sample_id'].unique().tolist()
	datasets = sorted(datasets)
	# Create a categorical variable with ordered categories
	dataset_copy = high_nzaf_outlie_filt.copy()
	dataset_copy['sample_id'] = pd.Categorical(dataset_copy['sample_id'], categories=datasets, ordered=True)
	# Sort 'sample_id' within each 'region_id'
	high_nzaf_outlie_filtsort = dataset_copy.groupby('region_id', sort=False, observed=False).apply(lambda x: x.sort_values('sample_id')).reset_index(drop=True)
	################ SAVEPOINT ################
	if threshold == 'iqr':
		if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regions_passing_IQR_threshold_sorted_by_FPS_var_datatable.tsv'):
			print(f'Saving data table of {motif_id} regions with high AF (>0.5) passing IQR threshold and sorted by FPS variance...')
			# save file
			high_nzaf_outlie_filtsort.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regions_passing_IQR_threshold_sorted_by_FPS_var_datatable.tsv', sep='\t', index=False)
		else:
			print(f'{motif_id} data table of regions passing IQR threshold already exists. Skipping...')
	elif threshold == 'median':
		if not os.path.exists(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regions_above_median_FPS_sorted_by_FPS_var_datatable.tsv'):
			print(f'Saving data table of {motif_id} regions with high AF (>0.5) above median FPS and sorted by FPS variance...')
			# save file
			high_nzaf_outlie_filtsort.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_high_AF_regions_above_median_FPS_sorted_by_FPS_var_datatable.tsv', sep='\t', index=False)
		else:
			print(f'{motif_id} data table of regions above median FPS already exists. Skipping...')
	################ SAVEPOINT ################
	######### FIND MAXIMA AND MINIMA #########
	# Find the index of the max FPS_scaled value for each region_id
	idx_max = high_nzaf_outlie_filtsort.groupby('region_id', observed=True)['FPS_scaled'].idxmax()
	idx_min = high_nzaf_outlie_filtsort.groupby('region_id', observed=True)['FPS_scaled'].idxmin()
	# Select the corresponding rows
	max_fps_scaled = high_nzaf_outlie_filtsort.loc[idx_max]
	min_fps_scaled = high_nzaf_outlie_filtsort.loc[idx_min]
	# do the same for AF
	idx_max_af = high_nzaf_outlie_filtsort.groupby('region_id', observed=True)['AF'].idxmax()
	idx_min_af = high_nzaf_outlie_filtsort.groupby('region_id', observed=True)['AF'].idxmin()
	# Select the corresponding rows
	max_af_raw = high_nzaf_outlie_filtsort.loc[idx_max_af]
	min_af_raw = high_nzaf_outlie_filtsort.loc[idx_min_af]
	# create masks for the inverse of max and min values
	mask = ~high_nzaf_outlie_filtsort.index.isin(idx_max)
	max_fps_scaled_inv = high_nzaf_outlie_filtsort[mask]
	mask = ~high_nzaf_outlie_filtsort.index.isin(idx_min)
	min_fps_scaled_inv = high_nzaf_outlie_filtsort[mask]
	# do the same for AF
	mask = ~high_nzaf_outlie_filtsort.index.isin(idx_max_af)
	max_af_raw_inv = high_nzaf_outlie_filtsort[mask]
	mask = ~high_nzaf_outlie_filtsort.index.isin(idx_min_af)
	min_af_raw_inv = high_nzaf_outlie_filtsort[mask]
	print(f'{motif_id} processed matrix has been extensively filtered and sorted. Entering plotting phase...')
	if plot == True:
		################ define color palettes ################
		springpastel = ["#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a"]
		dutchfield_colordict = {'S6R691V_her2': "#e60049", 'ANAB5F7_basal': "#0bb4ff", '98JKPD8_lumA': "#87bc45", 'PU24GB8_lumB': "#ef9b20", '2GAMBDQ_norm': "#b33dc6"}
		gray = 'lightgray'
		gray_colordict = {'S6R691V_her2': gray, 'ANAB5F7_basal': gray, '98JKPD8_lumA': gray, 'PU24GB8_lumB': gray, '2GAMBDQ_norm': gray}

		################ PLOT ASIDE ################
		# plot scatter plot of AF_var vs FPS_scaled_var unfiltered
		print(f'Plotting {motif_id} scatter plot of AF_var vs FPS_scaled_var (unfiltered)...')
		g = sns.jointplot(data=afps_stats_ms_var, x="FPS_scaled_var", y="AF_var", height=10, ratio=5, color='darkslateblue' )
		g.fig.suptitle(f'AF vs FPS var for {motif_id} (unfiltered: {len(afps_stats_ms_var)} regions)')
		# check existence of output directory
		os.makedirs(f'{output_path}/output-data/plots/{motif_id}', exist_ok=True)
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-unfilt.pdf'):
			print(f'Saving {motif_id} scatter plot of AF and scaled FPS variance (unfiltered)...')
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-unfilt.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
			print('Plot space closed.')
		else:
			print(f'{motif_id} unfiltered variance scatterplot already exists. Skipping...')
			plt.close()
		# now plot the filtered scatterplot
		print(f'Plotting {motif_id} scatter plot of AF_var vs FPS_scaled_var (filtered)...')
		g = sns.jointplot(data=test_df, x="FPS_scaled_var", y="AF_var", height=10, ratio=5, color='darkslateblue' )
		g.fig.suptitle(f'AF vs FPS var for {motif_id} (filtered: {len(test_df)} regions)')
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-filt.pdf'):
			print(f'Saving {motif_id} scatter plot of AF and scaled FPS variance (filtered)...')
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AF_vs_FPS-scaled_variance_scatterplot-filt.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
			print('Plot space closed.')
		else:
			print(f'{motif_id} filtered variance scatterplot already exists. Skipping...')
			plt.close()
		################ END PLOT ASIDE ################

  		################ START box plot of high AF regions (>0.5) ################
		print(f'Plotting {motif_id} boxplot for sites with AF > 0.5...')
		plt.figure(figsize=(10, 10), dpi=300)
		sns.boxplot(x='region_id', y='AF', data=high_af, color='lightgray', linecolor='black', linewidth=1.5, showfliers=False)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='AF', data=high_af, hue='sample_id', size=6, jitter=True, palette=springpastel, linewidth=0.5, edgecolor='black')
		# plot legend outside of the plot
		plt.legend(bbox_to_anchor=(1.225, 1),borderaxespad=0, markerscale=1, fontsize=10)
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (with AF > 0.5 in all subtypes)', fontsize=10)
		plt.ylabel('Allele frequency (AF)', fontsize=8)
		
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_sites_AF_morethan_0.5-boxplot.pdf'):
			print(f'Saving {motif_id} boxplot for sites with AF > 0.5...')
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_sites_AF_morethan_0.5-boxplot.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
			print('Plot space closed.')
		else:
			print(f'{motif_id} filtered high AF boxplot already exists. Skipping...')
			plt.close()
		################ END box plot of high AF regions (>0.5) ################
		
		################ START box plot of high AF set across raw AF and scaled FPS data ################
		print(f'Plotting {motif_id} boxplots for sites with AF > 0.5 with FPS scaled data...')
		plt.figure(figsize=(10, 10), dpi=300)
		plt.subplot(2, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_af, color='lightgray', linecolor='black', linewidth=1.5, showfliers=False)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='AF', data=high_af, hue='sample_id', size=6, jitter=True, palette=springpastel, linewidth=0.5, edgecolor='black')
		# plot legend outside of the plot
		plt.legend(bbox_to_anchor=(1.225, 1),borderaxespad=0, markerscale=1, fontsize=10)
		plt.xlabel('')
		plt.ylabel('Allele frequency (AF)', fontsize=8)
		# then plot the scaled FPS values for these sites
		plt.subplot(2, 1, 2)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_af, color='white', linecolor='black', linewidth=1.5, showfliers=False)
		plt.xticks(rotation=90, fontsize=6)
		sns.stripplot(x='region_id', y='FPS_scaled', data=high_af, hue='sample_id', size=6, jitter=True, palette=springpastel, legend=False, linewidth=0.5, edgecolor='black')
		plt.xlabel(f'{motif_id} TF binding sites with allelic variants (AF > 0.5)', fontsize=10)
		plt.ylabel(f'Footprint scores (FPS) (min-max scaled)', fontsize=8)
		plt.subplots_adjust(hspace=0.6)
		
		if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_sites_AF_morethan_0.5-with-FPS-scaled-boxplot.pdf'):
			print(f'Saving {motif_id} boxplot for sites with AF > 0.5...')
			# save the plot
			plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_sites_AF_morethan_0.5-with-FPS-scaled-boxplot.pdf', dpi=300, bbox_inches="tight")
			# close the plot
			plt.close()
			print('Plot space closed.')
		else:
			print(f'{motif_id} filtered high AF boxplot with FPS data already exists. Skipping...')
			plt.close()
		################ END box plot of high AF set across raw AF and scaled FPS data ################

		################# START box plot of AF distr on filtered sorted sites as well as the substype-hued stripplot ################
		print(f'Plotting {motif_id} boxplot of AF distribution on filtered sites sorted by FPS variance...')
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=high_nzaf_outlie_filtsort, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		plt.ylabel('AF Distribution Per Site', fontsize=10)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)
		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=10)
		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=high_nzaf_outlie_filtsort, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		plt.ylabel('FPS Scaled Distribution Per Site', fontsize=10)
		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		if threshold == 'iqr':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=10)
		elif threshold == 'median':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and scaled FPS above median)', fontsize=10)
		plt.ylabel('FPS Variance', fontsize=10)
		plt.subplots_adjust(hspace=0.6)

		# check existence of output directory
		if threshold == 'iqr':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-IQR.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-IQR.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		elif threshold == 'median':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-FPSmedian.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_by_FPS_var_boxplot-FPSmedian.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		################# END box plot of AF distr on filtered sorted sites as well as the substype-hued stripplot ################

		###################################### START box plot of AF distr on filtered sorted sites (with maxima) ################
		print(f'Plotting {motif_id} boxplot of AF distribution on filtered sorted sites, highlighting maxima...')
		import textwrap
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=max_af_raw_inv, hue='sample_id', palette=gray_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=max_af_raw, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		ylabel = textwrap.fill('AF Per Site (maxima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)
		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=8)
		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps_scaled_inv, hue='sample_id', palette=gray_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=max_fps_scaled, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS Per Site (maxima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		if threshold == 'iqr':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=8)
		elif threshold == 'median':
			plt.xlabel(f'{motif_id} binding sites with allelic variants (AF > 0.5 and scaled FPS above median)', fontsize=8)
		plt.ylabel('FPS Variance', fontsize=8)
		plt.subplots_adjust(hspace=0.6)
		# save the plot
		# check existence of output directory
		if threshold == 'iqr':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima-IQR.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima-IQR.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		elif threshold == 'median':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima-FPSmedian.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-maxima-FPSmedian.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		################### END box plot of AF distr on filtered sorted sites (with maxima) ################

		################### START box plot of AF distr on filtered sorted sites (with minima) ################
		import textwrap
		plt.figure(figsize=(10, 10), dpi=300)
		# specify subplot
		plt.subplot(4, 1, 1)
		sns.boxplot(x='region_id', y='AF', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='AF', data=min_af_raw_inv, hue='sample_id', palette=gray_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='AF', data=min_af_raw, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		ylabel = textwrap.fill('AF Per Site (minima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		# place legend outside of the plot
		plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0, markerscale=2, fontsize=10)
		plt.subplot(4, 1, 2)
		sns.barplot(x='region_id', y='AF_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		plt.ylabel('AF Variance', fontsize=8)
		plt.subplot(4, 1, 3)
		sns.boxplot(x='region_id', y='FPS_scaled', data=high_nzaf_outlie_filtsort, color='whitesmoke', linecolor='black', showfliers=False)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps_scaled_inv, hue='sample_id', palette=gray_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='dimgray', alpha=0.8)
		sns.stripplot(x='region_id', y='FPS_scaled', data=min_fps_scaled, hue='sample_id', palette=dutchfield_colordict, size=4, jitter=True, legend=False, linewidth=0.5, edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		plt.xlabel('')
		ylabel = textwrap.fill('Scaled FPS Per Site (minima)', width=20)
		plt.ylabel(ylabel, fontsize=8)
		plt.subplot(4, 1, 4)
		sns.barplot(x='region_id', y='FPS_scaled_var', data=high_nzaf_outlie_filtsort, color='darkslateblue', edgecolor='black')
		plt.xticks(rotation=90, fontsize=3)
		if threshold == 'iqr':
			plt.xlabel(f'{motif_id} TF binding sites with allelic variants (AF > 0.5 and passing IQR threshold)', fontsize=8)
		elif threshold == 'median':
			plt.xlabel(f'{motif_id} TF binding sites with allelic variants (AF > 0.5 and scaled FPS above median)', fontsize=8)
		plt.ylabel('FPS Variance', fontsize=8)
		plt.subplots_adjust(hspace=0.6)
		# save the plot
		# check existence of output directory
		if threshold == 'iqr':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-minima-IQR.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-minima-IQR.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		elif threshold == 'median':
			if not os.path.exists(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-minima-FPSmedian.pdf'):
				print(f'Saving {motif_id} boxplot of AF distribution on filtered sorted sites...')
				# save the plot
				plt.savefig(f'{output_path}/output-data/plots/{motif_id}/{motif_id}_AFdist_per_site_AFvar_filtsorted_with_FPS_boxplot-minima-FPSmedian.pdf', dpi=300, bbox_inches="tight")
				# close the plot
				plt.close()
				print('Plot space closed.')
			else:
				print(f'{motif_id} box plots of AF distribution on filtered sites already exists. Skipping...')
				plt.close()
		#################### END box plot of AF distr on filtered sorted sites (with minima) ################
	else:
		print('Skipping plotting phase...')

	if threshold == 'iqr':
		print('Now quantifying the number of sites where a subtype has maximum FPS_scaled value...')
		# quantify the number of filtered sites per sample_id
		max_af_raw_subset = max_af_raw[['region_id', 'sample_id']]
		max_af_raw_df = max_af_raw_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
		max_fps_scaled_subset = max_fps_scaled[['region_id', 'sample_id']]
		max_fps_scaled_df = max_fps_scaled_subset.groupby('sample_id', observed=False)['region_id'].count().to_frame()
		max_af_raw_subset.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_af_region-ids_unique.tsv', sep='\t', index=False)
		max_af_raw_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_af_regions_count.tsv', sep='\t', index=True)
		max_fps_scaled_subset.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_fps-scaled_region-ids_unique.tsv', sep='\t', index=False)
		max_fps_scaled_df.to_csv(f'{output_path}/output-data/tables/{motif_id}/{motif_id}_maximum_fps-scaled_regions_count.tsv', sep='\t', index=True)
		print('Region counts have been performed and saved to files.')
		print('Data analysis complete. Exiting...')
	elif threshold == 'median':
		print('Data analysis complete. Exiting...')


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
	#for target_file in inputs:
	#	process_data(target_file, output_dir, True)

	# uncomment this to run in parallel
	with cf.ProcessPoolExecutor(max_workers=8) as executor:
		executor.map(process_data, inputs, it.repeat(output_dir), it.repeat(True), it.repeat('iqr'))

	print ("Pipeline finished! All footprint matrices have been processed.")