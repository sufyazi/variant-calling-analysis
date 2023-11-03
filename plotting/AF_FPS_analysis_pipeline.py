#!/#!/usr/bin/env python3

####################
# import libraries #
####################

import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

####################
# define functions #
####################

#############
# load args #
#############
af_fps_raw_df = pd.read_csv(sys.argv[1], sep='\t')
#############
# load data #
#############

# first, load AF data
af_region_df = af_fps_raw_df.filter(regex='_AF$|_id$')
af_df_indexed = af_region_df.set_index('region_id')
print(af_df_indexed.head())

# next, load FPS data
fps_region_df = af_fps_raw_df.filter(regex='_FPS$|_id$')
fps_df_indexed = fps_region_df.set_index('region_id')
print(fps_df_indexed.head())

######################

# next, convert the dataframes to numpy arrays for easier manipulation