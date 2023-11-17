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

def process_data(target_file, output_path, plot=False):
    pass


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

##################
# define globals #
##################
output_path = '/scratch/users/ntu/suffiazi/outputs/tmp/plotting/'

if __name__ == '__main__':
	inputs = process_input_tsv(root_dir)
	# uncomment this to run serially
	# for target_file in inputs:
	# 	process_data(target_file, output_path)
 
	# uncomment this to run in parallel
	with cf.ProcessPoolExecutor(max_workers=8) as executor:
		executor.map(process_data, inputs, it.repeat(output_path), it.repeat(False))
	print ("Pipeline finished! All footprint matrices have had the site variances calculated.")