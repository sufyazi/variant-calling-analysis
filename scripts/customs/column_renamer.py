#!/usr/bin/env python3

import sys
import os
import polars as pl

# Get the path of the output noheader vcf file (.txt)
path = sys.argv[1]

def extract_column_names(path):
    """
    Extracts the column names from a csv file
    """
    df = pl.read_csv(path, separator="\t")

    column_names = df.columns
    column_names = [col for col in column_names if "aliquot" in col]
    
    return column_names

def rename_columns(path):
    