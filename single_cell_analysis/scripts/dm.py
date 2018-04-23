"""
Qiwen Hu 2018
scripts/dm.py

The script is used to perform dimension reduction for single cell data 
based on ZIFA (https://github.com/epierson9/ZIFA) and umap (https://github.com/lmcinnes/umap)

Usage:
  python scripts/dm.py -f input_file -m [ZIFA/umap]

Output:
  file contains features projected into 2 dimensional space	

"""
import sys
import os
import numpy as np
import pandas as pd
import random
from ZIFA import ZIFA
from ZIFA import block_ZIFA
import umap
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-f', '--input_file', help='input_file')
parser.add_argument('-m', '--method', help='perfom umap or ZIFA dimension reduction')
args = parser.parse_args()

input_file = args.input_file
method = args.method

# Set seed
np.random.seed(32)

# Load expression data
rnaseq_file = os.path.join('data', input_file)

rnaseq_df = pd.read_table(rnaseq_file, index_col = 0)
rnaseq_df = rnaseq_df.T
rnaseq_exp = rnaseq_df.as_matrix()

# Perform uMAP dimension reduction on expression data
if method == "umap":
    embedding = umap.UMAP(n_neighbors=10, min_dist=0.1, metric='correlation').fit_transform(rnaseq_exp)
    umap_out = pd.DataFrame(embedding, columns=['1', '2'])
    umap_out.index = rnaseq_df.index
    umap_out.index.name = 'id'
    umap_out_file = os.path.join('../features', 
                             input_file + '_rnaseq_umap_features.tsv')
    umap_out.to_csv(umap_out_file, sep='\t')


# Perform ZIFA dimension reduction on expression data
elif method == "ZIFA": 
    k = 2
    Zhat, params = block_ZIFA.fitModel(rnaseq_exp, k)
    zifa_out = pd.DataFrame(Zhat, columns=['1', '2'])
    zifa_out.index = rnaseq_df.index
    zifa_out.index.name = 'id'
    zifa_out_file = os.path.join('../features', 
                             input_file + '_rnaseq_ZIFA_features.tsv')
    zifa_out.to_csv(zifa_out_file, sep='\t')

