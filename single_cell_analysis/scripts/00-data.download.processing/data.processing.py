# Qiwen Hu 2018
# Data processing by z-score and zeroone norm
# Usage: python data.processing.py -i input_dir -f input_file -o out_dir
#

import os
import requests
import numpy as np
import pandas as pd
import argparse
from sklearn import preprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--input_file',
                    help='input_file')
parser.add_argument('-i', '--input_dir',
                    help='input dir')
parser.add_argument('-o', '--output_dir',
                    help='output dir')

args = parser.parse_args()

input_dir = args.input_dir
input_file = args.input_file
output_dir = args.output_dir

# Processing RNAseq data by z-score and zeroone norm
rna_out_zeroone_file = os.path.join(output_dir, str(input_file) + 'scaled_zeroone_rnaseq.tsv')

rnaseq_df = pd.read_table(os.path.join(input_dir, input_file), index_col=0)
rnaseq_df = rnaseq_df.fillna(0)
rnaseq_df = rnaseq_df.T

# Scale RNAseq data using zero-one normalization
rnaseq_scaled_zeroone_df = preprocessing.MinMaxScaler().fit_transform(rnaseq_df)
rnaseq_scaled_zeroone_df = pd.DataFrame(rnaseq_scaled_zeroone_df,
                                        columns=rnaseq_df.columns,
                                        index=rnaseq_df.index)
rnaseq_scaled_zeroone_df = rnaseq_scaled_zeroone_df.T
rnaseq_scaled_zeroone_df.to_csv(rna_out_zeroone_file, sep='\t')
