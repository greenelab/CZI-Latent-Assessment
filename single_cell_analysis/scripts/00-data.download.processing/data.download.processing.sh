#!/bin/bash

# Pipeline to generate/download simulated and real single cell datasets and do data processing
#
# Usage: bash scripts/00-data.download.processing/data.download.processing.sh
#
# Output: log(or TPM) - normalized count files and files contain celltype information. 
#         files are located at data/real_data for real datasets and data/simulate_data for simulated datasets

# 1. Generate simulated datasets
  Rscript scripts/00-data.download.processing/splatter.simulation.R

# 2. download and process real datasets
  Rscript scripts/00-data.download.processing/real.datasets.download.R

# 3. zero-one normalization

# for simulated datasets
data_dir="data/simulate_data"

for file in `ls $data_dir/*count.matrix.txt`; do
	python data.processing.py -f ${file/*\//} -i $data_dir -o $data_dir
done

# for real datasets
data_dir="data/real_data"

for file in `ls $data_dir/*exp.matrix.txt`; do
	python scripts/00-data.download.processing/data.processing.py -f ${file/*\//} -i $data_dir -o $data_dir
done
