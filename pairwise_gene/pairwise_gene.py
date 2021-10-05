# ---------------------------------------------------------
# This software is Copyright (c) 2021 Oregon State University.
# All Rights Reserved. You may use, modify and redistribute
# this software under the terms of the GNU Affero General
# Public License version 3 published by the Free Software
# Foundation. If you need any other terms of use, please
# contact the author or advantage@oregonstate.edu.
# ---------------------------------------------------------
# pairwise_gene.py
# Source: https://github.com/johnjare
# Created by Jared Johnson
#
# Calculates pairwise gene differences between isolates in
# a supplied presence/absence matrix.
# ---------------------------------------------------------

import pandas as pd
import itertools
import argparse

###### Arguments
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('-i', dest='input', type=str, help='path to gene matrix')
parser.add_argument('-o', dest='out', default = "pair_genes.csv", type=str, help='path for output file')
args = parser.parse_args()

##### Import Data
# import gene matrix
df = pd.read_csv(args.input, sep = "\t")

##### Calculation
# count number of samples in matrix
samples: int = len(df.columns)-1
# calculate gene presence across all samples
df["sum"] = df.sum(axis=1)
# create new matrix with only accessory genes and drop the 'sum' column
df_acc = df[df["sum"] < samples]
df_acc = df_acc.drop("sum", axis = 1)

# extract sample names
names: list = list(df_acc.columns[1:len(df_acc.columns)])
# create list of sample pairs
pairs: list = list(itertools.combinations(names,2))
# create empty dataframe
df_pairs = pd.DataFrame(columns=['iso1', 'iso2', 'genes'])
# empty counter
count: int = 0
# calculate pairwise gene differences from list
for pair in pairs:
    count = count+1
    
    iso1: str = pair[0]
    iso2: str = pair[1]
    
    genes: str = str(df_acc[df_acc[iso1] != df_acc[iso2]].shape[0])
    
    new_genes: list = [iso1, iso2, genes]
    
    df_pairs.loc[count] = new_genes

##### Write Output
df_pairs.to_csv(args.out, index = False)
    
