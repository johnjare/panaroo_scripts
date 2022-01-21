# ---------------------------------------------------------
# This software is Copyright (c) 2021 Oregon State University.
# All Rights Reserved. You may use, modify and redistribute
# this software under the terms of the GNU Affero General
# Public License version 3 published by the Free Software
# Foundation. If you need any other terms of use, please
# contact the author or advantage@oregonstate.edu.
# ---------------------------------------------------------
# pair_diff.py
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
parser.add_argument('-g', dest='genes', type=str, help='path to gene matrix')
parser.add_argument('-s', dest='SNPs', type=str, help='path to SNP VCF file')
parser.add_argument('-o', dest='out', type=str, help='path for output file')
args = parser.parse_args()

##### Import Data
# import VCF file and reformat
df_snp = pd.read_csv(args.SNPs, sep = "\t", skiprows=3)
columns = list(range(9,len(df_snp.columns)))
df_snp = df_snp.iloc[:,columns]
# import gene matrix file and reformat
df_gene = pd.read_csv(args.genes, sep = "\t")
columns = list(range(1,len(df_gene.columns)))

##### Subset only variable regions
# count number of samples in matrix
samples: int = len(df_gene.columns)
# calculate gene/SNP presence across all samples
df_snp["sum"] = df_snp.sum(axis=1)
df_gene["sum"] = df_gene.sum(axis=1)
# create new matrix with only accessory genes/core SNPs and drop the 'sum' column
df_snp_acc = df_snp[df_snp["sum"] < samples]
df_snp_acc = df_snp_acc.drop("sum", axis = 1)
df_gene_acc = df_gene[df_gene["sum"] < samples]
df_gene_acc = df_gene_acc.drop("sum", axis = 1)

##### Calculate pairwise differences
# extract sample names
names: list = list(df_gene_acc.columns[1:len(df_gene_acc.columns)])
# create list of sample pairs
pairs: list = list(itertools.combinations(names,2))
# create empty dataframe
df_pairs = pd.DataFrame(columns=['iso1', 'iso2', 'SNPs', 'genes'])
# empty counter
count: int = 0
# calculate pairwise gene/SNP differences from list
for pair in pairs:
    count = count+1
    
    iso1: str = pair[0]
    iso2: str = pair[1]
    
    # counts number of rows after subsetting for dissimilar values between two isolates
    SNPs: str = str(df_snp_acc[df_snp_acc[iso1] != df_snp_acc[iso2]].shape[0])
    gene: str = str(df_gene_acc[df_gene_acc[iso1] != df_gene_acc[iso2]].shape[0])
    
    new_value: list = [iso1, iso2, SNPs, gene]
    
    df_pairs.loc[count] = new_value

##### Write Output
df_pairs.to_csv(args.out, index = False)
