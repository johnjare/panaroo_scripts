# ---------------------------------------------------------
# This software is Copyright (c) 2021 Oregon State University.
# All Rights Reserved. You may use, modify and redistribute
# this software under the terms of the GNU Affero General
# Public License version 3 published by the Free Software
# Foundation. If you need any other terms of use, please
# contact the author or advantage@oregonstate.edu.
# ---------------------------------------------------------
# accessory_binary_genes.py
# Source: https://github.com/johnjare
# Created by Jared Johnson
#
# Creates a binary alignment file of accessory genes from a
# gene presence/absence matrix. The resulting alignment file
# can be used to generate a phylogenetic tree.
# ---------------------------------------------------------

import pandas as pd
import argparse
import os

###### Arguments
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('-i', dest='input', type=str, help='path to gene matrix')
parser.add_argument('-o', dest='out', default = "accessory_binary_genes.binary", type=str, help='path for output file')
args = parser.parse_args()

###### Import Data
df = pd.read_csv(args.input, sep = "\t")

###### Determine Accessory Genome
# calculate number of samples
samples: int = len(df.columns)-1
# calculate gene presence across entire dataset
df["sum"] = df.sum(axis=1)
# subset only accessory genes - i.e., those not present in all samples
df_acc = df[df["sum"] < samples]
df_acc = df_acc.drop("sum", axis = 1)

###### Build and Write Output File
# create output file and remove any previous versions
os.system("rm ./"+args.out)
os.system("touch ./"+args.out)

# build file
for sample in df_acc.columns[1:len(df_acc.columns)]:
    header: str = ">"+sample
    sequence: str = df_acc[sample].to_string(header=False, index=False).replace("\n", "")
    f = open(args.out, "a")
    f.write(header+"\n")
    f.write(sequence+"\n")
    f.close()
# write file
f = open(args.out, "r")
