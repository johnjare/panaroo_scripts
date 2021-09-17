from typing import IO
import argparse
import pandas as pd
import io

###### Arguments
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('-i', dest='input', type=str, help='path to gene matrix')
parser.add_argument('-p1', dest='p1', type=str, help='path to isolate list for population 1')
parser.add_argument('-p2', dest='p2', type=str, default = None, help='path to isolate list for population 2; if empty includes all isolates not in p1')
parser.add_argument('-o1', dest='o1', default = "p1_genes", type=str, help='output for population 1')
parser.add_argument('-o2', dest='o2', default = "p2_genes", type=str, help='output for population 2')
parser.add_argument('-a', dest='accessory', type=bool, default = False, help='returns list of accessory genes if "True"')
args = parser.parse_args()

print(args.p1)

##### Data
# import gene matrix and population lists
df = pd.read_csv(args.input, sep = "\t")
df_p1 = pd.read_csv(args.p1, sep = "\t", header = None)
p1: list = list(df_p1[0])

# creates list of all isolates not in population 1 if the 'p2' flag is left empty
if args.p2 != None:
    df_p2 = pd.read_csv(args.p2, sep = "\t", header = None)
    p2: list = list(df_p2[0])
else:
    df_p2 = df.drop(p1, axis=1)
    p2: list = list(df_p2.drop('Gene', axis=1))
    
##### Calculation
# calculate gene presence in population 1
df["p1_sum"] = df[p1].sum(axis=1)
# calculate gene presence in population 2
df["p2_sum"] = df[p2].sum(axis=1)

# genes present in all members of population 1 but none of population 2
p1_genes = df[(df["p1_sum"] == len(p1)) & (df["p2_sum"] == 0)]["Gene"]
# genes present in all members of population 2 but none of population 1
p2_genes = df[(df["p2_sum"] == len(p2)) & (df["p1_sum"] == 0)]["Gene"]

##### Save gene lists
p1_genes.to_csv(args.o1+".lst", index = False, header = False)
p2_genes.to_csv(args.o2+".lst", index = False, header = False)

# returns list of accessory genes if requested
if args.accessory == True:
    samples: int = len(df.columns)-1

    df["sum"] = df.sum(axis=1)

    accessory_genes = df[df["sum"] < samples]["Gene"]
    
    accessory_genes.to_csv("accessory_genes.lst", index = False, header = False)
