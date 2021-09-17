#  accessory_tree.py
#  
#
#  Created by Jared Johnson on 8/8/21.
#  

import pandas as pd
import os

df = pd.read_csv("gene_presence_absence.Rtab", sep = "\t")

samples: int = len(df.columns)-1

df["sum"] = df.sum(axis=1)

df_acc = df[df["sum"] < samples]
df_acc = df_acc.drop("sum", axis = 1)

os.system("rm ./accessory_genes.binary")
os.system("touch ./accessory_genes.binary")

for sample in df_acc.columns[1:len(df_acc.columns)]:
    header: str = ">"+sample
    sequence: str = df_acc[sample].to_string(header=False, index=False).replace("\n", "")
    f = open("accessory_genes.binary", "a")
    f.write(header+"\n")
    f.write(sequence+"\n")
    f.close()

f = open("accessory_genes.binary", "r")
