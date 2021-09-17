import pandas as pd
import itertools

# import gene matrix
df = pd.read_csv("gene_presence_absence.Rtab", sep = "\t")
# count number of samples in matrix
samples: int = len(df.columns)-1
# calculate gene presence across all samples
df["sum"] = df.sum(axis=1)
# create new matrix with only accessory genes and drop the 'sum' column
df_acc = df[df["sum"] < samples]
df_acc = df_acc.drop("sum", axis = 1)

# extract sample names
names: list = list(df_acc.columns[1:len(df_acc.columns)])

pairs: list = list(itertools.combinations(names,2))

df_pairs = pd.DataFrame(columns=['iso1', 'iso2', 'genes'])

count: int = 0

for pair in pairs:
    count = count+1
    
    iso1: str = pair[0]
    iso2: str = pair[1]
    
    genes: str = str(df_acc[df_acc[iso1] != df_acc[iso2]].shape[0])
    
    new_genes: list = [iso1, iso2, genes]
    
    df_pairs.loc[count] = new_genes
    
df_pairs.to_csv('pair_genes.csv', index = False)
    
