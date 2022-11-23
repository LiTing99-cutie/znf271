# %%%
import pandas as pd
import argparse
import os

# %%
parser = argparse.ArgumentParser()
parser.add_argument("--vs",help="ensembl id and gene symbol map")
parser.add_argument("--file",help="File to add gene symbl")
args=parser.parse_args()

# %% 
vs=pd.read_csv(args.vs,header=None,sep='\t')
vs.columns=["gene_id","gene_type","gene_name"]
vs

# %%
file=pd.read_csv(args.file,sep='\t')
file
# %%
m=pd.merge(vs,file,on="gene_name",how="inner")
# %%
m.to_csv(os.path.splitext(args.file)[0]+".with_geneName.txt",sep="\t",index=0)
