# %%
import numpy as np
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

# %%
import pandas as pd 
import pybedtools
from pybedtools import BedTool
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# %%
parser = argparse.ArgumentParser()
parser.add_argument("--t_e",help="Unprocessed terminal exon")
parser.add_argument("--iso_anno",help="Iso_seq PA usage bed6+")
parser.add_argument("--t_e_out",help="Processed terminal exon file")
args=parser.parse_args()
# currentpath=os.path.abspath(os.path.dirname(__file__))
path="/home/user/data2/lit/project/ZNF271/02-APA/"

# %%
data=pd.read_csv(path+"/"+args.t_e,sep='\t',header=None)
data.columns=["chr","strand","start","end","gene_id"]
data.head()

# %%
pa_anno=pd.read_csv(path+"/"+args.iso_anno,sep='\t',header=None)
pa_anno.columns=["chr","start","end","gene_id","cnt","stand","pa_usage"]
pa_anno=pa_anno[pa_anno.cnt>1]

# %%
data_reverse=data[data.strand=='-']
data_forward=data[data.strand=='+']
data_reverse.head()

# %%
tmp = [data_reverse, data_forward]
for i in range(len(tmp)):
    if i == 0:
        col = 'end' 
    else:
        col = 'start'
        
    cnt = pd.DataFrame(tmp[i].gene_id.value_counts()).reset_index()
    cnt.columns = ['gene_id','cnt1']
    tmp[i] = pd.merge(tmp[i], cnt, on='gene_id',how='left')

    cnt = pd.DataFrame(tmp[i][['gene_id', col]].value_counts()).reset_index()
    cnt.columns = ['gene_id',col,'cnt2']
    tmp[i] = pd.merge(tmp[i], cnt, on=['gene_id',col],how='left')
    
    tmp[i]['freq'] = tmp[i].cnt2 / tmp[i].cnt1
    tmp[i] = tmp[i][(tmp[i]['freq'] > 0.5) & (tmp[i]['cnt2'] > 1)]
    
    if i == 0:
        tmp[i] = tmp[i].groupby(['gene_id','chr','strand']).agg({"start":["max","max","min"], "end":"max"}).reset_index()
        tmp[i].columns = ['gene_id','chr','strand',"proximal_start","distal_end","distal_start",'proximal_end']
        tmp[i] = tmp[i][['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']]
        # tmp[i] = tmp[i][(tmp[i].proximal_end-tmp[i].proximal_start >= 100) & (tmp[i].distal_end-tmp[i].distal_start >= 100)]
        
    else:
        tmp[i] = tmp[i].groupby(['gene_id','chr','strand']).agg({"start":["min"], "end":["min","min","max"]}).reset_index()
        tmp[i].columns = ['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']
        tmp[i] = tmp[i][['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']]
        # tmp[i] = tmp[i][(tmp[i].proximal_end-tmp[i].proximal_start >= 100) & (tmp[i].distal_end-tmp[i].distal_start >= 100)]
    
    
        

# %%
tmp_concat=pd.concat([tmp[1],tmp[0]],axis=0)
tmp_concat

# %%
encode_anno_1=pd.DataFrame(data={
    'chr':tmp[1]['chr'],
    'start':tmp[1]['proximal_start'],
    'end':tmp[1]['distal_end']+30,
    'name':tmp[1]['gene_id'],
    'score':".",
    'strand':tmp[1]['strand']
})

encode_anno_0=pd.DataFrame(data={
    'chr':tmp[0]['chr'],
    'start':tmp[0]['distal_start']-30,
    'end':tmp[0]['proximal_end'],
    'name':tmp[0]['gene_id'],
    'score':".",
    'strand':tmp[0]['strand']
})


encode_anno=pd.concat([encode_anno_0,encode_anno_1],axis=0)

# %%
encode_bed=BedTool.from_dataframe(encode_anno)
iso_bed=BedTool.from_dataframe(pa_anno)
len(encode_bed)
# intersect
iso_and_encode=iso_bed.intersect(encode_bed,wa=True,wb=True)
iso_and_encode_df=iso_and_encode.to_dataframe(header=None)
# some genes may overlap
iso_and_encode_df=iso_and_encode_df[iso_and_encode_df.iloc[:,3]==iso_and_encode_df.iloc[:,10]]

iso_and_encode_df.to_csv(path+"/"+args.t_e_out,sep="\t",header=0,index=0)