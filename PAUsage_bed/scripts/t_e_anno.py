# Import libraries
import numpy as np
import argparse
import os
import warnings
warnings.filterwarnings('ignore')
import pandas as pd 
import pybedtools
from pybedtools import BedTool
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# Parse augments
parser = argparse.ArgumentParser()
parser.add_argument("--t_e",help="Unprocessed terminal exon")
parser.add_argument("--iso_anno",help="Iso_seq PA usage bed6+")
parser.add_argument("--out_path",help="Output path")
args=parser.parse_args()

# Read data in
te=pd.read_csv(args.t_e,sep='\t',header=None)
te.columns=["chr","strand","start","end","gene_id"]
pa_anno=pd.read_csv(args.iso_anno,sep='\t',header=None)
pa_anno.columns=["chr","start","end","gene_id","cnt","stand","pa_usage"]

# Generate terminal exon annotation
te_reverse=te[te.strand=='-']
te_forward=te[te.strand=='+']
tmp = [te_reverse, te_forward]
for i in range(len(tmp)):
    if i == 0:
        col = 'end' 
    else:
        col = 'start'
    # Extract terminal exons with most abundant 5' start site of isoforms
    tmp[i]["id_start"]=tmp[i]["gene_id"]+"_"+tmp[i][col].apply(str)
    id_start_c=tmp[i]["id_start"].value_counts().reset_index()
    id_start_c.columns=['id_start','c']
    id=id_start_c.id_start.str.split("_",expand=True).iloc[:,0]
    id_start_c['id']=id
    id_start_c_max=id_start_c.groupby(["id"],as_index=False).apply(lambda t:t[t.c==t.c.max()])
    # most abundant start as stable terminal exon
    AbFirstOrNot=id_start_c_max.id.value_counts().reset_index()
    AbFirstOrNot.columns=['id','c']
    AbFirstOrNot=AbFirstOrNot[AbFirstOrNot.c==1]['id']
    # merge
    fil=id_start_c_max.merge(AbFirstOrNot,on="id",how="inner").id_start
    tmp[i]=pd.merge(tmp[i],fil,on="id_start",how="inner").iloc[:,:5]  
    if i == 0:
        tmp[i] = tmp[i].groupby(['gene_id','chr','strand']).agg({"start":["max","max","min"], "end":"max"}).reset_index()
        tmp[i].columns = ['gene_id','chr','strand',"proximal_start","distal_end","distal_start",'proximal_end']
        tmp[i] = tmp[i][['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']]
    else:
        tmp[i] = tmp[i].groupby(['gene_id','chr','strand']).agg({"start":["min"], "end":["min","min","max"]}).reset_index()
        tmp[i].columns = ['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']
        tmp[i] = tmp[i][['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']]

# concat
tmp_concat=pd.concat([tmp[1],tmp[0]],axis=0)
tmp_concat
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
num=encode_anno._get_numeric_data()
num[num<0]=1
# intersect with PAS
encode_bed=BedTool.from_dataframe(encode_anno)
iso_bed=BedTool.from_dataframe(pa_anno)
len(encode_bed)
# intersect strand sensitive
iso_and_encode=iso_bed.intersect(encode_bed,wa=True,wb=True,s=True)
iso_and_encode_df=iso_and_encode.to_dataframe(header=None)
# some genes may overlap; make sure the same gene
iso_and_encode_df=iso_and_encode_df[iso_and_encode_df.iloc[:,3]==iso_and_encode_df.iloc[:,10]]
# subset columns
iso_and_encode_df=iso_and_encode_df.iloc[:,[0,1,2,3,5]]
iso_and_encode_df.columns=['chr','start','end','gene_id',"strand"]
len(iso_and_encode_df)

# make sure 2 PAS are presented in gene
cnt=pd.DataFrame(iso_and_encode_df.iloc[:,[3]].value_counts()).reset_index()       
m=pd.merge(iso_and_encode_df, cnt[cnt.iloc[:,1]>1], on='gene_id',how='right')
m_g=m.groupby(['gene_id']).agg({"end":["min","max"]}).reset_index()

# number of PAS site per gene and pas isoforms in te
cnt["new"]=(cnt.iloc[:,1]
.mask(cnt.iloc[:,1]>=4,"4+")
.mask(cnt.iloc[:,1]==2,"2")
.mask(cnt.iloc[:,1]==3,"3")
.mask(cnt.iloc[:,1]==1,"1"))
cnt.new.value_counts().to_csv(args.out_path+"/"+"pas_per_gene_fre_te.txt",header=0,sep='\t')
fo=open(args.out_path+"/"+"n_pa_isoforms_te.txt","w")
fo.write(str(len(iso_and_encode_df)))
fo.close()
# merge
pd.merge(m_g,tmp[1],on='gene_id',how='inner')
m_g_1=pd.merge(m_g,tmp[1],on='gene_id',how='inner').iloc[:,[0,4,5,6,2,2,3]]
m_g_1.columns = ['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']
m_g_1
pd.merge(m_g,tmp[0],on='gene_id',how='inner')
m_g_0=pd.merge(m_g,tmp[0],on='gene_id',how='inner').iloc[:,[0,4,5,3,7,2,3]]
m_g_0.columns = ['gene_id','chr','strand',"proximal_start","proximal_end","distal_start",'distal_end']
m_g_0

# proximal and distal lengh >100 bp
m_g_m=pd.concat([m_g_1,m_g_0])
m_g_m = m_g_m[(m_g_m.proximal_end-m_g_m.proximal_start >= 100) & (m_g_m.distal_end-m_g_m.distal_start >= 100)]

# write
m_g_m.to_csv(args.out_path+"/"+"terminal_exon_annotation.txt",sep="\t",header=0,index=0)

# how many isoforms left to be studied
tmp=pd.merge(iso_and_encode_df,m_g_m.iloc[:,0],on='gene_id',how='inner')
tmp.to_csv(args.out_path+"/"+"pa_isoforms_study.txt",sep="\t",header=0,index=0)