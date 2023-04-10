#!/usr/bin/env python
# coding: utf-8

# # import packages

# In[40]:


import pandas as pd
import scipy.stats as ss
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity='all'
import scipy
scipy.__version__
import argparse

# # read rpkm calculated by stringtie

# In[41]:

parser = argparse.ArgumentParser()
parser.add_argument("--rpkm",help="RPKM file")
parser.add_argument("--output",help="Wilcox test result")
args=parser.parse_args()
rpkm=pd.read_csv(args.rpkm,sep='\t',header=None)


# # read metadata,fragmentation score,period

# In[42]:


metadata=pd.read_csv("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.txt",sep='\t')
frag_score=pd.read_csv("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt",sep='\t')
period=pd.read_csv("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period.txt",sep='\t')
metadata.head()
frag_score.head()
period.head()

# # rpkm_proximal and rpkm_distal

# In[44]:


rpkm.columns=["type","gene_id","rpkm","sample"]
rpkm_pro=rpkm.query("type.str.contains('proximal')")
rpkm_pro.columns=["type","gene_id","rpkm_pro","sample"]
rpkm_dis=rpkm.query("type.str.contains('distal')")
rpkm_dis.columns=["type","gene_id","rpkm_dis","sample"]
len(rpkm_pro)
len(rpkm_dis)


# # same order

# In[45]:


rpkm_pro_s=rpkm_pro.sort_values(by=["gene_id","sample"]).reset_index()
len(rpkm_pro_s)
rpkm_dis_s=rpkm_dis.sort_values(by=["gene_id","sample"]).reset_index()
len(rpkm_dis_s)


# # pau

# In[46]:


pau=pd.DataFrame(
    data={
        "pau" : rpkm_dis_s["rpkm_dis"] / rpkm_pro_s["rpkm_pro"],
        "sample" : rpkm_dis_s["sample"],
        "gene_id" : rpkm_dis_s["gene_id"]
    }
)
pau.head()


# # filter samples based on fragmentation score 

# In[47]:


tmp=pd.merge(metadata,frag_score,on="sample",how='right')
md_fil=tmp[tmp.frag_score>0.885]
md_fil_p=pd.merge(md_fil,period,on="Developmental_Stage")
md_fil_p.head()


# # merge and filter lowly expressed samples

# In[48]:


m1=pd.merge(rpkm_dis_s,rpkm_pro_s,on=["gene_id","sample"])
m2=pd.merge(m1,pau,on=["gene_id","sample"])
final=pd.merge(m2,md_fil_p,on="sample",how='right')
final=final[final.rpkm_pro>1]
len(final)
final.head()


# # wilcox test

# In[49]:


sta_res=[]
array=rpkm.gene_id.drop_duplicates().values
for i in array:
    pau_b=final[(final.gene_id==i) & (final["period.abbre.abbre"]=="Before birth")].pau
    pau_a=final[(final.gene_id==i) & (final["period.abbre.abbre"]=="After birth")].pau
    rpkm_pro_b=final[(final.gene_id==i) & (final["period.abbre.abbre"]=="Before birth")].rpkm_pro
    rpkm_pro_a=final[(final.gene_id==i) & (final["period.abbre.abbre"]=="After birth")].rpkm_pro    
    if len(pau_a)>=3 and len(pau_b)>=3:
        sta_res.append([i,ss.mannwhitneyu(pau_b,pau_a,use_continuity=True).pvalue,pau_a.median(),pau_b.median(),
        ss.mannwhitneyu(rpkm_pro_b,rpkm_pro_a,use_continuity=True).pvalue,rpkm_pro_a.mean(),rpkm_pro_b.mean()])
# 2022-11-23
# sta_res=pd.DataFrame(sta_res,columns=["gene_id","pvalue","pau_a","pau_b","p_value_1","rpkm_pro_a","rpkm_pro_b"])
sta_res=pd.DataFrame(sta_res,columns=["gene_name","pvalue","pau_a","pau_b","p_value_1","rpkm_pro_a","rpkm_pro_b"])
sta_res.head()


# # write to file

# In[50]:


sta_res.to_csv(args.output,sep='\t',header=True,index=0)

