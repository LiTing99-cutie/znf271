#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests, sys
import pandas as pd
import os
from optparse import OptionParser


# In[ ]:


# Parse command line arguments
parser = OptionParser()
parser.add_option("-a", "--archive", dest="archive_name", help="archive monthYear",default='Jul2022')
parser.add_option("-s", "--species", dest="species_name", help="species name or alias",default='homo_sapiens')
(options, args) = parser.parse_args() 
print(options.archive_name)


# In[ ]:


cdna_co=pd.read_csv("lncRNA_orf_prot_co_filter.txt",sep='\t',header=None)
cdna_co


# In[ ]:


server = "https://"+options.archive_name+".rest.ensembl.org"
f=open('lncRNA_orf_prot_co_filter_genomic.txt','w')
print("transcript_id",'\t',"start",'\t',"end",'\t',"length",'\t',"strand",file=f)
for i in range(len(cdna_co.iloc[:,0])):
    ext = "/map/cdna/"+cdna_co.iloc[:,0][i]+"/"+str(cdna_co.iloc[:,2][i])+".."+str(cdna_co.iloc[:,3][i])+"?species="+options.species_name
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    
    list=[]
    for key in range(len(decoded['mappings'])):  
        list.append(decoded['mappings'][key]['start'])
        list.append(decoded['mappings'][key]['end'])
    print(cdna_co.iloc[:,1][i],'\t',min(list),'\t',max(list),cdna_co.iloc[:,4][i],'\t',cdna_co.iloc[:,5][i],file=f)
f.close()

