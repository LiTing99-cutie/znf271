#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests, sys
import pandas as pd
import os
from optparse import OptionParser


# In[5]:


# Parse command line arguments
parser = OptionParser()
parser.add_option("-i", "--input", dest="input_file", help="gene_id to predict",default='/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict/to_predict_gene_id.txt')
parser.add_option("-o", "--output", dest="output_path", help="Output file path",default='/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict/cdna/')
parser.add_option("-a", "--archive", dest="archive_name", help="archive monthYear",default='Jul2022')
parser.add_option("-s", "--species", dest="species_name", help="species name or alias",default='homo_sapiens')
(options, args) = parser.parse_args()
print(options.archive_name)


# In[ ]:


list=pd.read_csv(options.input_file,header=None)
list=list.iloc[:,0].tolist()
list


# In[ ]:


server = "https://"+options.archive_name+".rest.ensembl.org"
ext = "/info/data/?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))


# In[ ]:


if not os.path.exists(options.output_path):
  os.makedirs(options.output_path)
for gene in list:
    prefix=os.path.splitext(gene)[0]
    ext = "/sequence/id/"+prefix+"?species="+options.species_name+";type=cdna;multiple_sequences=1"

    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()


    # print(r.text)
    r.text

    fo = open(options.output_path+"/"+gene+".fa", "w")
    fo.write(r.text)
    fo.close()

