

rm(list=ls())
library(tidyfst)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

# PA site location bed 
PA_site <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/PAusage.type_II.fil.bed6+")
colnames(PA_site) <- c("chr","start","end","gene_id","count","strand")

# cds length
cds=fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/cds_l_transcript_l_ensembl_107.txt",header = TRUE)
colnames(cds)=c("gene_id","transcript_id","transcript_length","cds_length")

cds %>% distinct(transcript_id,gene_id,cds_length,transcript_length) %>% na.omit() ->cds.uniq

# gene_id and gene_symbol map
map <- read.table("/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/map/ensembl_gene_id_type_symbol.txt")
colnames(map) <- c("gene_id","type","gene_name")

# cds start and end for each transcript
info <- read.table("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/t_id_chr_str_t_s_e_cds_s_e_gene_name.txt")
colnames(info) <- c("transcript_id","chr","strand","t_s","t_e","cds_s","cds_e","gene_name")

merge(distinct(select(PA_site,gene_id)),cds.uniq,"gene_id") %>% merge(info,"transcript_id") ->  tmp

merge(map[,c(1,3)],PA_site,"gene_id") %>% select(-gene_id)-> m1

merge(tmp[,c("gene_name","cds_length","cds_s","cds_e")],m1,"gene_name",allow.cartesian=TRUE) %>% 
  mutate(pas_id=paste0(gene_name,"_",start,"_",end))-> m2

m2 %>% filter(strand=="+") %>% filter(end>=cds_e) %>% 
  group_dt(pas_id,filter_dt(cds_length==max(cds_length))) %>% 
  distinct_dt(pas_id,.keep_all = TRUE) %>% 
  group_dt(gene_name,mutate_dt(cds_l_max=max(cds_length))) %>% 
  mutate(cds_type=case_when(
    cds_length==cds_l_max~"primary",
    cds_length!=cds_l_max~"alternative",
  ))->tmp.1

m2 %>% filter(strand=="-") %>% filter(end<=cds_s) %>% 
  group_dt(pas_id,filter_dt(cds_length==max(cds_length))) %>% 
  distinct_dt(pas_id,.keep_all = TRUE) %>% 
  group_dt(gene_name,mutate_dt(cds_l_max=max(cds_length))) %>% 
  mutate(cds_type=case_when(
    cds_length==cds_l_max~"primary",
    cds_length!=cds_l_max~"alternative",
  ))->tmp.2

rbind(tmp.1,tmp.2) %>% group_by(cds_type,gene_name) %>% summarise(count_sum=sum(count))->tmp
# 
# write.table(unique(tmp$gene_name),file="/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/cds_type_II.gene.gene_name.lst",
#             col.names = FALSE,row.names = FALSE,quote = FALSE)






