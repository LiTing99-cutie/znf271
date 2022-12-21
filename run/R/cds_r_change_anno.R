# args <- commandArgs(TRUE)
# 1 PAusage.type_II.fil.bed6+

rm(list=ls())
suppressMessages(library(tidyfst))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

# PA site location bed 
PA_site <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/r_c_anno/PAusage.type_II_III.fil.bed6+")
colnames(PA_site) <- c("chr","start","end","gene_id","count","strand")

# cds length
cds=fread("/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/rhesus/cds_l_transcript_l_ensembl_97_r.txt",header = TRUE)
colnames(cds)=c("gene_id","transcript_id","transcript_length","cds_length")

cds %>% distinct(transcript_id,gene_id,cds_length,transcript_length) %>% na.omit() ->cds.uniq
cds.uniq$gene_id <- stringr::str_split_fixed(cds.uniq$gene_id,'\\.',2)[,1]
cds.uniq$transcript_id <- stringr::str_split_fixed(cds.uniq$transcript_id,'\\.',2)[,1]
# cds start and end for each transcript
info <- read.table("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/r_c_anno/t_id_chr_str_t_s_e_cds_s_e_gene_name.txt")
colnames(info) <- c("transcript_id","chr","strand","t_s","t_e","cds_s","cds_e","gene_id")

merge(distinct(select(PA_site,gene_id)),cds.uniq,by = "gene_id") %>% merge(.,info[,c(1:7)],by="transcript_id") ->tmp

merge(tmp[,c("gene_id","cds_length","cds_s","cds_e")],PA_site,"gene_id",allow.cartesian=TRUE) %>% 
  mutate(pas_id=paste0(gene_id,"_",start,"_",end))-> m2

m2 %>% filter(strand=="+") ->tmp.1
tmp.1[end<cds_e]$cds_length=0
tmp.1 %>% 
  group_dt(pas_id,filter_dt(cds_length==max(cds_length))) %>% 
  distinct_dt(pas_id,.keep_all = TRUE) %>% 
  group_dt(gene_id,mutate_dt(cds_l_max=max(cds_length))) %>% 
  filter(cds_l_max!=0) %>%
  mutate(cds_type=case_when(
    cds_length==cds_l_max~"primary",
    cds_length!=cds_l_max~"alternative",
  ))->tmp.1

m2 %>% filter(strand=="-") ->tmp.2
tmp.2[end>cds_s]$cds_length=0
tmp.2%>% 
  group_dt(pas_id,filter_dt(cds_length==max(cds_length))) %>% 
  distinct_dt(pas_id,.keep_all = TRUE) %>% 
  group_dt(gene_id,mutate_dt(cds_l_max=max(cds_length))) %>% 
  filter(cds_l_max!=0) %>%
  mutate(cds_type=case_when(
    cds_length==cds_l_max~"primary",
    cds_length!=cds_l_max~"alternative",
  ))->tmp.2

rbind(tmp.1,tmp.2) %>% dplyr::group_by(cds_type,gene_id) %>% dplyr::summarise(count_sum=sum(count))->tmp

fwrite(tmp,file ="/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/r_c_anno//cds_type_pa_count.txt" )

