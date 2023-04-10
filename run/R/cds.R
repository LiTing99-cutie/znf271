

#### Library packages ####
rm(list=ls())
library(tidyfst)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

#### Filter results ####
# 6813
res <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/final_res_ref_based_all_1_loose.with_geneName.txt")
# 6811
res <- na.omit(res)
# 4379
res[res$pau_a<=1 & res$pau_a>=-1 & res$pau_b>=-1 & res$pau_b<=1] -> res



#### Subset expression not change but pa change list ####
# annotate
res %>% mutate(diff_pa_usage=case_when(
  pvalue>0.05~"No change",
  pvalue<=0.05 & pau_a>pau_b~"Lengthening",
  pvalue<=0.05 & pau_a<pau_b~"Shortening"
)) %>% 
  mutate(diff_rpkm=case_when(
    p_value_1>0.05~"No change",
    p_value_1<=0.05 & rpkm_pro_a>rpkm_pro_b~"Up",
    p_value_1<=0.05 & rpkm_pro_a<rpkm_pro_b~"Down"
  )) -> res.mutate
# filter
filter(res.mutate,diff_pa_usage!="No change",diff_rpkm=="No change") -> expr_n_c_pa_c_lst
fwrite(expr_n_c_pa_c_lst,file = "output/final_list/expr_n_c_pa_c_lst.txt")

#### Merge expr_n_c_pa_c_lst with terminal exon annotation ####
terminal_exon=fread("annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt")
colnames(terminal_exon)=c("gene_name","chr","strand","proximal_start","proximal_end","distal_start","distal_end")
merge(expr_n_c_pa_c_lst[,c("gene_name","gene_id","gene_type")],terminal_exon,by="gene_name")->expr_n_c_pa_c_lst_m
filter(expr_n_c_pa_c_lst_m,gene_type!="protein_coding" & gene_type!="TEC") %>% select(gene_id) %>% 
  fwrite("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/orf_predict/loose_n_l/to_predict_gene_id.txt",col.names = F)

#### Prepare CDS annotation ####

cds_l=fread("annotation/cds_l_transcript_l_ensembl_107.txt",header = TRUE)
colnames(cds_l)=c("gene_id","transcript_id","transcript_length","cds_length")

cds_s_e =fread("annotation/gencode.v41.basic.annotation.cds.s_e.txt",header = FALSE)
colnames(cds_s_e)=c("transcript_id","strand","cds_start","cds_end","gene_name")

#### Identify cds type ####

cds_l %>% distinct(transcript_id,gene_id,cds_length,transcript_length) %>% na.omit() %>% 
  merge(.,expr_n_c_pa_c_lst_m,"gene_id") %>% merge(.,cds_s_e[,c("cds_start","cds_end","transcript_id")],"transcript_id") -> cds.na_omit.uniq.m
readRDS(file = "annotation/cds_lncRNA.rds") -> cds_lncRNA
rbind(cds_lncRNA,cds.na_omit.uniq.m) -> cds

# forward strand
cds %>% filter(strand=="+") %>% 
  filter(cds_end<=distal_end) %>% 
  group_dt(gene_id,mutate_dt(cds_end_max=max(cds_end),cds_end_min=min(cds_end),cds_length_max=max(cds_length))) %>% 
  filter(cds_length==cds_length_max) %>% 
  mutate(distal_type=case_when(
  distal_end>=cds_end~"coding",
  distal_end<cds_end_min~"non_coding",
  distal_end<cds_end & distal_end>=cds_end_min~"another protein product"),
  proximal_type=case_when(
    proximal_end>=cds_end~"coding",
    proximal_end<cds_end_min~"non_coding",
    proximal_end<cds_end & proximal_end>=cds_end_min ~ "another protein product"
)) -> tmp.1

# reverse strand
cds %>% filter(strand=="-") %>% 
  filter(cds_start>=distal_start) %>% 
  group_dt(gene_id,mutate_dt(cds_start_max=max(cds_start),cds_start_min=min(cds_start),cds_length_max=max(cds_length))) %>% 
  filter(cds_length==cds_length_max) %>% 
  mutate(distal_type=case_when(
    distal_start<=cds_start~"coding",
    distal_start>cds_start_max~"non_coding",
    distal_start>cds_start & distal_start<=cds_start_max~"another protein product"),
    proximal_type=case_when(
      proximal_start<=cds_start~"coding",
      proximal_start>cds_start_max~"non_coding",
      proximal_start>cds_start & proximal_start<=cds_start_max~"another protein product"
    )) -> tmp.2

distinct(rbind(tmp.1[,c("gene_id","proximal_type","distal_type")],tmp.2[,c("gene_id","proximal_type","distal_type")])) ->cds_type

merge(expr_n_c_pa_c_lst_m,cds_type,by="gene_id") -> expr_n_c_pa_c_lst_m_cds_type



#### lncRNA cds to verify in ribo_seq ####
cds_lncRNA %>% filter(strand=="+") %>% 
  filter(cds_end<=distal_end) %>% 
  group_dt(gene_id,filter_dt(cds_length==max(cds_length))) %>% 
  select(gene_id,chr,strand,cds_start,cds_end) %>% distinct_dt(gene_id,.keep_all = TRUE)->cds.1
cds_lncRNA %>% filter(strand=="-") %>% 
  filter(cds_start>=distal_start) %>% 
  group_dt(gene_id,filter_dt(cds_length==max(cds_length))) %>% 
  select(gene_id,chr,strand,cds_start,cds_end) %>% distinct_dt(gene_id,.keep_all = TRUE)->cds.2
rbind(cds.1,cds.2)->terminal_longest_cds
expr_n_c_pa_c_lst_m_cds_type %>% select(gene_id,gene_name,gene_type,proximal_type) %>% 
  merge(.,terminal_longest_cds,by="gene_id") -> coding_orf
write.table(coding_orf,file="output/final_list/lncRNA_cds.txt",quote=FALSE,row.names=FALSE)


#### Output for further analysis ####
expr_n_c_pa_c_lst_m_cds_type %>% filter(proximal_type=="another protein product" | proximal_type=="non_coding") ->typeIIAndIII
typeIIAndIII %>% fwrite(.,file="output/final_list/typeIIAndIII.txt")


