#### Library packages ####
rm(list=ls())
suppressMessages(library(tidyfst))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

#### Read files in ####

map <- read.table("annotation/map/gene_id_vs_transcript_id_strand.txt",col.names = c("gene_id","transcript_id","strand"))
lncRNA_cds_anno <- read.table("analysis/orf_predict/loose_n_l/lncRNA_orf_prot_co_filter_genomic.txt",header = T)
expr_n_c_pa_c_lst <- fread("output/final_list/expr_n_c_pa_c_lst.txt")
terminal_exon=fread("annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt")
colnames(terminal_exon)=c("gene_name","chr","strand","proximal_start","proximal_end","distal_start","distal_end")

#### Output ####
merge(expr_n_c_pa_c_lst[,c("gene_name","gene_id","gene_type")],terminal_exon,by="gene_name")->expr_n_c_pa_c_lst_m
merge(map,lncRNA_cds_anno,by="transcript_id") %>% filter(strand.x==strand.y) %>% merge(.,expr_n_c_pa_c_lst_m,"gene_id") %>% 
  rename(cds_length=length,cds_start=start,cds_end=end)%>% 
  mutate(transcript_length=rep(".",nrow(.))) %>% 
  .[,c("transcript_id","gene_id","transcript_length","cds_length","gene_name","gene_type","chr","strand",
       "proximal_start","proximal_end","distal_start","distal_end","cds_start","cds_end")] -> cds_lncRNA

saveRDS(cds_lncRNA,file = "annotation/cds_lncRNA.rds")
fwrite(cds_lncRNA,file = "annotation/cds_lncRNA.txt",sep = '\t')
