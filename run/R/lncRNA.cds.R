

rm(list=ls())
library(tidyfst)
library(tibble)
library(dplyr)

map <- read.table("/home/user/data2/lit/project/ZNF271/02-APA/annotation/map/gene_id_vs_transcript_id.txt",col.names = c("gene_id","transcript_id","strand"))

lncRNA_cds_anno <- read.table("/home/user/data2/lit/project/ZNF271/02-APA/annotation/lncRNA_orf_prot_co_filter_genomic.txt",header = T)

merge(map,lncRNA_cds_anno,by="transcript_id") %>% filter(strand.x==strand.y) %>% merge(.,terminal_exon,"gene_id") %>% 
  rename(cds_length=length,cds_start=start,cds_end=end)%>% 
  mutate(transcript_length=rep(".",nrow(.))) %>% 
  .[,c("transcript_id","gene_id","transcript_length","cds_length","chr","strand",
       "proximal_start","proximal_end","distal_start","distal_end","cds_start","cds_end")] -> cds_lncRNA

saveRDS(cds_lncRNA,file = "annotation/cds_lncRNA.rds")
