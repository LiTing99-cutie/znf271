

# library packages
rm(list=ls())
library(tidyfst)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("/home/user/data2/lit/project/ZNF271/02-APA")


# read files in 
final_list=fread("/home/user/data2/lit/project/ZNF271/02-APA/output/final_list/diff.lst",header = TRUE)
terminal_exon=fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1.txt")
colnames(terminal_exon)=c("gene_id","chr","strand","proximal_start","proximal_end","distal_start","distal_end")
cds=fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/cds_l_transcript_l_ensembl_107.txt",header = TRUE)
colnames(cds)=c("gene_id","transcript_id","transcript_length","cds_length")
cds_s_e =fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/gencode.v41.annotation.cds.s_e.txt",header = FALSE)
colnames(cds_s_e)=c("transcript_id","strand","cds_start","cds_end","gene_id")
# map=fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/map/ensembl_gene_id_type_symbol.txt",header = F)
# colnames(map) <- c("gene_id","gene_type","gene_symbol")

# identify orf type

## cds annotation(genomic start;genomic end;length)
cds %>% distinct(transcript_id,gene_id,cds_length,transcript_length) %>% na.omit() %>% 
  merge(.,terminal_exon,"gene_id") %>% merge(.,cds_s_e[,c("cds_start","cds_end","transcript_id")],"transcript_id") -> cds.na_omit.uniq.m
readRDS(file = "annotation/cds_lncRNA.rds") -> cds_lncRNA
rbind(cds_lncRNA,cds.na_omit.uniq.m) -> cds

cds %>% filter(strand=="+") %>% 
  filter(cds_end<=distal_end) %>% 
  group_dt(gene_id,mutate_dt(cds_end_max=max(cds_end),cds_end_min=min(cds_end),cds_length_max=max(cds_length))) %>% 
  filter(cds_length==cds_length_max) %>% 
  merge(final_list,by="gene_id") %>% 
  mutate(distal_type=case_when(
  distal_end>=cds_end~"coding",
  distal_end<cds_end_min~"non_coding",
  distal_end<cds_end & distal_end>=cds_end_min~"another protein product"),
  proximal_type=case_when(
    proximal_end>=cds_end~"coding",
    proximal_end<cds_end_min~"non_coding",
    proximal_end<cds_end & proximal_end>=cds_end_min ~ "another protein product"
)) -> tmp.1

cds %>% filter(strand=="-") %>% 
  filter(cds_start>=distal_start) %>% 
  group_dt(gene_id,mutate_dt(cds_start_max=max(cds_start),cds_start_min=min(cds_start),cds_length_max=max(cds_length))) %>% 
  filter(cds_length==cds_length_max) %>% 
  merge(final_list,by="gene_id") %>% 
  mutate(distal_type=case_when(
    distal_start<=cds_start~"coding",
    distal_start>cds_start_max~"non_coding",
    distal_start>cds_start & distal_start<=cds_start_max~"another protein product"),
    proximal_type=case_when(
      proximal_start<=cds_start~"coding",
      proximal_start>cds_start_max~"non_coding",
      proximal_start>cds_start & proximal_start<=cds_start_max~"another protein product"
    )) -> tmp.2

  
distinct(rbind(tmp.1[,c("gene_id","proximal_type","distal_type")],tmp.2[,c("gene_id","proximal_type","distal_type")])) ->final_type

merge(final_list,final_type,by="gene_id") -> final


write.table(file="/home/user/data2/lit/project/ZNF271/02-APA/output/final_list/lncRNA_cds.txt",quote=FALSE,row.names=FALSE)

# alternative protein product 
final %>% select(gene_id,proximal_type,gene_type) %>% filter(proximal_type=="another protein product") %>% 
  merge(select(cds,gene_id,proximal_start,proximal_end,cds_start,cds_end,strand),.,by="gene_id") ->tmp
filter(tmp,strand=="+") %>% filter(proximal_end>=cds_end) -> tmp.1
filter(tmp,strand=="-") %>% filter(proximal_start<=cds_start) -> tmp.2
rbind(tmp.1,tmp.2)->proximal_alternative_orf

# terminal longset cds
cds %>% filter(strand=="+") %>% 
  filter(cds_end<=distal_end) %>% 
  group_dt(gene_id,filter_dt(cds_length==max(cds_length))) %>% 
  select(gene_id,chr,strand,cds_start,cds_end) %>% distinct_dt(gene_id,.keep_all = TRUE)->cds.1
cds %>% filter(strand=="-") %>% 
  filter(cds_start>=distal_start) %>% 
  group_dt(gene_id,filter_dt(cds_length==max(cds_length))) %>% 
  select(gene_id,chr,strand,cds_start,cds_end) %>% distinct_dt(gene_id,.keep_all = TRUE)->cds.2
rbind(cds.1,cds.2)->terminal_longest_cds

final %>% select(gene_id,gene_type) %>% merge(.,terminal_longest_cds,by="gene_id") -> coding_orf

# write 
filter(final,gene_type=="lncRNA" | gene_type=="transcribed_unitary_pseudogene") %>% 
  write.table(file="/home/user/data2/lit/project/ZNF271/02-APA/output/final_list/lncRNA_pro_non_coding.txt",quote=FALSE,row.names=FALSE)
filter(coding_orf,gene_type=="lncRNA" | gene_type=="transcribed_unitary_pseudogene") %>% 
  write.table(file="/home/user/data2/lit/project/ZNF271/02-APA/output/final_list/lncRNA_cds.txt",quote=FALSE,row.names=FALSE)
write.table(final_type,file="/home/user/data2/lit/project/ZNF271/02-APA/output/final_list/cds_type.txt",quote=FALSE,row.names=FALSE)
# plot
p <- as.data.frame(table(final$proximal_type)) %>% mutate(type=c("Another orf","Longest orf","Non coding")) %>% 
  ggpie("Freq",
        fill = "type", color = "white",
        palette = c("#00AFBB", "#E7B800", "#FC4E07"),
        lab.pos = "out", lab.font = "white")+
theme(legend.title = element_blank(),plot.title = element_text(hjust=0.5,vjust = -80))
  
  ggsave("output/R/cds_effect.pdf",p,width = 5,height = 5)
