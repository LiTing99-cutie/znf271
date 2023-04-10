# arguments
args <- commandArgs(T)
# 1 file_path
# 2 output_path
# packages
library("moments")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")
# rm previous data
# rm(list=ls())


# read data in to get longest exon
h <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed")
m <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.bed")
r <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.bed")
# ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/hrm.add_lnc.txt",header = F)
cal_len <- function(df){
  colnames(df) <- c("gene_name","len","cnt")
  df %>% group_by(gene_name) %>% summarise(max_len=max(len)) ->res
  return(res)
}

# change s_2 gene_name with ortholog name in s_1
merge(cal_len(m),ortho,by.x = "gene_name",by.y = "V3") %>% mutate(gene_name=NULL,V2=NULL) %>% rename(gene_name=V1) -> tmp_1
merge(cal_len(r),ortho,by.x = "gene_name",by.y = "V2") %>% mutate(gene_name=NULL,V3=NULL) %>% rename(gene_name=V1) -> tmp_2

# merge and calculate the length of longest exon in three species
merge(tmp_1,tmp_2,"gene_name") %>% merge(cal_len(h),"gene_name") %>% column_to_rownames("gene_name")  ->tmp
apply(tmp, 1, max) %>% data.frame() %>% rownames_to_column("gene_name") %>% rename(len=".")-> max_len

# fwrite
fwrite(max_len,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/max_len.txt",sep = '\t')
