# correction on library size

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

# read scale factor in
scale <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hm.sf.txt")
h$V3 <- h$V3/as.numeric(scale[1,1])
h %>% filter(V3>=2) -> h_s
m$V3 <- m$V3/as.numeric(scale[2,1])
m %>% filter(V3>=2) -> m_s
# ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/hrm.add_lnc.txt",header = F)
cal_len <- function(df){
  colnames(df) <- c("gene_name","len","cnt")
  df %>% group_by(gene_name) %>% summarise(max_len=max(len)) ->res
  return(res)
}

# change s_2 gene_name with ortholog name in s_1
merge(cal_len(m_s),ortho,by.x = "gene_name",by.y = "V3") %>% mutate(gene_name=NULL,V2=NULL) %>% rename(gene_name=V1) -> tmp_1
merge(cal_len(r),ortho,by.x = "gene_name",by.y = "V2") %>% mutate(gene_name=NULL,V3=NULL) %>% rename(gene_name=V1) -> tmp_2

# merge and calculate the length of longest exon in three species
merge(tmp_1,tmp_2,"gene_name") %>% merge(cal_len(h_s),"gene_name") %>% column_to_rownames("gene_name")  ->tmp
apply(tmp, 1, max) %>% data.frame() %>% rownames_to_column("gene_name") %>% rename(len=".")-> max_len

# fwrite
fwrite(max_len,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/max_len_s.txt",sep = '\t')

# rename 

rn <- function(df){
  colnames(df) <- c("gene_name","len","cnt")
  return(df)
}
# change m and r gene_name with ortholog name in h
merge(rn(m_s),ortho,by.x = "gene_name",by.y = "V3") %>% mutate(gene_name=NULL,V2=NULL) %>% rename(gene_name=V1) %>% select(gene_name,len,cnt)-> m_rn
# output
output_path_m="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m_s"
output_path_h="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h_s"
ifelse(!dir.exists(output_path_m),dir.create(output_path_m),FALSE)
ifelse(!dir.exists(output_path_h),dir.create(output_path_h),FALSE)
fwrite(m_rn,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m_s/m.cnt.rn.bed")
fwrite(rn(h_s),"/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h_s/h.cnt.rn.bed")
