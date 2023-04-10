
# packages
library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")
# rm previous data
# rm(list=ls())

# read data in
h <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed")
m <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.bed")
r <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.bed")
# ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/hrm.add_lnc.txt",header = F)
rn <- function(df){
  colnames(df) <- c("gene_name","len","cnt")
  return(df)
}
# change m and r gene_name with ortholog name in h
merge(rn(m),ortho,by.x = "gene_name",by.y = "V3") %>% mutate(gene_name=NULL,V2=NULL) %>% rename(gene_name=V1) %>% select(gene_name,len,cnt)-> m_rn
merge(rn(r),ortho,by.x = "gene_name",by.y = "V2") %>% mutate(gene_name=NULL,V3=NULL) %>% rename(gene_name=V1) %>% select(gene_name,len,cnt)-> r_rn
# output
fwrite(m_rn,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.rn.bed")
fwrite(r_rn,"/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.rn.bed")
fwrite(rn(h),"/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.rn.bed")
