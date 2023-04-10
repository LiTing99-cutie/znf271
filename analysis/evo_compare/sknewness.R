# arguments
args <- commandArgs(T)
# 1 file_path
# 2 output_path
# packages
library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")
# rm previous data
# rm(list=ls())

# read test data in 
test <- fread(args[1])
# test <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed")
colnames(test) <- c("gene_name","len","cnt")
test %>% group_dt(gene_name,mutate_dt(per=len/max(len))) ->test

# PAS number >2 gene
# for those PAS number=1; skewness will be NA because the denominator is 0
test %>% count(gene_name) %>% filter(n>1) %>% .$"gene_name" -> gene_lst

# function 
GetSknewness <- function(gene){
  test %>% filter(gene_name==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"len"],times=input[i,"cnt"])) %>% unlist()->s
  }
  return(skewness(s))
}

# calculate skewness for each gene
sapply(gene_lst,GetSknewness) %>% unlist() %>% data.frame() %>% rownames_to_column()-> res
colnames(res) <- c("gene_name","skewness")

# fwrite
fwrite(res,sep = '\t',file = args[2])

