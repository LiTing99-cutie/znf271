# correction on library size

# arguments
args <- commandArgs(T)
# 1 file_path
# 2 output_path
# packages
library("e1071")
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

# read data in 
h <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h_s/h.cnt.rn.bed")
r <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.rn.bed")
m <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m_s/m.cnt.rn.bed")
max_len <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/max_len_s.txt")

# gene_lst
intersect(unique(h$gene_name),max_len$gene_name) %>% intersect(.,unique(m$gene_name)) %>% intersect(.,unique(r$gene_name)) ->gene_lst
# rounding
h$cnt <- round(h$cnt)
m$cnt <- round(m$cnt)
# function 
GetSknewness <- function(df,gene){
  df %>% filter(gene_name==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"len"],times=input[i,"cnt"])) %>% unlist()->s
  }
  filter(max_len,gene_name==gene) %>% .$"len" -> gene_max_len
  s_add <- c(1,gene_max_len,s)
  return(e1071::skewness(s_add))
}

# calculate skewness for each gene
final <- function(df){
  sapply(gene_lst,GetSknewness,df=df) %>% unlist() %>% data.frame() %>% rownames_to_column()-> res
  colnames(res) <- c("gene_name","skewness")
  return(res)
}

# apply to each species
s <- list(h,r,m)
lapply(s,final) ->res

# merge
merge(res[[1]],res[[2]],by="gene_name") %>% merge(res[[3]],by="gene_name") -> res_df

# fwrite
fwrite(res_df,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw_s.txt")
