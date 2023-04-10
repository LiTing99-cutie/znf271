
# packages
library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")

# read data in 
res_df <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw.txt")

# select genes
res_df %>% column_to_rownames("gene_name") ->res_df_1
res_df_1>0 ->tmp
tmp <- as.data.frame(tmp)
# human-specific positive skewness
filter(tmp,skewness.x & !skewness.y & !skewness) %>% rownames() -> h_s_p_skew_g
# human-specific negative skewness
filter(tmp,!skewness.x & skewness.y & skewness) %>% rownames() -> h_s_n_skew_g
# all
h_s <- c(h_s_p_skew_g,h_s_n_skew_g)
# number
length(h_s_p_skew_g)+length(h_s_n_skew_g)
apply(tmp, 1, sum) %>% data.frame() %>% rename(positive_skewness_n=".") -> tmp_1

tmp_1 %>% rownames_to_column("gene_name") %>% mutate(diff=case_when(positive_skewness_n==1 | positive_skewness_n==2 ~ "Diff",
                                                                    positive_skewness_n==0 | positive_skewness_n==3 ~ "NotDiff")) -> evo

evo$h_s_or_not <- if_else(evo$gene_name %in% h_s,"h_s","n_h_s")
evo$h_s_p_or_not <- if_else(evo$gene_name %in% h_s_p_skew_g,"h_s","n_h_s")
evo$h_s_n_or_not <- if_else(evo$gene_name %in% h_s_n_skew_g,"h_s","n_h_s")
# fwrite
fwrite(evo,sep = '\t',file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo.txt")

# take cases as example
filter(tmp_1,positive_skewness_n==0) %>% head(1)
filter(tmp_1,positive_skewness_n==1) %>% head(1)
filter(tmp_1,positive_skewness_n==2) %>% head(1)
filter(tmp_1,positive_skewness_n==3) %>% head(1)
# human-specific positive skewness
filter(tmp,skewness.x & !skewness.y & !skewness) %>% head(1)
# human-specific negative skewness
filter(tmp,!skewness.x & skewness.y & skewness) %>% head(1)

GetTab <- function(df,gene){
  df %>% filter(gene_name==gene) ->input
  return(input)
}

GetHist <- function(df,gene){
  df %>% filter(gene_name==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"len"],times=input[i,"cnt"])) %>% unlist()->Terminal_exon_length
  }
  hist(Terminal_exon_length)
}

GetSknewness_old <- function(df,gene){
  df %>% filter(gene_name==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"len"],times=input[i,"cnt"])) %>% unlist()->Terminal_exon_length
  }
  skewness(Terminal_exon_length)
}

case <- function(df,gene){
  GetTab(df,gene) -> res_1
  # GetHist(df,gene) ->res_2
  GetSknewness_old(df,gene) ->res_2
  GetSknewness(df,gene) ->res_3
  l <- list(res_1,res_2,res_3)
  return(l)
}

GetSknewness <- function(df,gene){
  df %>% filter(gene_name==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"len"],times=input[i,"cnt"])) %>% unlist()->s
  }
  filter(max_len,gene_name==gene) %>% .$"len" -> gene_max_len
  s_add <- c(1,gene_max_len,s)
  return(skewness(s_add))
}
genes <- c("A4GALT","AACS","AAMDC","AAGAB","ACBD7","ZNF271P")
lapply(s,case,genes[6])


