# arguments
args <- commandArgs(T)
# 1 s_1 file path
# 2 s_2 file path
# 3 ortholog(s_1_gene_name,s_2_gene_name) file path
# 4 output file path
# packages
library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")

# read test data in 
s_1 <- fread(args[1])
s_2 <- fread(args[2])

# rename columns and generate normalized length
rn_add_c <- function(df){
  colnames(df) <- c("gene_name","len","cnt")
  df %>% group_dt(gene_name,mutate_dt(per=len/max(len))) ->df
  return(df)
}

s_1 <- rn_add_c(s_1)
s_2 <- rn_add_c(s_2)

# change s_2 gene_name with ortholog name in s_1
# ortholog <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hm.txt",header = F)
ortholog <- fread(args[3],header = F)
merge(s_2,ortholog,by.x = "gene_name",by.y = "V2") %>% mutate(gene_name=NULL) %>% rename(gene_name=V1) ->s_2_1

# gene_lst
intersect(unique(s_1$gene_name),unique(s_2_1$gene_name)) -> gene_lst

# function
expand <- function(df,col){
  arr <- c()
  for (i in 1:nrow(df)){
    arr <- c(arr,rep(df[i,..col],times=df[i,"cnt"])) %>% unlist()->s
  }
  return(s)
}
GetPvalue <- function(gene,col){
  s_1 %>% filter(gene_name==gene) ->i_1
  s_2_1 %>% filter(gene_name==gene) ->i_2
  expand(i_1,col) -> x
  expand(i_2,col) -> y
  wilcox.test(x,y) -> res_1
  ks.test(x,y) -> res_2
  return(c(res_1$p.value,res_2$p.value))
}

# calculate
cal <- function(col){
  gene_lst %>% sapply(.,GetPvalue,col=col) %>% t() %>% data.frame() %>% rownames_to_column() ->res
  colnames(res) <- c("gene_name","p_value_wc","p_value_ks")
  return(res)
}

cal("per") ->res_1
cal("len") ->res_2

# significant gene number and intersect number
n_sig <- function(df){
  filter(df,p_value_wc<0.05) %>% nrow() ->n_1
  filter(df,p_value_ks<0.05) %>% nrow() ->n_2
  filter(df,p_value_ks<0.05 & p_value_wc<0.05) %>% nrow() ->n_3
  return(c(n_1,n_2,n_3))
}
print("gene_number")
length(gene_lst)
print("per")
n_sig(res_1)
print("len")
n_sig(res_2)
# intersect genes (absolute length and relative length)
filter(res_1,p_value_ks<0.05 & p_value_wc<0.05) %>% select(gene_name)->sig_1
filter(res_2,p_value_ks<0.05 & p_value_wc<0.05) %>% select(gene_name)->sig_2
intersect(sig_1,sig_2)->res
print("sig_gene_number")
nrow(res)

df <- data.frame(gene_name=gene_lst)
df$diff_1[df$gene_name %in% sig_1$gene_name] <- 'Diff'
df$diff_1[!df$gene_name %in% sig_1$gene_name] <- 'NotDiff'
df$diff_2[df$gene_name %in% sig_2$gene_name] <- 'Diff'
df$diff_2[!df$gene_name %in% sig_2$gene_name] <- 'NotDiff'
# fwrite
# fwrite(res,sep = '\t',file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hm.sig.h_g.txt")
fwrite(df,sep = '\t',file =args[4])
