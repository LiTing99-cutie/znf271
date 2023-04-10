args <- commandArgs(T)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

# function
add_g_n <- function(df){
  df$gene_name <- str_split_fixed(df$V4,":",2)[,1]
  return(df)
}

case <- function(gene){
  filter(df_1,gene_name==gene) -> i1
  filter(df_2,gene_name==gene) -> i2
  filter(common,gene_name==gene) -> i3
  return(list(i1,i2,i3))
}


# human
df_1 <- fread("h/PAS.te.PAusage.scale.bed") %>% add_g_n(.)

# rhesus with orthologous name in h
df_2 <- fread("r/PAS.te.PAusage.scale.rn.bed") %>% add_g_n(.)

# human and rhesus intersect
unique(df_1$gene_name) %in% unique(df_2$gene_name) %>% table()
unique(df_2$gene_name) %in% unique(df_1$gene_name) %>% table()
intersect(df_1$gene_name,df_2$gene_name) %>% length()

# orthologous genes with all PASs in human reciprocol mapped (rhesus at least 1 PAS)
fread("h/PAS.te.cnt.scale.ToRheMac8.forCmp.bed6") %>% separate(V4,sep = ":",into = c("gene_name","pos")) %>% distinct(gene_name) -> gene_lst
intersect(gene_lst,distinct(df_2,gene_name)) -> gene_lst

# common
common <- fread("h/PAS.te.PAusage.scale.CommonPA.H.Rpos.tsv") %>% add_g_n(.)

## delta pau sum
common %>% group_by(gene_name) %>% summarise(delta_pau_sum=sum(abs(V11-V5))) -> delta_c

# species only
!df_1$V4 %in% common$V4 -> idx
df_1_o <- df_1[idx]
!df_2$V4 %in% common$V10 -> idx
df_2_o <- df_2[idx]

## calculate delta pau sum
df_1_o %>% group_by(gene_name) %>% summarise(delta_pau_sum=sum(V5)) -> delta_o_1
df_2_o %>% group_by(gene_name) %>% summarise(delta_pau_sum=sum(V5)) -> delta_o_2

# merge
merge(delta_o_1,delta_o_2,"gene_name",all = T) %>% merge(delta_c,"gene_name",all = T)%>% replace(is.na(.),0) %>% 
  rowwise() %>% mutate(delta_pau_sum_all=sum(c_across(2:4))) -> res

# only include genes in gene_lst
res$gene_name %in% gene_lst$gene_name -> idx
res[idx,] -> res

# inspection
filter(res,delta_pau_sum_all==1) %>% filter(delta_pau_sum.x==1 && delta_pau_sum.y==0 && delta_pau_sum==0) %>% nrow()
filter(res,delta_pau_sum_all==1) %>% nrow()
filter(res,delta_pau_sum_all==1) %>% filter(!(delta_pau_sum.x==1 && delta_pau_sum.y==0 && delta_pau_sum==0)) %>% nrow()

# output
fwrite(res,"output/hr.delta_pau_sum.txt",sep = '\t')

