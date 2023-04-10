

library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")

# read data in 
res_df <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw.txt")

# select genes base in delta skewness
res_df %>% column_to_rownames("gene_name") ->res_df_1
range <- function(x){
  max(x)-min(x)
}
apply(res_df_1, 1, range) %>% data.frame() %>% rename(range=".") %>% arrange(range) %>% summary(.$range)
apply(res_df_1, 1, range) %>% data.frame() %>% rownames_to_column("gene_name")  %>% rename(range=".") %>% arrange(range) -> tmp
tmp %>% tail(nrow(res_df_1)/10) -> tmp_1
tmp$diff_based_on_delta_skewness <- if_else(tmp$gene_name %in% tmp_1$gene_name,"Diff","NotDiff")

# fwrite
fwrite(tmp,sep = '\t',file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_1_evo.txt")
