
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
PN <- function(cutoff){
  res_df$diff <- if_else((res_df$skewness.x<=cutoff & res_df$skewness.x>=-cutoff & res_df$skewness.y<=cutoff & res_df$skewness.y>=-cutoff & 
                            res_df$skewness<=cutoff & res_df$skewness>=-cutoff) | (res_df$skewness.x <= -cutoff & res_df$skewness.y <= -cutoff & 
                                                                         res_df$skewness <= -cutoff) | (res_df$skewness.x >= cutoff & res_df$skewness.y >= cutoff & res_df$skewness >= cutoff),"NotDiff","Diff")
  return(res_df)
}

evo <- list(PN(0.1),PN(0.2),PN(1))

# fwrite
saveRDS(evo,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_2_evo.rds")



