# Use a different differential res

rm(list=ls())
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")
library("clusterProfiler")
library("DOSE")
library("stringr")

res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

res %>% filter(p_pau<0.05) %>% select(gene_name) -> l_1
res %>% filter(abs(delta_pau)>=0.2) %>% select(gene_name) -> l_2

table(l_2$gene_name %in% l_1$gene_name) -> inter_1_2

res %>% filter(p_expr<0.05) %>% select(gene_name) -> l_3
res %>% filter(abs(fc)>=1) %>% select(gene_name) -> l_4
table(l_4$gene_name %in% l_3$gene_name) -> inter_3_4

pau <- data.frame(s_1_n=nrow(l_1),s_2_n=nrow(l_2),inter_1_2=inter_1_2[2]/(inter_1_2[1]+inter_1_2[2]))
expr <- data.frame(s_1_n=nrow(l_3),s_2_n=nrow(l_4),inter_1_2=inter_3_4[2]/(inter_3_4[1]+inter_3_4[2]))

GO <- function(df){
  enrich <- enrichGO(gene=df,
                     OrgDb="org.Hs.eg.db",
                     ont="ALL",
                     keyType = "SYMBOL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  dotplot(enrich)
}
GO(l_1$gene_name)
GO(l_2$gene_name)
GO(l_3$gene_name)
GO(l_4$gene_name)

l_1[!l_1$gene_name %in% l_3$gene_name,] ->l_1_3
l_2[!l_2$gene_name %in% l_4$gene_name,] ->l_2_4

GO(l_1_3)
GO(l_2_4)

table(l_1_3 %in% l_2_4)
table(l_2_4 %in% l_1_3)


res %>% filter(p_pau<0.001) %>% select(gene_name) -> l_5
res %>% filter(p_expr<0.001) %>% select(gene_name) -> l_6

table(l_2$gene_name %in% l_5$gene_name)
table(l_4$gene_name %in% l_6$gene_name)
l_5[!l_5$gene_name %in% l_6$gene_name,] ->l_5_6
GO(l_5_6)
