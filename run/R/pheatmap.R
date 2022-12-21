
rm(list=ls())
library(tidyfst)
library(tibble)
library(pheatmap)
library(dplyr)
library(magrittr)

setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

res=fread("output/final_list/final_res_ref_based_all_1.with_geneName.txt")

pdf("output/R/pa_usage_dist_a.pdf")
hist(res$pau_a,xlab = "PA usage after birth",ylab="Count",main=NULL,col = "#2b83ba")
dev.off()
pdf("output/R/pa_usage_dist_b.pdf")
hist(res$pau_b,xlab = "PA usage before birth",ylab="Count",main=NULL,col = "#2b83ba")
dev.off()
res=filter(res,pvalue<0.001)
res=filter(res,pau_a <2 &  pau_b<2)
res$diff=res$pau_a-res$pau_b
res.o=res[order(res$diff),]
res.o$gene_type[grep("pseudogene",res.o$gene_type)]="Pseudogene"
res.o$gene_type[grep("protein_coding",res.o$gene_type)]="Protein_coding"
rownames(res.o)=NULL
res.o %<>% distinct_dt(gene_name,.keep_all = T)
ph.input=res.o[,c("gene_name","pau_b","pau_a")] %>% column_to_rownames("gene_name") %>% t() %>% data.frame()


annotation_col = data.frame(
  Diff = if_else(res.o$diff>0,"Lengthening","Shortening"),
  GeneType=res.o$gene_type,
  row.names = colnames(ph.input))

ann_colors = list(
  Diff = c(Shortening= '#fdae61',Lengthening='#abdda4'),
  GeneType = c(lncRNA='#7fc97f',Protein_coding='#beaed4',Pseudogene='#fdc086'))

col_number=match("ZNF271P",colnames(ph.input))
labels_col = rep("", ncol(ph.input))
labels_col[c(col_number)]=c("ZNF271P")
labels_row=c("Before birth","After birth")


ph.input %>% pheatmap(cluster_rows = FALSE,color = colorRampPalette(c("#ffffd1", "#751428"))(50),labels_col = labels_col,labels_row=labels_row,
                      cluster_cols = FALSE,  annotation_col = annotation_col,annotation_colors = ann_colors) ->p

pdf("output/R/pheatmap.pdf")
p
dev.off()
