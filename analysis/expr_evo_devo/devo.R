suppressMessages(library(tidyfst,quietly = T))
suppressMessages(library(tibble,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(ggpubr,quietly = T))
suppressMessages(library(ggrepel,quietly = T))
suppressMessages(library(data.table,quietly = T))
suppressMessages(library(parallel,quietly = T))
#### Read rpkm ####
rpkm_path <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/expr_evo_devo/rpkm/","human","/stringtie.rpkm.txt")
h <- fread(rpkm_path)
colnames(h) <- c("gene_name","rpkm","sample")
#### Read metadata ####
metadata=fread("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.txt",sep='\t')
period=fread("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period.txt",sep='\t')

merge(metadata,period,"Developmental_Stage") %>% merge(h,"sample") ->res

# unique(res$gene_name) -> gene_lst
evo <- readRDS("output/anova.rds")
evo$gene_name -> gene_lst 

GetPvalue <- function(res,gene){
  library(dplyr)
  library(magrittr)
  res %>% filter(gene_name==gene) ->input
  input_b <- filter(input,period.abbre.abbre=="Before birth")
  input_a <- filter(input,period.abbre.abbre=="After birth")
  p <- wilcox.test(input_b$rpkm,input_a$rpkm)$p.value
  fc <- log2(mean(input_a$rpkm)/mean(input_b$rpkm))
  out <- list("p_expr"=p,"expr_embryo"=mean(input_b$rpkm),"expr_postnatal"=mean(input_a$rpkm))
  return(out)
}

sapply(gene_lst,GetPvalue,res=res) %>% apply(MARGIN = 2,unlist) %>% t() %>% data.frame()%>%
  rownames_to_column(.,"gene_name") -> res_f

fwrite(res_f,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/expr_evo_devo/output/devo.txt",sep = '\t')
