args <- commandArgs(TRUE)
# 1 species
# 2 gene
#### Library packages ####
suppressMessages(library(tidyfst,quietly = T))
suppressMessages(library(tibble,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(ggpubr,quietly = T))
suppressMessages(library(ggrepel,quietly = T))
suppressMessages(library(data.table,quietly = T))

#### Read rpkm ####
rpkm_path <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare/",args[1],"/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt")
rpkm <- fread(rpkm_path)
colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]

#### Merge ####
metadata_path <- paste0("/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/",args[1],".metadata.clean.txt")
metadata <- fread(file=metadata_path)

merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(metadata,"sample") %>%
  select(gene_id,rpkm.x,rpkm.y,stage) %>% filter(rpkm.x>1) %>% mutate(pau=rpkm.y/rpkm.x)-> res

#### Prepare gene list ####
unique(res$gene_id) -> gene_lst

GetPvalue <- function(gene){
  res %>% filter(gene_id==gene) ->input
  input_b <- filter(input,stage=="embryo")
  input_a <- filter(input,stage=="postnatal")
  if(nrow(input_a)>=8 && nrow(input_b)>=8){
    p_1 <- wilcox.test(input_b$rpkm.x,input_a$rpkm.x)$p.value
    p_2 <- wilcox.test(input_b$pau,input_a$pau)$p.value
  }else{
    p_1 <- NA
    p_2 <- NA
  }
  fc <- log2(mean(input_a$rpkm.x)/mean(input_b$rpkm.x))
  delta_pau <- median(input_a$pau)-median(input_b$pau)
  pau_a=median(input_a$pau)
  pau_b=median(input_b$pau)
  out <- list("p_expr"=p_1,"p_pau"=p_2,"fc"=fc,"delta_pau"=delta_pau,"pau_a"=pau_a,"pau_b"=pau_b)
  return(out)
}

GetPvalue(args[2])
# sapply(gene_lst,GetPvalue) %>% apply(MARGIN = 2,unlist) %>% t() %>% data.frame()%>%
#   rownames_to_column(.,"gene_name") %>% na.omit() -> res_f
# 
# saveRDS(res_f,file = "output/R/wilcox.rds")