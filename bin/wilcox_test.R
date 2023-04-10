args <- commandArgs(TRUE)
#### Library packages ####
suppressMessages(library(tidyfst,quietly = T))
suppressMessages(library(tibble,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(ggpubr,quietly = T))
suppressMessages(library(ggrepel,quietly = T))
suppressMessages(library(data.table,quietly = T))

#### Read rpkm ####
rpkm <- fread(args[1])
rpkm <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt")
colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]
print("Gene number")
nrow(rpkm)/2

#### Read metadata ####
metadata=fread("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.txt",sep='\t')
frag_score=fread("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt",sep='\t')
period=fread("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period.txt",sep='\t')

#### Merge ####
merge(metadata,frag_score,"sample") %>% filter(frag_score>0.885) %>%
  merge(period,"Developmental_Stage")->md
merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(md,"sample") %>%
  select(gene_id,rpkm.x,rpkm.y,period.abbre.abbre) %>% filter(rpkm.x>1) %>% mutate(pau=rpkm.y/rpkm.x)-> res

#### Prepare gene list ####
unique(res$gene_id) -> gene_lst

GetPvalue <- function(gene){
  res %>% filter(gene_id==gene) ->input
  input_b <- filter(input,period.abbre.abbre=="Before birth")
  input_a <- filter(input,period.abbre.abbre=="After birth")
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

sapply(gene_lst,GetPvalue) %>% apply(MARGIN = 2,unlist) %>% t() %>% data.frame()%>%
  rownames_to_column(.,"gene_name") %>% na.omit() -> res_f

saveRDS(res_f,file = "output/R/wilcox.rds")

###### ZNF271P #####
# res %>% filter(gene_id=="ZNF271P") ->input
# input_b <- filter(input,period.abbre.abbre=="Before birth")
# input_a <- filter(input,period.abbre.abbre=="After birth")
# median(input_b$pau)
# # 0.9684989
# median(input_a$pau)
# # 0.7109515