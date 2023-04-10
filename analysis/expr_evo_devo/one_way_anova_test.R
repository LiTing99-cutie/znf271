suppressMessages(library(tidyfst,quietly = T))
suppressMessages(library(tibble,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(ggpubr,quietly = T))
suppressMessages(library(ggrepel,quietly = T))
suppressMessages(library(data.table,quietly = T))

read_rpkm <- function(spe){
  rpkm_path <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/expr_evo_devo/rpkm/",spe,"/stringtie.rpkm.txt")
  rpkm <- fread(rpkm_path)
  colnames(rpkm) <- c("gene_name","rpkm","sample")
  return(rpkm)
}

h <- read_rpkm("human")
r <- read_rpkm("rhesus")
m <- read_rpkm("mouse")

r <- r[!r$gene_name =="."]

ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)
# change m and r gene_name with ortholog name in h
merge(m,ortho,by.x = "gene_name",by.y = "V3") %>% mutate(gene_name=NULL,V2=NULL) %>% rename(gene_name=V1) %>% select(gene_name,rpkm,sample)-> m_rn
merge(r,ortho,by.x = "gene_name",by.y = "V2",allow.cartesian = T) %>% mutate(gene_name=NULL,V3=NULL) %>% rename(gene_name=V1) %>% select(gene_name,rpkm,sample)-> r_rn

intersect(h$gene_name,m_rn$gene_name) %>% intersect(.,r_rn$gene_name) -> gene_lst

GetPvalue <- function(gene){
   h %>% filter(gene_name==gene) ->h_i
  r_rn %>% filter(gene_name==gene) ->r_i
  m_rn %>% filter(gene_name==gene) ->m_i
  df_1 <- data.frame(rpkm=h_i$rpkm,spe="human")
  df_2 <- data.frame(rpkm=r_i$rpkm,spe="rhesus")
  df_3 <- data.frame(rpkm=m_i$rpkm,spe="mouse")
  df <- rbind(df_1,df_2,df_3)
  oneway.test(rpkm~spe,data=df,var.equal = T) ->test
  l <- list(mean(df_1$rpkm),mean(df_2$rpkm),mean(df_3$rpkm),test$p.value)
  return(l)
}

sapply(gene_lst,GetPvalue) %>% apply(MARGIN = 2,unlist) %>% t() %>% data.frame()%>%
  rownames_to_column(.,"gene_name") %>% rename(h_mean=X1,r_mean=X2,m_mean=X3,p_value=X4) -> res

saveRDS(res,file = "output/anova.rds")
