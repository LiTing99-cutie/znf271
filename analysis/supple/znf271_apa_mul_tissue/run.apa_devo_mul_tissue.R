args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib()
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
# rpkm_path <- args[1]
rpkm_path <- "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue/rpkm/human.all_tissue.stringtie.rpkm.txt"
rpkm <- fread(rpkm_path)
colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]

#### Merge with metadata ####
metadata_path <- "/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/human.metadata.clean.txt"
metadata <- fread(file=metadata_path)

merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(metadata,"sample") %>%
  select(gene_id,rpkm.x,rpkm.y,stage,organ) %>% filter(rpkm.x>1) %>% mutate(pau=rpkm.y/rpkm.x)-> res

res %<>%
  mutate(organ = case_when(
    grepl("brain", organ) ~ "brain",
    TRUE ~ organ
  ))

filter(res,organ=="heart") -> res.heart
filter(res,organ=="kidney") -> res.kidney
filter(res,organ=="liver") -> res.liver
filter(res,organ=="testis") -> res.testis 
filter(res,organ=="brain") -> res.brain
l <- list(res.heart,res.kidney,res.liver,res.testis,res.brain)
#### Case delta pau and p_value in multiple tissues ####
args[3] <- "ZNF271P"
GetPvalue <- function(df,gene){
  df %>% filter(gene_id==gene) ->input
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

lapply(l,GetPvalue,gene=args[3]) %>% unlist() %>% .[c(4,10,16,22,28)] -> delta_pau
lapply(l,GetPvalue,gene=args[3]) %>% unlist() %>% .[c(2,8,14,20,26)] -> p_value

df <- data.frame(organ=c("Heart","Kidney","Liver","Testis","Brain"),delta_pau=delta_pau,p_value=p_value)

fwrite(df,file='/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue/output/mul_tissue_delta_pau_p_value.txt')
#### plot ####

firstup <- function(x){
  substr(x,1,1) <- toupper(substr(x,1,1))
  x
}

plot_custom <- function(organ,gene){
  organ_1 <- enquo(organ)
  gene_1 <- enquo(gene)
  res %>% filter(gene_id==UQ(gene_1)) %>% filter(organ==UQ(organ_1)) -> ggplot.input
  ggplot.input$organ <- firstup(ggplot.input$organ)
  ggplot.input$stage[ggplot.input$stage=="embryo"] <- "Embryo"
  ggplot.input$stage[ggplot.input$stage=="postnatal"] <- "Postnatal"
  ggplot.input$stage <- factor(ggplot.input$stage,levels = c("Embryo","Postnatal"))

  theme <- theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(size = 15),
                 axis.title = element_text(size = 16),
                 legend.title = element_text(size = 14),
                 legend.text=element_text(size = 13),
                 plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
                 plot.title=element_text(hjust=0.5)
                 # legend.position=c(0.9,0.9)
                 )

  compare <- function(p){
    my_comparisons <- list(c("Embryo","Postnatal"))
    res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"))
    res
  }

  output_path=paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue/output/",organ,"/")
  ifelse(!dir.exists(output_path),dir.create(output_path),FALSE)
  filename2=paste0(output_path,"rpkm.pro.2.pdf")
  filename3=paste0(output_path,"rpkm.dis.2.pdf")
  filename4=paste0(output_path,"pau.2.pdf")
  
  plot_cus_1_pau <- function(filename,y,ylab,ylim){
    p <- ggboxplot(ggplot.input,x="stage",y=y,fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
                   palette = c('#f8d396','#a3bdd8'))+
      geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
      labs(x="Developmental stage",y = ylab)+
      ggtitle(firstup(organ))+
      guides(fill=guide_legend(title = NULL))+
      theme+ylim(0,ylim)
    p1 <- compare(p)
    ggsave(filename,p1,height = 5,width = 5)
    return(p1)
  }
  
  plot_cus_1 <- function(filename,y,ylab){
    p <- ggboxplot(ggplot.input,x="stage",y=y,fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
                   palette = c('#f8d396','#a3bdd8'))+
      geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
      labs(x="Developmental stage",y = ylab)+
      ggtitle(firstup(organ))+
      guides(fill=guide_legend(title = NULL))+
      theme
    p1 <- compare(p)
    ggsave(filename,p1,height = 5,width = 5)
    return(p1)
  }

  plot_cus_1(filename2,"rpkm.x","Proximal RPKM")
  plot_cus_1(filename3,"rpkm.y","Distal RPKM")
  plot_cus_1_pau(filename4,"pau","PA usage",1.7) -> p
  
  return(p)
}

l_organ <- list("brain","heart","kidney","liver","testis")
lapply(l_organ,plot_custom,gene="ZNF271P") -> mul_organ_p

mul_organ_p[[1]]+mul_organ_p[[2]]+mul_organ_p[[3]]+mul_organ_p[[4]]+mul_organ_p[[5]]+plot_annotation(tag_levels="A") -> p

ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue/output/mul_tissue_delta_pau_p_value.pdf",
       width = 10,height = 10)
