
args <- commandArgs(TRUE)
# 1 gene 

suppressMessages(library(data.table,quietly=TRUE))
suppressMessages(library(dplyr,quietly=TRUE))
suppressMessages(library(magrittr,quietly=TRUE))
suppressMessages(library(ggplot2,quietly=TRUE))
suppressMessages(library(reshape2,quietly=TRUE))
suppressMessages(library(ggpubr,quietly=TRUE))

# read file in
metadata_path <- paste0("/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/",args[1],".metadata.clean.txt")
metadata <- fread(file=metadata_path)
rpkm_path <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare/",args[1],"/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.txt")
rpkm <- fread(rpkm_path)

colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]

# add %>% filter(gene_id==args[2]) %>% 
merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(metadata,"sample") %>% filter(gene_id==args[2]) %>% 
  select(gene_id,rpkm.x,rpkm.y,stage) %>% filter(rpkm.x>1) %>% mutate(pau=rpkm.y/rpkm.x)-> ggplot.input

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
               legend.position=c(0.9,0.9))

compare <- function(p){
  my_comparisons <- list(c("Embryo","Postnatal"))
  res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                        method = "wilcox.test",method.args = list(alternative = "two.sided"))
  res
}

output_path=paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/devo_compare/output/",args[1],"/")
ifelse(!dir.exists(output_path),dir.create(output_path),FALSE)
filename2=paste0(output_path,"rpkm.pro.2.pdf")
filename3=paste0(output_path,"rpkm.dis.2.pdf")
filename4=paste0(output_path,"pau.2.pdf")

pdf(file = filename2,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="stage",y="rpkm.x",fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "Proximal RPKM")+
  guides(fill=guide_legend(title = NULL))+theme
compare(p)
dev.off()

pdf(file = filename3,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="stage",y="rpkm.y",fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "Distal RPKM")+
  guides(fill=guide_legend(title = NULL))+theme
compare(p)
dev.off()

pdf(file = filename4,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="stage",y="pau",fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "PA usage")+
  guides(fill=guide_legend(title = NULL))+theme
compare(p)
dev.off()