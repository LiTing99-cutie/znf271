
args <- commandArgs(TRUE)
# 1 gene 
# rm(list=ls())

# setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

suppressMessages(library(data.table,quietly=TRUE))
suppressMessages(library(dplyr,quietly=TRUE))
suppressMessages(library(magrittr,quietly=TRUE))
suppressMessages(library(ggplot2,quietly=TRUE))
suppressMessages(library(reshape2,quietly=TRUE))
suppressMessages(library(ggpubr,quietly=TRUE))

# read metadata 
metadata <- readRDS(file="/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.rds")
period <- readRDS(file="/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period_description.rds")

# 23/3/27 all_1 to all_1_loose
file=paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/stringtie_whole_genome_ref_based_all_1_loose/stringtie.rpkm.",args[1],".txt")
rpkm <- fread(file)

rpkm.proximal <- rpkm[grep("proximal",rpkm$V1),]
colnames(rpkm.proximal) <- c("type","gene_id","rpkm_proximal","sample")
rpkm.distal <- rpkm[grep("distal",rpkm$V1),]
colnames(rpkm.distal) <- c("type","gene_id","rpkm_distal","sample")

pau <- data.frame(pau=rpkm.distal$rpkm_distal/rpkm.proximal$rpkm_proximal,sample=rpkm.distal$sample)

merge(rpkm.proximal[,3:4],rpkm.distal[,3:4],"sample") %>% merge(.,pau,"sample") %>%
  merge(.,metadata,"sample") %>% merge(.,period,"Developmental_Stage")->ggplot.input

theme <- theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_text(size = 15),
               axis.title = element_text(size = 16),
               legend.title = element_text(size = 14,face = "bold"),
               legend.text=element_text(size = 13),
               plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
               legend.position=c(0.9,0.9))

output_path=paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/",args[1],"/")
# system("[ -d output_path ] || mkdir -p output_path")
filename1=paste0(output_path,"pau.5.pdf")
filename2=paste0(output_path,"rpkm.pro.2.pdf")
filename3=paste0(output_path,"rpkm.dis.2.pdf")
filename4=paste0(output_path,"pau.2.pdf")

pdf(file = filename1,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="period",y="pau",color="period",outlier.shape = NA,width = 0.5,lwd=1,
          palette = c('#99d8c9','#66c2a4','#41ae76','#238b45','#005824'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "PA usage")+
  guides(fill=guide_legend(title = NULL))

my_comparisons <- list(c("Embryonic","Early fetal"),c("Early fetal","Early mid-fetal"),c("Embryonic","Early mid-fetal"))
p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                      method = "wilcox.test",method.args = list(alternative = "two.sided"))+
  theme(plot.margin = unit(c(0.5,1.5,0.5,0.5),"cm"))
dev.off()

pdf(file = filename2,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="period.abbre.abbre",y="rpkm_proximal",fill="period.abbre.abbre",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "Proximal RPKM")+
  guides(fill=guide_legend(title = NULL))+theme

my_comparisons <- list(c("Before birth","After birth"))
p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                      method = "wilcox.test",method.args = list(alternative = "two.sided"))+theme
dev.off()

pdf(file = filename3,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="period.abbre.abbre",y="rpkm_distal",fill="period.abbre.abbre",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "Distal RPKM")+
  guides(fill=guide_legend(title = NULL))+theme

my_comparisons <- list(c("Before birth","After birth"))
p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                      method = "wilcox.test",method.args = list(alternative = "two.sided"))+theme
dev.off()

pdf(file = filename4,height = 5,width = 5)
p <- ggboxplot(ggplot.input,x="period.abbre.abbre",y="pau",fill="period.abbre.abbre",outlier.shape = NA,width = 0.6,lwd=1,
               palette = c('#f8d396','#a3bdd8'))+
  geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0))+
  labs(x="Developmental stage",y = "PA usage")+
  guides(fill=guide_legend(title = NULL))+theme

my_comparisons <- list(c("Before birth","After birth"))
p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                      method = "wilcox.test",method.args = list(alternative = "two.sided"))+theme
dev.off()


