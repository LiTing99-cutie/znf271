source("/home/user/data2/lit/bin/lit_utils.R")
lib()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/stringtie/stringtie.rpkm.txt"
  args[2] <- "/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/human.metadata.clean.txt"
  args[3] <- "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt"
  args[4] <- "ZNF271P"
  args[5] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}
#### Read rpkm ####
rpkm <- fread(args[1])
colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]
#### Read metadata ####
metadata=fread(args[2])
frag_score=fread(args[3],sep='\t')
colnames(frag_score)=c("sample","frag_score")
#### Merge ####
merge(metadata,frag_score) %>% filter(frag_score>0.885) -> md
merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(md,"sample") %>%
  select(gene_id,rpkm.x,rpkm.y,stage) %>% filter(rpkm.x>=1) %>% mutate(pau=rpkm.y/rpkm.x)-> res
filter(res,gene_id==args[4]) %>% as.data.frame()-> ggplot.input
ggplot.input$stage <- firstup(ggplot.input$stage)

theme <- theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_text(size = 15),
               axis.title = element_text(size = 16),
               legend.title = element_text(size = 14,face = "bold"),
               legend.text=element_text(size = 13),
               plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
               legend.position=c(0.9,0.9))

compare <- function(p){
  my_comparisons <- list(c("Embryo","Postnatal"))
  res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                               method = "wilcox.test",method.args = list(alternative = "two.sided"))
  res
}

filename2=paste0(args[5],"rpkm.pro.2.pdf")
filename3=paste0(args[5],"rpkm.dis.2.pdf")
filename4=paste0(args[5],"pau.2.pdf")

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


