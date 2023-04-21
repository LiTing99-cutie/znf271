source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib()
library(bedtoolsr)
args <- commandArgs(T)
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/human/PAusage.bed6+"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}

PAS <- fread(args[1])

gene_number=data.frame()
pas_number=data.frame()
for (i in 1:20){
  PAS[PAS$V5>=i] %>% count(V4) -> df
  gene_number=rbind(gene_number,data.frame(n=nrow(df),cutoff=i))
  pas_number=rbind(pas_number,data.frame(n=df$n,cutoff=i))
}


pas_number$cutoff <- factor(pas_number$cutoff)

filepath=paste0(args[2],"/cutoff_gene_number.pdf")
ggplot(gene_number,aes(x=cutoff,y=n))+geom_point()+xlab("Cutoff")+ylab("Gene Number")
ggsave(file=filepath,width = 5,height = 5)

filepath=paste0(args[2],"/cutoff_pas_number.pdf")
ggplot(pas_number,aes(x=cutoff,y=n,fill=cutoff))+geom_boxplot(outlier.shape=NA)+xlab("Cutoff")+ylab("PAS Number")+
  stat_summary(fun="mean",geom="point",size=2.5,color="red",fill="red")+
  scale_y_continuous(breaks=seq(0,12,2))+coord_trans(x = "identity", y = "identity", xlim = NULL, ylim = c(0,12))
ggsave(file=filepath,width = 5,height = 5)

