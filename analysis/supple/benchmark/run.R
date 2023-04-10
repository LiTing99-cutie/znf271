
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

res <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/intersect.txt")

data.frame(name=c("PolyASite","DaPars2"),per=c(res$inter_2_1/res$n_2,res$inter_2_3/res$n_2)) -> i

p <- ggplot(i,aes(x=name,y=per*100))+geom_col(fill=c("#f2a93b","78b6f9"),width = 0.5)+theme_1()+xlab(NULL)+ylab("Propotion (%)")+
  geom_text(aes(label=paste0(round(per,3)*100,"%")),size=4,vjust=-0.5)+
  ggtitle("Iso-seq PASs identified by other methods")
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/output/compare.pdf",width = 5,height = 5)

db <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/benchmark/atlas.clusters.2.0.GRCh38.96.bed")

summary(db$V5)
summary(db$V7)
summary(db$V8)
