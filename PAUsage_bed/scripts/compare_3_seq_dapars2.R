args=commandArgs(T)
# 1 input file
# 2 output path
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

res <- fread(args[1])

data.frame(name=c("PolyASite","DaPars2"),per=c(res$inter_2_1/res$n_2,res$inter_2_3/res$n_2)) -> i

p <- ggplot(i,aes(x=name,y=per*100))+geom_col(fill=c("#f2a93b","78b6f9"),width = 0.5)+theme_1()+xlab(NULL)+ylab("Propotion (%)")+
  geom_text(aes(label=paste0(round(per,3)*100,"%")),size=4,vjust=-0.5)+
  ggtitle("Iso-seq PASs identified by other methods")
filepath=paste0(args[2],"/compare.pdf")
ggsave(p,filename = filepath,width = 5,height = 5)
