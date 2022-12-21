
rm(list=ls())
library(tidyfst)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(stringr)

setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

# res=fread("output/final_list/final_res_ref_based_all_1.with_geneName.txt")
# 2022_12_19 remove genes whose pau is wierd
res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rds")

res %>% mutate(diff_pa_usage=case_when(
  pvalue>0.001~"No change",
  pvalue<=0.001 & pau_a>pau_b~"Lengthening",
  pvalue<=0.001 & pau_a<pau_b~"Shortening"
)) %>% 
  mutate(diff_rpkm=case_when(
    p_value_1>0.001~"No change",
    p_value_1<=0.001 & rpkm_pro_a>rpkm_pro_b~"Up",
    p_value_1<=0.001 & rpkm_pro_a<rpkm_pro_b~"Down"
  )) -> res.mutate


tab=as.data.frame(table(res.mutate$diff_pa_usage,res.mutate$diff_rpkm))

tab %>% arrange(Var1,desc(Var2)) %>% group_by(Var1) %>% mutate(ylabel=cumsum(Freq)-0.5*Freq) -> tab_1

pdf("output/R/hist.pdf")
tab_1 %>% ggplot(mapping=aes(x=Var1,y=Freq,fill=Var2))+geom_col(width = 0.6)+scale_fill_tableau()+
  theme_pubr()+theme(axis.text = element_text(size = 15),
                     axis.title = element_text(size = 16,face = "bold"),
                     legend.title = element_text(size = 14,face = "bold"),
                       legend.text=element_text(size = 13))+xlab(NULL)+ylab("Count")+labs(fill="Expression")+
  geom_text(mapping = aes(y=ylabel,label=paste(Freq)),size=2.5)
dev.off()

filter(res.mutate,diff_pa_usage!="No change",diff_rpkm=="No change") -> final.list

write.table(final.list[,c("gene_name","gene_id","gene_type","diff_pa_usage","diff_rpkm")],file = "output/final_list/diff.lst",
            quote = FALSE,row.names = FALSE,sep = '\t')

# hist(as.data.frame(table(final.list$gene_type)),xlab = "PA usage after birth",ylab="Count",main=NULL,col = "#2b83ba")
tmp <- as.data.frame(table(final.list$gene_type))

pseu_or_not=c()
for (i in 1:length(tmp$Var1)){
  target=as.character(tmp$Var1[i])
  pseu_or_not=c(pseu_or_not,str_detect(target,"pseudogene"))
}

pseu_sum=sum(tmp[pseu_or_not,2])

tmp <- rbind(tmp[!pseu_or_not,],data.frame("Var1"="Pseudogene","Freq"=pseu_sum))

p <- tmp %>% mutate(type=c("lncRNA","Coding","Pseudogene")) %>% 
  ggplot(aes(x=type,y=Freq))+
  geom_bar(stat = 'identity',width=0.5)+labs(y = "Count",x='') + guides(fill = "none")+
  geom_text(aes(label = Freq),size=5,vjust=-0.3)+
  theme_pubr()+theme(axis.text.x = element_text(size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 16,face = "bold"),
                     legend.title = element_text(size = 14,face = "bold"),
                     legend.text=element_text(size = 13))


ggsave("output/R/final_gene_type.pdf",p,width = 5,height = 5)
