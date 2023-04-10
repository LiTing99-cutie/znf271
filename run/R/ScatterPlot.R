
library(data.table)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(tidyr)
library(dplyr)

spe_type <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt")
res <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/final_res_ref_based_all_1_loose.with_geneName.txt")
res$x=res$pau_a-res$pau_b
res$y=log2(res$rpkm_pro_a/res$rpkm_pro_b)
res %<>% mutate(.,sig=case_when(pvalue<=0.05~"Sig",pvalue>0.05~"NotSig"))
res %<>% mutate(.,sig_1=case_when(p_value_1<=0.05~"Sig",p_value_1>0.05~"NotSig"))
res <- unite(res,"SigAll",sig,sig_1,remove = FALSE)
res <- na.omit(res)
res[res$pau_a<=1 & res$pau_a>=-1 & res$pau_b>=-1 & res$pau_b<=1] -> res

res$label[res$gene_name %in% spe_type$gene_name] <- res[res$gene_name %in% spe_type$gene_name]$gene_name
res$type[res$SigAll == "Sig_NotSig"] <- "Sig_NotSig"
res$type[res$SigAll != "Sig_NotSig"] <- "Other"
res$type[res$gene_name %in% spe_type$gene_name] <- "Sig_NotSig_Disrupt CDS"
theme <-theme_pubr()+theme(axis.text = element_text(size = 12),
                   axis.title = element_text(size = 14))
####### sig or not ######
res %>% 
  ggplot(aes(x=x,y=y,color=SigAll))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),"#857fb6", "#bd6a9a", "#857fb6"))+
  xlab("PA usage change")+
  ylab("Gene expression fold change")+guides(color = "none")+geom_label_repel(aes(label=res$label),size=3,max.overlaps=50,
                                                                              min.segment.length = 0)->p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.pdf",width = 6,height = 6)

res %>% 
  ggplot(aes(x=x,y=y,color=SigAll))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),colorspace::lighten("gray",0.3), "#bd6a9a", colorspace::lighten("gray",0.3)))+
  xlab("PA usage change")+
  ylab("Gene expression fold change")+guides(color = "none")+geom_label_repel(aes(label=res$label),size=3,max.overlaps=50,
                                                                              min.segment.length = 0)->p

ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.expr_c_no_c.pdf",width = 6,height = 6)

saveRDS(res,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rds")

####### type ######
res %>% 
  ggplot(aes(x=x,y=y,color=type))+
  geom_point(size=1.5,alpha=0.5)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.7),"#bd6a9a"))+
  xlab("PA usage change")+
  ylab("Gene expression fold change")+guides(color = "none")

hist(filter(res,type=="Sig_NotSig")$y)
hist(filter(res,type=="Sig_NotSig")$x)
