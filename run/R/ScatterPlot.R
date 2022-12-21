
library(data.table)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(tidyr)
library(dplyr)

# 4877
res <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/final_res_ref_based_all_1.with_geneName.txt")
spe_type <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt")

res$x=res$pau_a-res$pau_b
res$y=log2(res$rpkm_pro_a/res$rpkm_pro_b)
res %<>% mutate(.,sig=case_when(pvalue<=0.001~"Sig",pvalue>0.001~"NotSig"))
res %<>% mutate(.,sig_1=case_when(p_value_1<=0.001~"Sig",p_value_1>0.001~"NotSig"))

res <- unite(res,"SigAll",sig,sig_1,remove = FALSE)

# 4875
res <- na.omit(res)

# 3215
res[res$pau_a<=1 & res$pau_a>=-1 & res$pau_b>=-1 & res$pau_b<=1] -> res

res$label[res$gene_name %in% spe_type$gene_name] <- res[res$gene_name %in% spe_type$gene_name]$gene_name

theme <-theme_pubr()+theme(axis.text = element_text(size = 12),
                   axis.title = element_text(size = 14))

res %>% 
  ggplot(aes(x=x,y=y,color=SigAll))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),"#857fb6", "#bd6a9a", "#857fb6"))+
  xlab("PA usage change")+
  ylab("Gene expression fold change")+guides(color = "none")+geom_label_repel(aes(label=res$label),size=3,max.overlaps=50,
                                                                              min.segment.length = 0)->p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.pdf",width = 5,height = 5)
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.svg",width = 5,height = 5)
# not supported
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.eps",width = 5,height = 5)
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.tiff",width = 5,height = 5)

pdf(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.1.pdf",width = 5,height = 5)
p
dev.off()


res %>% 
  ggplot(aes(x=x,y=y,color=SigAll))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),colorspace::lighten("gray",0.3), "#bd6a9a", colorspace::lighten("gray",0.3)))+
  xlab("PA usage change")+
  ylab("Gene expression fold change")+guides(color = "none")+geom_label_repel(aes(label=res$label),size=3,max.overlaps=50,
                                                                              min.segment.length = 0)->p
pdf(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/ScatterPlot.2.pdf",width = 5,height = 5)
p
dev.off()

table(res$sig,res$sig_1)

saveRDS(res,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rds")



