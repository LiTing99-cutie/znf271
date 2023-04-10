# Use a different differential res

options(warn = -1)
# rm(list=ls())
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")
suppressMessages(library("clusterProfiler"))
suppressMessages(library("DOSE"))
suppressMessages(library("stringr"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
typeIIAndIII <- fread(file="output/final_list/typeIIAndIII.txt")
###### expression change significantly during development ######
res_fil <- filter(res,diff_expr=="UP" | diff_expr=="DOWN")
enrich <- enrichGO(gene=res_fil$gene_name,
                   OrgDb="org.Hs.eg.db",
                   ont="ALL",
                   keyType = "SYMBOL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
pdf(file = "output/R/GO_c/dotplot.GO.pdf",width = 8,height = 6)
dotplot(enrich)
invisible(dev.off())
write.table(enrich,file = "output/R/GO_c/GO.txt",sep = '\t',quote = F,row.names = F)
###### PA change significantly during development ######
res_fil <- filter(res,diff_pau=="Diff")
enrich <- enrichGO(gene=res_fil$gene_name,
                   OrgDb="org.Hs.eg.db",
                   ont="ALL",
                   keyType = "SYMBOL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
pdf(file = "output/R/GO_c/dotplot.GO.PA_S.pdf",width = 8,height = 6)
dotplot(enrich)
invisible(dev.off())
write.table(enrich,file = "output/R/GO_c/GO.PA_S.txt",sep = '\t',quote = F,row.names = F)
###### [chosen to display] PA change significantly during development (expression do not change significantly) ######
res_fil <- filter(res,diff_expr=="NotDiff" & diff_pau=="Diff")
enrich <- enrichGO(gene=res_fil$gene_name,
                   OrgDb="org.Hs.eg.db",
                   ont="ALL",
                   keyType = "SYMBOL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
pdf(file = "output/R/GO_c/dotplot.GO.N+S.pdf",width = 6,height = 5)
dotplot(enrich)
invisible(dev.off())
write.table(enrich,file = "output/R/GO_c/GO.N+S.txt",sep = '\t',quote = F,row.names = F)
###### Disrupt CDS ######
# genelist <- bitr(typeIIAndIII$gene_name,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
# head(genelist,2)
# enrich <- enrichGO(gene=typeIIAndIII$gene_name,
#                    OrgDb="org.Hs.eg.db",
#                    ont="ALL",
#                    keyType = "SYMBOL",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff = 0.05,
#                    qvalueCutoff = 0.05,
#                    readable = TRUE)
# pdf(file = "output/R/GO_c/dotplot.GO.dis.pdf",width = 8,height = 6)
# dotplot(enrich)
# dev.off()
# write.table(enrich,file = "output/R/GO_c/GO.dis.txt",sep = '\t',quote = F,row.names = F)


###### Custom Plot ######
GO <- fread("output/R/GO_c/GO.N+S.txt")
GO <- head(GO,10)
gr1 <- as.numeric(str_split(GO$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(GO$GeneRatio,"/",simplify = T)[,2])
GO$ratio <- gr1/gr2
GO <- dplyr::arrange(GO,ratio)
GO$Description <- factor(GO$Description,levels = GO$Description,ordered = T)
settings.1 <- element_text(size = 14,color = "black")
settings.2 <- element_text(size = 12,color = "black")
settings.3 <- element_text(size = 10,color = "black")
p <- GO %>% 
  ggplot(aes(x=ratio,y=Description))+
  geom_point(aes(color=p.adjust,size=Count))+
  scale_color_gradient(low="#b46793",high = "#377fae")+
  scale_size_continuous(range = c(4,8))+
  xlab("GeneRatio")+
  ylab(NULL)+
  theme_bw()+
  theme(
    axis.title = settings.1,
    axis.text = settings.2,
    legend.title = settings.2,
    legend.text = settings.3
  )+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 30))
pdf(file = "output/R/GO_c/dotplot.GO.N+S.cus.pdf",width = 6.5,height = 4.5)
p
invisible(dev.off())

GO <- fread("output/R/GO_c/GO.txt")
GO <- head(GO,10)
gr1 <- as.numeric(str_split(GO$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(GO$GeneRatio,"/",simplify = T)[,2])
GO$ratio <- gr1/gr2
GO <- dplyr::arrange(GO,ratio)
GO$Description <- factor(GO$Description,levels = GO$Description,ordered = T)
settings.1 <- element_text(size = 14,color = "black")
settings.2 <- element_text(size = 12,color = "black")
settings.3 <- element_text(size = 10,color = "black")
p <- GO %>% 
  ggplot(aes(x=ratio,y=Description))+
  geom_point(aes(color=p.adjust,size=Count))+
  scale_color_gradient(low="#b46793",high = "#377fae")+
  scale_size_continuous(range = c(4,8))+
  xlab("GeneRatio")+
  ylab(NULL)+
  theme_bw()+
  theme(
    axis.title = settings.1,
    axis.text = settings.2,
    legend.title = settings.2,
    legend.text = settings.3
  )+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 30))
pdf(file = "output/R/GO_c/dotplot.GO.cus.pdf",width = 6.5,height = 4.5)
p
invisible(dev.off())

