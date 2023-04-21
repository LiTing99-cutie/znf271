source("/home/user/data2/lit/bin/lit_utils.R")
lib()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output//wilcox.p_adjust.diff.txt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

options(warn = -1)

res=fread(args[1])


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
filepath=paste0(args[2],"/GO.txt")
write.table(enrich,file = filepath,sep = '\t',quote = F,row.names = F)

###### Custom Plot ######
GO <- fread(filepath)
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
filepath=paste0(args[2],"/GO.pdf")
pdf(file = filepath,width = 6.5,height = 4.5)
p
invisible(dev.off())

