
# Use a different differential res

#### Library packages ####
rm(list=ls())
suppressMessages(library(tidyfst,quietly = T))
suppressMessages(library(tibble,quietly = T))
suppressMessages(library(dplyr,quietly = T))
suppressMessages(library(ggplot2,quietly = T))
suppressMessages(library(ggpubr,quietly = T))
suppressMessages(library(ggrepel,quietly = T))
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1")

#### Read wilcox res ####
readRDS(file = "output/R/wilcox.rds") -> res
print("Gene number after filtering based on proximal rpkm >1 and sample number>=8")
nrow(res)

#### Remove outliers of pau_a and pau_b (remove dots like 221.929195) ####
fun.outlier<- function(x,time.iqr=1.5) {
  outlier.low<-quantile(x,probs=c(0.25))-IQR(x)*time.iqr
  outlier.high<-quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  x[which(x>outlier.high | x<outlier.low)]<-NA
  x
}
res$pau_a <- fun.outlier(res$pau_a)
res$pau_b <- fun.outlier(res$pau_b)
res <- na.omit(res)
print("Gene number after filtering pau_a or pau_b outliers")
nrow(res)

#### FDR ####
p <- res$p_expr
p.adjust(p,"BH")->p.adust
mutate(res,p_expr.adjust=p.adust) ->res
res %>% mutate(diff_pau=case_when(abs(delta_pau)>=0.2~"Diff",
                                  abs(delta_pau)<0.2~"NotDiff"),
               diff_expr=case_when(fc>=1 & p_expr.adjust<0.05~"UP",
                                   !(abs(fc)>=1 & p_expr.adjust<0.05)~"NotDiff",
                                   fc<=-1 & p_expr.adjust<0.05~"DOWN"))->res
saveRDS(res,file = "output/final_list/res.rm_outlier.rds")
write.table(res,file = "output/final_list/res.rm_outlier.txt",sep = '\t',quote = F,row.names = F,col.names = T)
#### Remove outliers of delta_pau to draw plots ####
res$delta_pau <- fun.outlier(res$delta_pau)
res <- na.omit(res)
print("Gene number after filtering delta pau (only filter while ploting scatter plot)")
nrow(res)
res$diff_expr <- factor(res$diff_expr,levels = c("UP","NotDiff","DOWN"))
res %>% mutate(type=case_when(
  diff_expr !="NotDiff" ~ "expr_c",
  diff_expr =="NotDiff" & diff_pau=="Diff" ~"expr_nc_pau_c",
  diff_expr =="NotDiff" & diff_pau=="NotDiff" ~"expr_nc_pau_nc"
)) ->res
res$type <- factor(res$type,levels = c("expr_c","expr_nc_pau_c","expr_nc_pau_nc"))
res$label[res$gene_name=="ZNF271P"] <- "ZNF271P"
res %>% 
  ggplot(aes(x=delta_pau,y=fc,color=type))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),"#bd6a9a",colorspace::lighten("gray",0.3)))+
  xlab("PA usage change")+
  ylab("Log2(fold change)")+guides(color = "none")+
  xlim(-0.6,0.6)+
  geom_hline(aes(yintercept=c(1)),linetype="dashed")+
  geom_hline(aes(yintercept=c(-1)),linetype="dashed")+
  geom_vline(aes(xintercept=c(-0.2)),linetype="dashed")+
  geom_vline(aes(xintercept=c(0.2)),linetype="dashed") -> p
  # geom_label_repel(aes(label=res$label),min.segment.length = 0,color="black",size=3) 
  

ggsave(p,filename = "output/R/ScatterPlot.c.pdf",width = 5,height = 5)

             