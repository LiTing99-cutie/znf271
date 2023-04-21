source("/home/user/data2/lit/bin/lit_utils.R")
lib()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/wilcox.txt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/disrupt_cds_pro_pa.txt"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}

#### Read wilcox res ####
fread(args[1]) -> res
typeII <- fread(args[2])

#### FDR ####
p <- res$p_expr
p.adjust(p,"BH")->p.adust
mutate(res,p_expr.adjust=p.adust) ->res
res %>% mutate(diff_pau=case_when(abs(delta_pau)>=0.2~"Diff",
                                  abs(delta_pau)<0.2~"NotDiff"),
               diff_expr=case_when(abs(fc)>=1 & p_expr.adjust<=0.05~"Diff",
                                   (abs(fc)<1 | p_expr.adjust>0.05)~"NotDiff"))->res
res$cds[res$gene_name %in% typeII$gene_name] <- "Disrupt cds"
res$cds[!res$gene_name %in% typeII$gene_name] <- "Do dot disrupt cds"
#### write ####
filepath <- paste0(args[3],"/wilcox.p_adjust.diff.cds.txt")
fwrite_c(res,filepath)

#### Remove outliers of delta_pau to draw plots ####
fun.outlier<- function(x,time.iqr=1.5) {
  outlier.low<-quantile(x,probs=c(0.25))-IQR(x)*time.iqr
  outlier.high<-quantile(x,probs=c(0.75))+IQR(x)*time.iqr
  x[which(x>outlier.high | x<outlier.low)]<-NA
  x
}
res$delta_pau <- fun.outlier(res$delta_pau)
res <- na.omit(res)
res %>% mutate(type=case_when(
  diff_expr =="Diff" ~ "expr_c",
  diff_expr =="NotDiff" & diff_pau=="Diff" ~"expr_nc_pau_c",
  diff_expr =="NotDiff" & diff_pau=="NotDiff" ~"expr_nc_pau_nc"
)) ->res
res$type <- factor(res$type,levels = c("expr_c","expr_nc_pau_c","expr_nc_pau_nc"))
res %>% 
  ggplot(aes(x=delta_pau,y=fc,color=type))+
  geom_point(size=1.5,alpha=1)+theme_pubr()+
  scale_color_manual(values=c(colorspace::lighten("gray",0.3),"#bd6a9a",colorspace::lighten("gray",0.3)))+
  xlab("PA usage change")+
  ylab("Log2 (fold change)")+guides(color = "none")+
  xlim(-0.6,0.6)+
  geom_hline(aes(yintercept=c(1)),linetype="dashed")+
  geom_hline(aes(yintercept=c(-1)),linetype="dashed")+
  geom_vline(aes(xintercept=c(-0.2)),linetype="dashed")+
  geom_vline(aes(xintercept=c(0.2)),linetype="dashed") -> p
  
filepath=paste0(args[3],"/ScatterPlot.pdf")
ggsave(p,filename = filepath,width = 5,height = 5)



             