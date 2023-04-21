source("/home/user/data2/lit/bin/lit_utils.R")
lib()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/wilcox.p_adjust.diff.cds.txt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

# read data in 
res=fread(args[1])

df1 <- data.frame(table(res$diff_expr))
df1$Var1 <- factor(df1$Var1,levels = c("Diff","NotDiff"))
df1$per <- df1$Freq/sum(df1$Freq)
labs <- paste0(df1$Freq," (",round(df1$per*100,2),"%",")")
p <- df1 %>% 
  ggpie("Freq",
        fill = "Var1", color = "black",
        palette = c("#857fb6", colorspace::lighten("gray",0.3)),
        lab.pos = "in", lab.font = c("bold","black"),
        label = labs)+
  theme(plot.title = element_text(hjust=0.5,vjust = -80))+
  labs(fill="Expression Change")
filepath=paste0(args[2],"/pie1.pdf")
ggsave(p,filename=filepath,width = 5,height = 5)

df2 <- data.frame(table(res[res$diff_expr=="NotDiff",]$diff_pau))
df2$Var1 <- factor(df2$Var1,levels = c("NotDiff","Diff"))
df2$per <- df2$Freq/sum(df2$Freq)
labs <- paste0(df2$Freq," (",round(df2$per*100,2),"%",")")
p <- df2 %>% 
  ggpie("Freq",
        fill = "Var1", color = "black",
        palette = c(colorspace::lighten("gray",0.3),"#bd6a9a"),
        lab.pos = "in", lab.font = c("bold","black"),
        label = labs)+
  theme(plot.title = element_text(hjust=0.5,vjust = -80))+
  labs(fill="PAU Change")
filepath=paste0(args[2],"/pie2.pdf")
ggsave(p,filename=filepath,width = 5,height = 5)

df <- data.frame(table(res[res$diff_expr=="NotDiff" & res$diff_pau=="Diff",]$cds))
df$Var1 <- factor(df$Var1,levels = c("Disrupt cds","Do dot disrupt cds"))
df$per <- df$Freq/sum(df$Freq)
labs <- paste0(df$Freq," (",round(df1$per*100,2),"%",")")
p <- df %>% 
  ggpie("Freq",
        fill = "Var1", color = "black",
        palette = c("#bd6a9a","#817bae","#b5b6b6"),
        lab.pos = "in", lab.font = c("bold","black"),
        label = labs)+
  theme(plot.title = element_text(hjust=0.5,vjust = -80))+
  labs(fill="Effect On Coding")
filepath=paste0(args[2],"/pie3.pdf")
ggsave(p,filename=filepath,width = 5,height = 5)
