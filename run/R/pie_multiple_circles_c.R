
# rm(list=ls())

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))


res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

if(F){
  df <- data.frame(table(res$sig,res$sig_1))
  
  
  add <- data.frame(table(res$SigAll))
  
  merge(df,add,by="Freq")->df
  
  df$Freq2 <- df$Freq
  df$Freq1 <- df$Freq
  df$Freq <- NULL
  df$Var1 <- factor(df$Var1.y,levels=c("NotSig_Sig","Sig_Sig","NotSig_NotSig","Sig_NotSig"))
  df$Var2 <- factor(df$Var2,levels=c("Sig","NotSig"))
  df$Var1.x <- NULL
  df$Var1.y <- NULL
  df$y <- 'a'
  df$x <- 'b'
  
  # 分别求所占百分比
  dat1 = aggregate(df$Freq1, by = list(df$Var1), FUN = sum)
  dat1$per1 = dat1$x / sum(dat1$x)
  
  # for循环构建标签的相对位置
  for (i in seq(nrow(dat1), 1)) {
    if (i == nrow(dat1)) {
      dat1$per.y1[i] = dat1$per1[i] / 2
    }else{
      dat1$per.y1[i] = sum(dat1$per1[(i + 1):nrow(dat1)]) + dat1$per1[i] / 2
    }
  }
  
  # 构建标签后合并数据
  dat1$label1 = paste(dat1$x,'(',round(dat1$per1*100, 2),'%',')', sep = '')
  df = merge(df, dat1[,c(1,3,4,5)], by.x = 'Var1', by.y = 'Group.1')
  
  
  # 重复操作
  dat2 = aggregate(df$Freq2, by = list(df$Var2), FUN = sum)
  dat2$per2 = dat2$x / sum(dat2$x)
  
  for (i in seq(nrow(dat2), 1)) {
    if (i == nrow(dat2)) {
      dat2$per.y2[i] = dat2$per2[i] / 2
    }else{
      dat2$per.y2[i] = sum(dat2$per2[(i + 1):nrow(dat2)]) + dat2$per2[i] / 2
    }
  }
  
  dat2$label2 = paste(dat2$x,'(',round(dat2$per2*100, 2),'%',')', sep = '')
  df = merge(df, dat2[,c(1,3,4,5)], by.x = 'Var2', by.y = 'Group.1')
  
  
  # 绘图
  p <- ggplot(df) +
    # 绘制柱状图
    geom_bar(aes(y, 
                 per2/2, 
                 fill = Var2),
             stat = 'identity', width = 1.2) +
    # 添加标签
    geom_text(aes(1, as.numeric(per.y2), 
                  label = label2),
              size =2.5, color = 'black') +
    # 绘制柱状图
    geom_bar(aes(x, per1, fill = Var1), 
             stat = 'identity', width = 1, color = 'white') +
    # 添加标签
    geom_text(aes(2, as.numeric(per.y1),label = label1),
              size = 2.5, color = 'black') +
    # 设置Y轴刻度
    scale_y_continuous(labels = scales::percent) +
    coord_polar(theta = "y") + # 转换坐标轴
    theme_void() +
    # scale_fill_igv() + # 设置填充色+=
    # scale_fill_manual(values = c("black","white","red","green","yellow","gray"))+
    scale_fill_manual(values = c("#b2acd1",colorspace::lighten("gray",0.3),colorspace::darken("gray",0.2),"#857fb6","#bd6a9a","#dc79ae"))+
    theme(legend.position = 'none') # 隐藏图例
}


if(T){
  typeIIAndIII <- fread(file="output/final_list/typeIIAndIII.txt")
  typeIIAndIII_n <- nrow(typeIIAndIII)
  typeII_n <- nrow(typeIIAndIII[proximal_type=="another protein product"])
  typeIII_n <- nrow(typeIIAndIII[proximal_type=="non_coding"])
  
  print("numbers of type II and type III genes")
  table(typeIIAndIII$proximal_type,typeIIAndIII$gene_type)
  
  typeI_n <- nrow(res[res$diff_expr=="NotDiff" & res$diff_pau=="Diff",])-typeIIAndIII_n
  
  res$cds[res$gene_name %in% typeIIAndIII$gene_name] <- "disrupt"
  res$cds[!res$gene_name %in% typeIIAndIII$gene_name] <- "not disrupt"
  
  res %>% mutate(diff_expr_1=case_when(diff_expr=="UP"~"Diff",diff_expr=="DOWN"~"Diff",diff_expr=="NotDiff"~"NotDiff")) ->res
  
  df1 <- data.frame(table(res$diff_expr_1))
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
  ggsave("/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/pie1.pdf",p,width = 5,height = 5)
  
  df2 <- data.frame(table(res[res$diff_expr_1=="NotDiff",]$diff_pau))
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
  ggsave("/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/pie2.pdf",p,width = 5,height = 5)
  
  
  df1 <- data.frame(Var1=c("Type I","Type II","Type III"),Freq=c(typeI_n,typeII_n,typeIII_n))
  df1$Var1 <- factor(df1$Var1,levels = c("Type I","Type II","Type III"))
  df1$per <- df1$Freq/sum(df1$Freq)
  labs <- paste0(df1$Freq," (",round(df1$per*100,2),"%",")")
  p <- df1 %>% 
    ggpie("Freq",
          fill = "Var1", color = "black",
          palette = c("#bd6a9a","#817bae","#b5b6b6"),
          lab.pos = "in", lab.font = c("bold","black"),
          label = labs)+
    theme(plot.title = element_text(hjust=0.5,vjust = -80))+
    labs(fill="Effect On Coding")
  ggsave("/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/pie3.pdf",p,width = 5,height = 5)
  
}
