source("/home/user/data2/lit/bin/lit_utils.R")

# res
res <- readRDS(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

# miRNA
apa_miRNA_pro <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/AGO2/hg38.pro.miRNA.txt")
apa_miRNA_dis <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/AGO2/hg38.dis.miRNA.txt")

# function
cnt_m_s <- function(df){
  df %>% count(V4) ->df
  return(df)
}
group_mean <- function(df,col){
  col <- enquo(col)
  df %>% group_by(UQ(col)) %>% summarise(n_mean=mean(n)) -> t_1
  df %>% group_by(UQ(col)) %>% summarise(n_pro_mean=mean(n_pro)) -> t_2
  df %>% group_by(UQ(col)) %>% summarise(n_dis_mean=mean(n_dis)) -> t_3
  return(list(t_1,t_2,t_3))
}
plot_cus <- function(df,x,y,xlab,ylab,compare_group){
  p <- ggboxplot(df,x=x,y=y,fill=x,outlier.shape = NA,width = 0.6,lwd=1,
                 palette = c('#f8d396','#a3bdd8'))+
    geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0.2))+
    labs(x=xlab,y = ylab)+
    guides(fill = "none")+theme_1()
  compare <- function(p){
    my_comparisons <- compare_group
    res <- p+ stat_compare_means(comparisons = my_comparisons,
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"))
    res
  }
  compare(p)
}

# merge
merge(cnt_m_s(apa_miRNA_pro),cnt_m_s(apa_miRNA_dis),"V4",all = T) -> tmp
tmp[is.na(tmp)] <- 0                                                                             
colnames(tmp) <- c("gene_name","n_pro","n_dis")
tmp %<>% mutate(n=n_pro+n_dis)
merge(tmp,res,"gene_name",all.y=T)  -> tmp_1
tmp_1[is.na(tmp_1)] <- 0 

#### pau+ vs pau- ####

# compare mean
group_mean(tmp_1,diff_pau)

# plot
tmp_1 %<>% mutate(diff_pau_rn=if_else(diff_pau=="Diff","Y","N"))
plot_cus(tmp_1,"diff_pau_rn","n","Whether abs value of \n delta PAU >= 0.2","Number of miRNA-binding sites on terminal exon",list(c("Y","N"))) -> p1
plot_cus(tmp_1,"diff_pau_rn","n_pro","Whether abs value of \n delta PAU >= 0.2","Number of miRNA-binding sites on proximal region",list(c("Y","N"))) -> p2
plot_cus(tmp_1,"diff_pau_rn","n_dis","Whether abs value of \n delta PAU >= 0.2","Number of miRNA-binding sites on distal region",list(c("Y","N"))) -> p3
p1+p2+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/output/pau_miRNA_ago.pdf",height = 6,width = 12)

#### TypeII and Type III vs TypeI ####
typeII_III <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt")
tmp_1 %>% filter(diff_pau=="Diff"&diff_expr=="NotDiff") -> tmp_3
tmp_3$gene_name %in%  typeII_III$gene_name -> condi
tmp_3$type[condi] <- "TypeII and Type III"
tmp_3$type[!condi] <- "TypeI"
tmp_3 %<>% mutate(dis_mi_per =n_dis/(n+0.01))

# compare mean
group_mean(tmp_3,type)
# plot
plot_cus(tmp_3,"type","n","Type","Number of miRNA-binding sites on terminal exon",list(c("TypeI","TypeII and Type III"))) -> p1 
plot_cus(tmp_3,"type","n_pro","Type","Number of miRNA-binding sites on proximal region",list(c("TypeI","TypeII and Type III"))) -> p2 
plot_cus(tmp_3,"type","n_dis","Type","Number of miRNA-binding sites on distal region",list(c("TypeI","TypeII and Type III"))) -> p3 

p1+p2+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/output/type_miRNA_ago.pdf",height = 6,width = 12)

#### pau+ expr- vs pau+ expr+ ####
tmp_1 %>% mutate(type_list = case_when(
  diff_pau=="Diff" & diff_expr=="NotDiff" ~ "PAU+ Expr-",
  diff_pau=="Diff" & diff_expr!="NotDiff" ~ "PAU+ Expr+"
)) -> tmp_2
na.omit(tmp_2) -> tmp_3

# compare mean
group_mean(tmp_3,type_list)

plot_cus(tmp_3,"type_list","n",NULL,"Number of miRNA-binding sites on terminal exon",list(c("PAU+ Expr-","PAU+ Expr+"))) -> p1
plot_cus(tmp_3,"type_list","n_pro",NULL,"Number of miRNA-binding sites on proximal region",list(c("PAU+ Expr-","PAU+ Expr+"))) -> p2
plot_cus(tmp_3,"type_list","n_dis",NULL,"Number of miRNA-binding sites on distal region",list(c("PAU+ Expr-","PAU+ Expr+"))) -> p3

p1+p2+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/miRNA/output/pau_expr_miRNA_ago.pdf",height = 6,width = 12)
