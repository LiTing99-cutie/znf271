

setwd("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/typeII-III/pau_skewness")
source("/home/user/data2/lit/bin/lit_utils.R")

# function
compare <- function(p,group){
  my_comparisons <- group
  res <- p+ stat_compare_means(comparisons = my_comparisons,
                               method = "wilcox.test",method.args = list(alternative = "two.sided"))
  res
}

plot_cus <- function(df,x,y,xlab,ylab){
  p <- ggboxplot(df,x=x,y=y,fill=x,outlier.shape = NA,width = 0.6,lwd=1,
                 palette = c('#f8d396','#a3bdd8'))+
    geom_jitter(size=1,color="gray40",alpha=1,position=position_jitter(width = 0.2,height = 0.2))+
    labs(x=xlab,y = ylab)+
    guides(fill = "none")+theme_1()
  return(p)
}

# read data in 
res_df <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw_s.txt")

# evo data
res_df %>% column_to_rownames("gene_name") ->res_df_1
range <- function(x){
  max(x)-min(x)
}
apply(res_df_1, 1, range) %>% data.frame() %>% rename(range=".") %>% rownames_to_column("gene_name") -> evo

# devo data
devo=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
devo %<>% mutate(diff_expr=case_when(
  diff_expr=="UP" | diff_expr=="DOWN" ~ "Diff",
  TRUE ~ diff_expr
))
devo %<>% mutate(delta_pau_abs=abs(delta_pau))
# 1176 list
devo %>% mutate(type_list = case_when(
  diff_pau=="Diff" & diff_expr=="NotDiff" ~ "List",
  diff_pau=="Diff" & diff_expr!="NotDiff" ~ "NotList"
)) -> devo_1
na.omit(devo_1) -> devo_2
p <- plot_cus(devo_2,"type_list","delta_pau_abs","Type","Delta PAU")
compare(p,list(c("NotList","List"))) -> p1
# 61 list
typeII_III <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/typeIIAndIII.txt")
devo %>% filter(diff_pau=="Diff" & diff_expr=="NotDiff") -> devo_list
devo_list$type <- if_else(devo_list$gene_name %in% typeII_III$gene_name,"TypeII and Type III","Type I")
p <- plot_cus(devo_list,"type","delta_pau_abs","Type","Delta PAU")
compare(p,list(c("TypeII and Type III","Type I"))) ->p2

# expr- vs expr+
p <- plot_cus(devo,"diff_expr","delta_pau_abs","Type","Delta PAU")
compare(p,list(c("Diff","NotDiff"))) ->p3

# merge
merge(evo,devo_1,"gene_name") %>% na.omit(.) -> devo_1
devo_1$type_list <- factor(devo_1$type_list,levels = c("NotList","List"))
p <- plot_cus(devo_1,"type_list","range","Type","Skewness Range")
compare(p,list(c("NotList","List"))) -> p4

merge(evo,devo_list,"gene_name") -> tmp
p <- plot_cus(tmp,"type","range","Type","Skewness Range")
compare(p,list(c("TypeII and Type III","Type I"))) ->p5


# pau+expr- vs pau-expr-
devo %>% mutate(type_list = case_when(
  diff_pau=="Diff" & diff_expr=="NotDiff" ~ "List",
  diff_pau!="Diff" & diff_expr=="NotDiff" ~ "NotList"
)) -> devo_1
merge(evo,devo_1,"gene_name") %>% na.omit(.) -> devo_1
devo_1$type_list <- factor(devo_1$type_list,levels = c("NotList","List"))
p <- plot_cus(devo_1,"type_list","range","Type","Skewness Range")
compare(p,list(c("NotList","List"))) -> p6

# expr- vs expr+
merge(evo,devo,"gene_name") ->df
p <- plot_cus(df,"diff_expr","range","Type","Skewness Range")
compare(p,list(c("Diff","NotDiff"))) -> p7


#### output ####
# Type I vs Type II and III
p2+p5 -> p 
ggsave(p,filename = "output.pdf",width = 8,height = 5)
