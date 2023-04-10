
library("e1071")
library("data.table")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyfst")
library("EnvStats")

# read data in 
res_df <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw.txt")

# evo data
res_df %>% column_to_rownames("gene_name") ->res_df_1
range <- function(x){
  max(x)-min(x)
}
apply(res_df_1, 1, range) %>% data.frame() %>% rename(range=".") %>% rownames_to_column("gene_name") -> evo

# process devo data
devo=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

# merge
merge(evo,devo,"gene_name") %>% filter(diff_expr=="NotDiff") %>% mutate(delta_pau_abs=abs(delta_pau),range_abs=abs(range)) %>% 
  mutate(range_abs_scale=scale(range_abs),delta_pau_abs_scale=scale(delta_pau_abs))-> df
merge(evo,devo,"gene_name") %>% mutate(delta_pau_abs=abs(delta_pau),range_abs=abs(range))-> df_1


# skewness range distribution
ggplot(evo,aes(x=range))+
  geom_histogram()+
  theme_1()+
  xlab("Skewness range")+ylab("Count") -> p1
ggplot(evo,aes(x=range/100))+
  geom_histogram()+
  theme_1()+
  xlab("Skewness range/100")+ylab("Count") -> p2
ggplot(devo,aes(x=abs(delta_pau)))+
  geom_histogram()+
  theme_1()+
  xlab("Delta PAU")+ylab("Count") -> p3

p1+p2+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/output/pau_abs_fc_dis.pdf",width = 12,height = 5)

res <- function(df){
  df %>% mutate(bin=cut_interval(df$delta_pau_abs,5)) ->df_binned
  df %>% mutate(bin=cut_number(df$delta_pau_abs,5)) ->df_binned_cut_number
  compare <- function(p,df){
    le <- levels(df$bin)
    my_comparisons <- list(c(le[1],le[2]),
                           c(le[2],le[3]),
                           c(le[3],le[4]),
                           c(le[4],le[5]))
    # my_comparisons <- list(c(le[1],le[2]))
    res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"),label.y = c(130,135,140,145))
    res
  }
  
  compare_lim <- function(p,df){
    le <- levels(df$bin)
    my_comparisons <- list(c(le[1],le[2]),
                           c(le[2],le[3]),
                           c(le[3],le[4]),
                           c(le[4],le[5]))
    # my_comparisons <- list(c(le[1],le[2]))
    res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"),label.y = 20)
    res
  }
  plot_cus <- function(df,ylim){
    p <- ggplot(df,mapping = aes(x=bin,y=range_abs,fill=bin))+geom_boxplot(outlier.shape = NA,width=0.6)+theme_classic()+
      stat_boxplot(geom = "errorbar",width=0.3)+coord_trans(ylim=c(0,ylim))+labs(x="Delta PAU",y = "Skewness Range")+
      stat_summary(fun.data = function(x) {
        return(data.frame(y=ylim, label = length(x)))}, 
        geom = "text", hjust = 0.5, color = "black", size = 4)+
      theme(axis.text.x = element_text(size = 12,angle=35,vjust = 0.5,hjust = 0.5),
            axis.title = element_text(size = 14))+
      guides(fill = "none")+
      scale_fill_tableau()
    return(p)
  }
  plot_cus_nolim <- function(df,anno_y_pos){
    p <- ggplot(df,mapping = aes(x=bin,y=range_abs,fill=bin))+
      # geom_boxplot(outlier.shape = NA)+
      geom_violin(trim = T)+
      geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+
      # geom_jitter(size=1,color="gray40",alpha=0.3,position=position_jitter(width = 0.2,height = 0.2))+
      theme_classic()+
      # stat_boxplot(geom = "errorbar",width=0.5)+
      labs(x="Delta PAU",y = "Skewness Range")+
      stat_summary(fun.data = function(x) {
        return(data.frame(y=anno_y_pos, label = length(x)))}, 
        geom = "text", hjust = 0.5, color = "black", size = 4)+
      theme(axis.text.x = element_text(size = 12,angle=35,vjust = 0.5,hjust = 0.5),
            axis.title = element_text(size = 14))+
      guides(fill = "none")+
      coord_trans(ylim=c(0,anno_y_pos))
    return(p)
  }
  p1 <- compare_lim(plot_cus(df_binned,30),df_binned)
  p2 <- compare_lim(plot_cus(df_binned_cut_number,30),df_binned_cut_number)
  p3 <- compare(plot_cus_nolim(df_binned,160),df_binned)
  p4 <- compare(plot_cus_nolim(df_binned_cut_number,160),df_binned_cut_number)
  l <- list(p1,p2,p3,p4)
  return(l)
}
res(df)
res(df_1)

cor.test(df$delta_pau_abs,df$range_abs)
cor.test(df_1$delta_pau_abs,df_1$range_abs)

ggscatter(df,x="delta_pau_abs",y="range_abs",add = "reg.line",conf.int = TRUE,cor.coef = TRUE,cor.method = "pearson")+
  labs(x="Delta PAU",y = "Skewness Range")
ggscatter(df_1,x="delta_pau_abs",y="range_abs",add = "reg.line",conf.int = TRUE,cor.coef = TRUE,cor.method = "pearson")+
  labs(x="Delta PAU",y = "Skewness Range")
