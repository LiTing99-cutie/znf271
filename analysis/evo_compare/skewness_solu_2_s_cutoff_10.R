
library("e1071")
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

# read data in 
res_df <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo_raw_s_cutoff_10.txt")

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
ggplot(devo,aes(x=abs(delta_pau)))+
  geom_histogram()+
  theme_1()+
  xlab("Delta PAU")+ylab("Count") -> p3

p1+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/output/pau_abs_fc_dis_s_cutoff_10.pdf",width = 8,height = 5)

res <- function(df){
  df %>% mutate(bin=cut_number(df$delta_pau_abs,5)) ->df_binned_cut_number
  compare <- function(p,df){
    le <- levels(df$bin)
    my_comparisons <- list(c(le[1],le[2]),
                           c(le[2],le[3]),
                           c(le[3],le[4]),
                           c(le[4],le[5]))
    res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"))
    res
  }
  plot_cus <- function(df,anno_y_pos){
    p <- ggplot(df,mapping = aes(x=bin,y=range_abs,fill=bin))+
      geom_violin(trim = T)+
      geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+
      theme_classic()+
      labs(x="Delta PAU",y = "Skewness Range")+
      # labs(x=NULL,y = NULL)+
      stat_summary(fun.data = function(x) {
        return(data.frame(y=anno_y_pos, label = paste0("n= ",length(x)) ))}, 
        geom = "text", hjust = 0.5, color = "black", size = 5)+
      theme_1()+
      theme(axis.text.x = element_text(size = 15,angle=35,vjust = 0.5,hjust = 0.5))+
      guides(fill = "none")
    return(p)
  }
  
  p <- compare(plot_cus(df_binned_cut_number,65),df_binned_cut_number)
  return(p)
}
res(df_1) +ggtitle("PAU +-")-> p1
res(df) +ggtitle("Expression- & PAU +-") ->p2

cor.test(df$delta_pau_abs,df$range_abs)
cor.test(df_1$delta_pau_abs,df_1$range_abs)

ggscatter(df_1,x="delta_pau_abs",y="range_abs",add = "reg.line",conf.int = TRUE,cor.coef = TRUE,cor.method = "pearson")+
  labs(x="Delta PAU",y = "Skewness Range")+
  # labs(x=NULL,y = NULL)+
  theme_1() -> p3
ggscatter(df,x="delta_pau_abs",y="range_abs",add = "reg.line",conf.int = TRUE,cor.coef = TRUE,cor.method = "pearson")+
  labs(x="Delta PAU",y = "Skewness Range")+
  # labs(x=NULL,y = NULL)+
  theme_1() -> p4

p1+p2+p3+p4 -> p

ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/output/pau_evo_devo_cutoff_10.pdf",width = 10,height = 8)

