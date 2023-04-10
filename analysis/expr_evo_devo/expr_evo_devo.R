
# 1. evo
evo <- readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/expr_evo_devo/output/anova.rds")
# filter lowly expressed genes in three species
# filter(evo,!(h_mean<1 & r_mean<1 & m_mean<1)) -> evo
filter(evo,h_mean>=1 & r_mean>=1 & m_mean>=1) -> evo
# p_value adjust
evo$adjust_p_value <- p.adjust(evo$p_value) 
# sig gene number
table(evo$p_value<0.05)
table(evo$adjust_p_value<0.05)
# log2 fc
select(evo,h_mean,r_mean,m_mean) ->tmp
max_min <- function(x){
  log2(max(x)/min(x))
}
apply(tmp,1,max_min) %>% as.data.frame %>% rename(fc=".") -> tmp_1
cbind(evo,tmp_1) -> tmp_2
filter(tmp_2,adjust_p_value<0.05 & abs(fc)>1) %>% nrow()

# 2. devo
devo <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/expr_evo_devo/output/devo.txt")
# filter lowly expressed genes in three species
# filter(devo,!(expr_embryo <1 & expr_postnatal<1)) -> devo
filter(devo,expr_embryo >=1 & expr_postnatal>=1) -> devo
# p_value adjust
devo$adjust_p_value <- p.adjust(devo$p_expr) 
# sig gene number
table(devo$p_expr<0.05)
table(devo$adjust_p_value<0.05)
# log2 fc
devo$fc <- log2(devo$expr_postnatal/devo$expr_embryo)
filter(devo,adjust_p_value<0.05 & abs(fc)>1) %>% nrow()

# 2.5 devo and evo distribution
# skewness range distribution
ggplot(tmp_2,aes(x=fc))+
  geom_histogram()+
  theme_1()+
  xlab("Absolute value of log2 (Fold Change) during evolution")+ylab("Count") -> p1
ggplot(devo,aes(x=abs(fc)))+
  geom_histogram()+
  theme_1()+
  xlab("Absolute value of log2 (Fold Change) during development")+ylab("Count") -> p2
p1+p2 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/output/expr_abs_fc_dis.pdf",width = 15,height = 3)

# 3. merge and test
tmp_2$diff <- if_else(tmp_2$adjust_p_value<0.05 & abs(tmp_2$fc)>1,"Diff","NotDiff")
devo$diff <- if_else(devo$adjust_p_value<0.05 & abs(devo$fc)>1,"Diff","NotDiff")
merge(tmp_2,devo,"gene_name") ->m
table(m$diff.x,m$diff.y) %>% fisher.test()

standardize <- function(x){
  (x-mean(x))/sd(x)
}
m %>% mutate(diff_devo_abs=abs(fc.y),diff_evo_abs=abs(fc.x)) %>% mutate(diff_devo_abs_scale=scale(diff_devo_abs),
                                                                        diff_evo_abs_scale=scale(diff_evo_abs)) -> df


res <- function(df){
  df %>% mutate(bin=cut_interval(df$diff_devo_abs,5)) ->df_binned
  df %>% mutate(bin=cut_number(df$diff_devo_abs,5)) ->df_binned_cut_number
  compare <- function(p,df){
    le <- levels(df$bin)
    my_comparisons <- list(c(le[1],le[2]),
                           c(le[2],le[3]),
                           c(le[3],le[4]),
                           c(le[4],le[5]))
    # my_comparisons <- list(c(le[1],le[2]))
    res <- p+ stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                                 method = "wilcox.test",method.args = list(alternative = "two.sided"))
    res
  }
  plot_cus <- function(df,ylim){
    p <- ggplot(df,mapping = aes(x=bin,y=diff_evo_abs,fill=bin))+geom_boxplot(outlier.shape = NA)+theme_classic()+
      stat_boxplot(geom = "errorbar",width=0.5)+coord_trans(ylim=c(0,ylim))+labs(x="Devo-log2 (Fold Change)",y = "Evo-log2 (Fold Change)")+
      stat_summary(fun.data = function(x) {
        return(data.frame(y=ylim, label = length(x)))}, 
        geom = "text", hjust = 0.5, color = "black", size = 4)+
      theme(axis.text.x = element_text(size = 12,angle=35,vjust = 0.5,hjust = 0.5),
            axis.title = element_text(size = 14))+
      guides(fill = "none")
    return(p)
  }
  plot_cus_nolim <- function(df,anno_y_pos){
    p <- ggplot(df,mapping = aes(x=bin,y=diff_evo_abs,fill=bin))+
      geom_violin(trim = T)+
      geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+
      # geom_jitter(size=1,color="gray40",alpha=0.3,position=position_jitter(width = 0.2,height = 0.2))+
      theme_classic()+
      # stat_boxplot(geom = "errorbar",width=0.5)+
      labs(x="Absolute value of log2 (Fold Change) during development",y = "Absolute value of log2 (Fold Change) during evolution")+
      stat_summary(fun.data = function(x) {
        return(data.frame(y=anno_y_pos, label = length(x)))}, 
        geom = "text", hjust = 0.5, color = "black", size = 4)+
      theme(axis.text.x = element_text(size = 12,angle=35,vjust = 0.5,hjust = 0.5),
            axis.title = element_text(size = 14))+
      guides(fill = "none")+
      coord_trans(ylim=c(0,anno_y_pos))+
      scale_fill_tableau()
      # scale_fill_aaas()
    return(p)
  }
  p1 <- plot_cus(df_binned,6)
  p2 <- plot_cus(df_binned_cut_number,4)
  p3 <- compare(plot_cus_nolim(df_binned,11),df_binned)
  p4 <- compare(plot_cus_nolim(df_binned_cut_number,11),df_binned_cut_number)
  l <- list(p1,p2,p3,p4)
  return(l)
}

res(df)

cor.test(df$diff_devo_abs,df$diff_evo_abs)

ggscatter(df,x="diff_devo_abs",y="diff_evo_abs",add = "reg.line",conf.int = TRUE,cor.coef = TRUE,cor.method = "pearson")+
  labs(x="Absolute value of log2 (Fold Change) during development",y = "Absolute value of log2 (Fold Change) during evolution")
