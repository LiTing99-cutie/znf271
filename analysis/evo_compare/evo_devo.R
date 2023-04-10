
# packages
library(ggplot2)
library(ggpubr)

# process devo data
res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
res$diff_expr[res$diff_expr=="DOWN" | res$diff_expr=="UP"] <- "Diff"
res$type[res$diff_pau=="Diff" & res$diff_expr=="NotDiff"] <- "TypeI"
res$type[res$diff_pau=="NotDiff" & res$diff_expr=="NotDiff"] <- "TypeII"
res$type[res$diff_pau=="NotDiff" & res$diff_expr=="Diff"] <- "TypeIII"
res$type[res$diff_pau=="Diff" & res$diff_expr=="Diff"] <- "TypeIIII"

modi_res <- function(cutoff){
  res$type[abs(res$delta_pau)>=cutoff & res$diff_expr=="NotDiff"] <- "TypeI"
  res$type[abs(res$delta_pau)<cutoff & res$diff_expr=="NotDiff"] <- "TypeII"
  res$type[abs(res$delta_pau)<cutoff & res$diff_expr=="Diff"] <- "TypeIII"
  res$type[abs(res$delta_pau)>=cutoff & res$diff_expr=="Diff"] <- "TypeIIII"
  return(res)
}

modi_res_no_abs <- function(cutoff){
  res$type[res$delta_pau>=cutoff & res$diff_expr=="NotDiff"] <- "TypeI"
  res$type[res$delta_pau<cutoff & res$diff_expr=="NotDiff"] <- "TypeII"
  res$type[res$delta_pau<cutoff & res$diff_expr=="Diff"] <- "TypeIII"
  res$type[res$delta_pau>=cutoff & res$diff_expr=="Diff"] <- "TypeIIII"
  return(res)
}

# read evo data in 
# evo <- fread(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hrm_diff.txt")
evo_1 <- fread(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo.txt")
evo_2 <- fread(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_1_evo.txt")
evo_3 <- readRDS(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_2_evo.rds")
# merge
m_evo_devo <- function(cutoff,modi_res,evo){
  merge(modi_res(cutoff),evo,"gene_name") -> m
  return(m)
}


# compare

# choose which col to compare
compare <- function(m,col){
  table(m$type,m[,col]) %>% head(2) %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> pau_10_expr_0
  table(m$type,m[,col]) %>% .[c(4,3),] %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> pau_10_expr_1
  table(m$type,m[,col]) %>% .[c(3,2),] %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> expr_10_pau_0
  table(m$type,m[,col]) %>% .[c(4,1),] %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> expr_10_pau_1
  table(m$diff_pau,m[,col]) %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> pau_c
  table(m$diff_expr,m[,col]) %>% as.matrix() %>% 
    fisher.test(alternative="two.sided") -> expr_c
  
  get_estimate <- function(df){
    return(df$estimate)
  }
  
  get_p_value <- function(df){
    return(df$p.value)
  }
  
  l <- list(pau_10_expr_0,pau_10_expr_1,expr_10_pau_0,expr_10_pau_1,pau_c,expr_c)
  
  lapply(l, get_estimate) %>% unlist() -> estimate
  lapply(l, get_p_value) %>% unlist() -> p_value
  
  res <- data.frame(compare_group=c("pau+expr- VS pau-expr-","pau+expr+ VS pau-expr+","expr+pau- VS expr-pau-","expr+pau+ VS expr-pau+","pau+ VS pau-","expr+ VS expr-"),
                    odds_ratio=estimate,
                    p_value=p_value)
  return(res)
}

plot_odds_p <- function(df){
  df %<>% mutate(.,p_signif=case_when(p_value<0.001~"***",
                                      p_value>=0.001 & p_value<0.01~"**",
                                      p_value>=0.01 & p_value<0.05~"*",
                                      p_value>=0.05~"N.S."))
  
  x <- factor(df$compare_group,levels = c("pau+expr- VS pau-expr-","pau+expr+ VS pau-expr+","expr+pau- VS expr-pau-","expr+pau+ VS expr-pau+","pau+ VS pau-","expr+ VS expr-"))
  p <- ggplot(df,aes(x=x,y=odds_ratio))+geom_point(color="#e96140",size=2)+geom_hline(yintercept = 1,linetype=2)+
    geom_text(aes(label = p_signif),size=4,vjust=-0.6)+theme_pubr()+
    theme(axis.text.x = element_text(size = 16,angle=35,vjust = 0.5,hjust = 0.5),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16))+xlab(NULL)
  return(p)
}

compare(m_evo_devo(0.2,modi_res,evo_1),"diff") %>% plot_odds_p()
compare(m_evo_devo(0.3,modi_res,evo_1),"diff") %>% plot_odds_p()

compare(m_evo_devo(0.2,modi_res,evo_2),"diff_based_on_delta_skewness") %>% plot_odds_p()
compare(m_evo_devo(0.3,modi_res,evo_2),"diff_based_on_delta_skewness") %>% plot_odds_p()

compare(m_evo_devo(0.2,modi_res,evo_3[[1]]),"diff") %>% plot_odds_p()
compare(m_evo_devo(0.25,modi_res,evo_3[[1]]),"diff") %>% plot_odds_p()
compare(m_evo_devo(0.3,modi_res,evo_3[[1]]),"diff") %>% plot_odds_p()
# compare(m_evo_devo(0.2,modi_res),"h_s_or_not") %>% plot_odds_p()

# two evo list
merge(evo_1,evo_2,"gene_name") -> merge_12
table(merge_12$"diff",merge_12$"diff_based_on_delta_skewness")
