library("patchwork")


res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

fread("/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/terminal_exon_annotation.ref_based_all_1_loose.txt") %>% 
  select(V1) %>% rename(gene_name=V1) -> type_0
res %>% select(gene_name) -> type_1
res %>% filter(diff_pau=="Diff" & diff_expr=="NotDiff") %>% select(gene_name) -> type_2
res %>% filter(diff_pau=="Diff" & diff_expr!="NotDiff") %>% select(gene_name) -> type_3
res %>% filter(diff_pau=="NotDiff" & diff_expr=="NotDiff") %>% select(gene_name) -> type_4
res %>% filter(diff_pau=="NotDiff" & diff_expr!="NotDiff") %>% select(gene_name) -> type_5
res %>% filter(diff_pau=="Diff") %>% select(gene_name) -> type_6
res %>% filter(diff_expr!="NotDiff") %>% select(gene_name) -> type_7

cp_pas_n_diff <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n_diff.txt")
cp <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp.txt")

rbind(filter(cp,p_adjust<0.05) %>% select("gene_name"),select(cp_pas_n_diff,"gene_name")) -> enrich
# select(cp_pas_n_diff,"gene_name") -> enrich
# filter(cp,p_adjust<0.05) %>% select("gene_name") ->enrich

res_hrm <- readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare/res_hrm.rds")
enrich_1 <- filter(res_hrm,p_adjust<0.05) %>% select(GN_H) %>% rename(gene_name=GN_H)

intersect(enrich_1$gene_name,enrich$gene_name) %>% data.frame() %>% rename(gene_name='.')-> enrich_2

test <- function(type_1,type_2,enrich,alternative){
  intersect(enrich$gene_name,type_1$gene_name) -> enrich_type_1
  enrich_type_1 %>% length()-> n_1
  intersect(type_2$gene_name,enrich_type_1) %>% length()-> n_2
  x_1=n_2
  y_1=n_1-n_2
  x_2=nrow(type_2)-n_2
  y_2=nrow(type_1)-nrow(type_2)-(n_1-n_2)
  diff=(x_1/x_2)/(y_1/y_2)
  fisher.test(matrix(c(x_1,y_1,x_2,y_2),nrow=2),alternative = alternative) %>% .$p.value->p_value
  data.frame(x_1,y_1,x_2,y_2,p_value,diff) -> df
  return(df)
}

l <- list(type_2,type_3,type_4,type_5,type_6,type_7)

sapply(l, test,type_1=type_0,enrich=enrich,alternative="two.sided") -> diff_back_fisher_0
sapply(l, test,type_1=type_1,enrich=enrich,alternative="two.sided") -> diff_back_fisher_1

sapply(l, test,type_1=type_0,enrich=enrich_2,alternative="two.sided") -> diff_back_fisher_0_1
sapply(l, test,type_1=type_1,enrich=enrich_2,alternative="two.sided") -> diff_back_fisher_1_1


######

permutation <- function(type_1,type_2,enrich){
  df <- data.frame()
  for (i in 1:10000){
    sample(type_1$gene_name,nrow(type_2),replace = F) -> tmp_1
    intersect(tmp_1,enrich$gene_name) %>% length()-> n_1
    df <- rbind(df,n_1)
  }
  # permutation for 10000 times
  df -> tmp
  colnames(tmp) <- "Overlap"
  # overlap in fact
  intersect(type_2$gene_name,enrich$gene_name) %>% length() ->n
  # p_value
  p_value <- sum(tmp$Overlap>n)/10000
  p <- ggplot(tmp,aes(x="Type",y=Overlap))+geom_violin(fill="#f19a94",color="#e93224",width=1)+
    geom_boxplot(outlier.shape = NA,fill=NA,color="black",width=0.5)+
    geom_point(x="Type",y=n)+xlab(NULL)
  l=list(p,p_value)
  return(l)
}

diff_back <- function(type_1,enrich){
  lapply(l, permutation,type_1=type_1,enrich=enrich) ->per_re
  
  df <- data.frame()
  for (i in 1:6){
    per_re[[i]][[2]]->p_value
    df <- rbind(df,p_value)
  }
  colnames(df) <- "p_value"
  
  for (i in 1:6){
    name <- paste0("p_",i)
    per_re[[i]][[1]] -> p
    assign(name,p)
  }
  
  p_1+p_2+p_3+p_4+p_5+p_6 ->p
  l <- list(p,df)
  
  return(l)
  
}

diff_back(type_0,enrich) -> diff_back_mc_0
diff_back(type_1,enrich) -> diff_back_mc_1


