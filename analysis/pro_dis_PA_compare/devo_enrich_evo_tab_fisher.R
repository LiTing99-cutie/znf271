
# evo list read in
cp_pas_n <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n.txt")
cp_pas_n_diff <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n_diff.txt")
cp <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp.txt")
union_hrm <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/union_hrm.txt")

all <- cp_pas_n[,"gene_name"]
res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
all <- merge(all,res,"gene_name")

# evo list
rbind(filter(cp,p_adjust<0.05) %>% select("gene_name"),select(cp_pas_n_diff,"gene_name")) -> evo_diff_1
select(cp_pas_n_diff,"gene_name") -> evo_diff_2
filter(cp,p_adjust<0.05) %>% select("gene_name") -> evo_diff_3
union_hrm -> evo_diff_4
rbind(union_hrm,select(cp_pas_n_diff,"gene_name")) -> evo_diff_5
# fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/evo_diff_pro_dis_region.lst") ->evo_diff_6
# l <- list(evo_diff_1,evo_diff_2,evo_diff_3,evo_diff_4,evo_diff_5,evo_diff_6)
l <- list(evo_diff_1,evo_diff_2,evo_diff_3,evo_diff_4,evo_diff_5)

all$evo <- if_else(all$gene_name %in% evo_diff$gene_name,"Diff","NotDiff")


all$type[all$diff_pau=="Diff" & all$diff_expr=="NotDiff"] <- "TypeI"
all$type[all$diff_pau=="NotDiff" & all$diff_expr=="NotDiff"] <- "TypeII"

all$type[all$diff_pau=="NotDiff" & all$diff_expr!="NotDiff"] <- "TypeIII"
all$type[all$diff_pau=="Diff" & all$diff_expr!="NotDiff"] <- "TypeIIII"

# delta pau <0.2, expr sig change or not
# "TypeIII" + "TypeII"
# expr no sig change,delta pau >0.2 or not
# "TypeI" + "TypeII"

evo_devo <- function(evo_diff){
  all$evo <- if_else(all$gene_name %in% evo_diff$gene_name,"Diff","NotDiff")
  test <- function(condition){
    all$type <- if_else(condition,"Type","NotType")
    all$type <- factor(all$type,levels = c("Type","NotType"))
    table(all$type,all$evo) %>% as.matrix() -> matri
    matri %>% fisher.test(alternative="two.sided") ->fisher
    l <- data.frame(p_value=as.numeric(fisher$p.value),odds_ratio=as.numeric(fisher$estimate))
    return(l)
  }
  rbind(test(all$diff_pau=="Diff" & all$diff_expr=="NotDiff"),
        test(all$diff_pau=="Diff" & all$diff_expr!="NotDiff"),
        test(all$diff_pau=="NotDiff" & all$diff_expr=="NotDiff"),
        test(all$diff_pau=="NotDiff" & all$diff_expr!="NotDiff"),
        test(all$diff_pau=="Diff"),
        test(all$diff_expr!="NotDiff")) -> type_all
  
  type_all %<>% mutate(.,p_signif=case_when(p_value<0.001~"***",
                                            p_value>=0.001 & p_value<0.01~"**",
                                            p_value>=0.01 & p_value<0.05~"*",
                                            p_value>=0.05~"N.S."))
  
  x <- factor(c("PAU_NOT_EXPR","PAU_EXPR","NOT_PAU_NOT_EXPR","NOT_PAU_EXPR","PAU","EXPR"),
              levels = c("PAU_NOT_EXPR","PAU_EXPR","NOT_PAU_NOT_EXPR","NOT_PAU_EXPR","PAU","EXPR"))
  p <- ggplot(type_all,aes(x=x,y=odds_ratio))+geom_point(color="#e96140")+geom_hline(yintercept = 1,linetype=2)+
    geom_text(aes(label = p_signif),size=3,vjust=-0.5)+theme_pubr()+
    theme(axis.text.x = element_text(size = 8,angle=35,vjust = 0.5,hjust = 0.5))+xlab(NULL)
  # filename <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/type_enrich/",prefix,".pdf")
  # ggsave(p,filename = filename,height = 4,
  #        width = 4) 
  return(p)
}
 
lapply(l,evo_devo) -> res

res[[1]]+res[[2]]+res[[3]]+res[[4]]+res[[5]]+res[[6]] -> p
filename <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/type_enrich/","all",".png")
ggsave(p,filename = filename,height = 8,width = 8,dpi=300)


       