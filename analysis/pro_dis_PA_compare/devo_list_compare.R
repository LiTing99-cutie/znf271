suppressMessages(library(clusterProfiler))
suppressMessages(library(DOSE))

# we want to know what GO terms our evo list enrich

GO_func <- function(enrich,file_name){
  GO <- enrichGO(gene=enrich$gene_name,
                 OrgDb="org.Hs.eg.db",
                 ont="ALL",
                 keyType = "SYMBOL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)
  path <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/",file_name,".pdf")
  path_1 <- paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/",file_name,".txt")
  pdf(file = path,width = 6,height = 5)
  print(dotplot(GO))
  dev.off()
  write.table(GO,file = path_1,sep = '\t',quote = F,row.names = F)
}

# IS
cp_pas_n_diff <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n_diff.txt")
cp <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp.txt")
rbind(filter(cp,p_adjust<0.05) %>% select("gene_name"),select(cp_pas_n_diff,"gene_name")) -> enrich

# Top PAS
res_hm <- readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare/res_hm.rds")
res_hr <- readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare/res_hr.rds")
res_hrm <- readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare/res_hrm.rds")
enrich_hr <- filter(res_hr,p_adjust<0.05) %>% select(GN_H)
enrich_hm <- filter(res_hm,p_adjust<0.05) %>% select(GN_H)
enrich_1 <- filter(res_hrm,p_adjust<0.05) %>% select(GN_H) %>% rename(gene_name=GN_H)

# intersect
intersect(enrich_1$gene_name,enrich$gene_name) %>% data.frame() %>% rename(gene_name='.')-> enrich_2

GO_func(enrich,"IS")
GO_func(enrich_1,"Top_PAS")
GO_func(enrich_2,"Inter")