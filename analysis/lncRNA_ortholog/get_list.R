

res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
map <- fread("/home/user/data2/lit/project/ZNF271/02-APA/annotation/map/ensembl_gene_id_type_symbol.txt",header = F)

merge(res[1],map,by.x="gene_name",by.y="V3") %>% filter(V2!="protein_coding")-> gn_gi_gt

fwrite(gn_gi_gt["gene_name"],"/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/gene_lst.txt",col.names =F)
