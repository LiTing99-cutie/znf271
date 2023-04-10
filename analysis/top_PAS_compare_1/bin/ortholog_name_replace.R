args <- commandArgs(T)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

# file of species to be renamed

df <- fread(args[1]) %>% separate(V4,sep = ":",into = c("gene_name","pos")) 
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hr.txt",header = F)
ortho_add <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/lncRNA_ortholog/hr.lnc.cl.txt",header = F)
ortho <- rbind(ortho,ortho_add)

merge(df,ortho,by.x = "gene_name",by.y = "V2") %>% mutate(gene_name=NULL) %>% rename(gene_name=V1.y) %>% 
  unite("gene_name_pos",gene_name,pos,sep = ":") -> res


fwrite(res,file = args[2],sep = '\t',col.names = F) 