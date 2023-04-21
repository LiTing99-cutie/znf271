source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/pro_dis_bin.txt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/map/ensembl_gene_id_type_symbol.txt"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict"
}
ifelse(!dir.exists(args[3]),dir.create(args[3]),"Path already exist")
# extract gene name
fread(args[1]) %>% .$gene -> gene_name
map <- fread(args[2],header = F)
map[map$V3 %in% gene_name] %>% filter(V2!="protein_coding") %>% select(V1) -> lncRNA_gene_id

# fwrite
filepath <- paste0(args[3],"/to_predict_gene_id.txt")
fwrite(lncRNA_gene_id,filepath,sep='\t',col.names = F)
