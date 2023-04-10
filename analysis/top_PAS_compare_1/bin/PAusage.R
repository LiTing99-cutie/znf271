args <- commandArgs(T)
source("/home/user/data2/lit/bin/lit_utils.R")
lib()
# file input
# file output

df <- fread(args[1])

df %>% separate(V4,sep = ":",into = c("gene_name","pos")) %>% group_dt(gene_name,mutate_dt(per=V5/sum(V5)))  -> tmp
tmp %>% select(V1,V2,V3,gene_name,pos,per,V6) %>% unite("gene_name_pos",gene_name:pos,sep = ":") -> res

fwrite(res,file = args[2],sep = '\t',col.names = F)
