source("/home/user/data2/lit/bin/lit_utils.R")
lib()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/stringtie/stringtie.rpkm.txt"
  args[2] <- "/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/metadata/human.metadata.clean.txt"
  args[3] <- "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt"
  args[4] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

#### Read rpkm ####
rpkm <- fread(args[1])
colnames(rpkm) <- c("type","gene_id","rpkm","sample")
index <- grep("proximal",rpkm$type)
rpkm_pro <- rpkm[index,]
index <- grep("distal",rpkm$type)
rpkm_dis <- rpkm[index,]

#### Read metadata ####
metadata=fread(args[2])
frag_score=fread(args[3],sep='\t')
colnames(frag_score) = c("sample","frag_score")
#### Merge ####
merge(metadata,frag_score) %>% filter(frag_score>0.885) -> md
# see if fragmentation score of embryo or postnatal groups have significant difference
md %>% group_by(stage) %>% summarise(mean=mean(frag_score),median=median(frag_score),n=n()) -> frag_score_group
compare_means(frag_score~stage,md)
ggboxplot(md,x="stage",y="frag_score",fill="stage",outlier.shape = NA,width = 0.6,lwd=1,
          palette = c('#f8d396','#a3bdd8'))+stat_compare_means()
## merge
merge(rpkm_pro,rpkm_dis,by=c("gene_id","sample")) %>% merge(md,"sample") %>%
  select(gene_id,rpkm.x,rpkm.y,stage) %>% filter(rpkm.x>=1) %>% mutate(pau=rpkm.y/rpkm.x)-> res
#### Prepare gene list ####
unique(res$gene_id) -> gene_lst

GetPvalue <- function(gene){
  res %>% filter(gene_id==gene) ->input
  input_b <- filter(input,stage=="embryo")
  input_a <- filter(input,stage=="postnatal")
  if(nrow(input_a)>=8 && nrow(input_b)>=8){
    p_1 <- wilcox.test(input_b$rpkm.x,input_a$rpkm.x)$p.value
  }else{
    p_1 <- NA
  }
  pau_b=median(input_b$pau)
  pau_a=median(input_a$pau)
  delta_pau <- median(input_a$pau)-median(input_b$pau)
  expr_b=mean(input_b$rpkm.x)
  expr_a=mean(input_a$rpkm.x)
  fc <- log2(expr_a/expr_b)
  out <- list("pau_b"=pau_b,"pau_a"=pau_a,"delta_pau"=delta_pau,"expr_b"=expr_b,"expr_a"=expr_a,"fc"=fc,"p_expr"=p_1)
  return(out)
}

sapply(gene_lst,GetPvalue) %>% apply(MARGIN = 2,unlist) %>% t() %>% data.frame()%>%
  rownames_to_column(.,"gene_name") -> res
res %>% na.omit(.) -> res_f

filepath=paste0(args[4],"/wilcox.txt")
fwrite_c(res_f,filepath)
filepath=paste0(args[4],"/frag_score_group.txt")
fwrite_c(frag_score_group,filepath)
