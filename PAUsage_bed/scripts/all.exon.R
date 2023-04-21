
args <- commandArgs(T)
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}
# 1 gene_pred
# 2 output path

genepred <- fread(args[1])
unique(genepred$V12) -> gene_lst

main <- function(genepred,gene){
  each_gene <- genepred[genepred$V12==gene] 
  lapply(each_gene$V1, each_transcript_exon_bed,each_gene=each_gene) -> l
  names(l) <- each_gene$V1
  return(l)
}

each_transcript_exon_bed <- function(each_gene,transcript){
  each_gene[each_gene$V1==transcript] -> each_transcript
  each_transcript$V9 %>% str_split(',') %>% unlist() %>% as.numeric() -> start
  each_transcript$V10 %>% str_split(',') %>% unlist() %>% as.numeric() -> end
  exon_bed <- data.frame(chr=each_transcript$V2,start=start,end=end,puppet_1=".",puppet_2=".",strand=each_transcript$V3) %>% na.omit()
  return(exon_bed)
}

lapply(gene_lst, main,genepred=genepred) -> all_exon
names(all_exon) <- gene_lst

# write
file_path=paste0(args[2],"/all_exon.rds")
saveRDS(all_exon,file = file_path)





