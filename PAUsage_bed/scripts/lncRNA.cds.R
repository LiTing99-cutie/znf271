source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/predict/lncRNA_orf_prot_co_filter_genomic.txt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/map/transcript_id_chr_strand_symbol.txt"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

#### Read files in ####
lncRNA <- read.table(args[1],header = T)
map <- read.table(args[2],col.names = c("transcript_id","chr","strand","gene_symbol"))
head(lncRNA)
head(map)

#### Output ####
merge(lncRNA,map,by = "transcript_id") %>% filter(strand.x==strand.y) %>% 
  select(chr,start,end,transcript_id,gene_symbol,strand.y) -> lncRNA_cds_s_e

filepath=paste0(args[3],"/lncRNA_cds_s_e.txt")
fwrite_c(lncRNA_cds_s_e,filepath)
