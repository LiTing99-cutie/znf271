source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
library(bedtoolsr)
args <- commandArgs(T)
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_per.bed"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te.bed"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/all_exon.rds"
  args[4] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}
# 1 te_per
# 2 te_collapse
# 3 all_exon
# 3 output path

do_te_per <- fread(args[1])
do_te <- fread(args[2])
all_exon <- readRDS(args[3])

# function
extract_bed_from_list <- function(gene,transcript){
  all_exon[[gene]][[transcript]] -> bed
  return(bed)
}
do_te_inter_other_trans <- function(gene){
  do_te[do_te$gene_id==gene] -> each_do_te
  do_te_per[do_te_per$gene_id==gene] %>% .$"transcript_id" -> transcript_id
  ! all_exon[[gene]] %>% names() %in% transcript_id -> inx
  names(all_exon[[gene]])[inx] -> other_trans
  # combine other transcript te
  other_trans_bed <- data.frame()
  for (i in other_trans){
    other_trans_bed <- rbind(other_trans_bed,extract_bed_from_list(gene,i))
  }
  bedtoolsr::bt.intersect(each_do_te,other_trans_bed,wo=T,s=T) -> inter
  return(data.frame(inter))
}

# test
do_te_inter_other_trans("ZNF271P")

# run
# gene_lst <- head(do_te$gene_id,100)
gene_lst <- do_te$gene_id

res <- data.frame()
i=1
system.time(for (gene in gene_lst){
  res <- rbind(res,do_te_inter_other_trans(gene))
  if (i %% 1000 ==0){print(i)}
  i=i+1
})

# write
file_path=paste0(args[4],"/do_te_inter_other_trans.bed")
fwrite_c(res,file_path)
