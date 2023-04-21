source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args <- commandArgs(T)
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_per.bed"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te.bed"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/exon.bed6"
  args[4] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}
# 1 te_per
# 2 te_collapse
# 3 all_exon
# 3 output path

#### read data in ####
do_te_per <- fread(args[1])
do_te <- fread(args[2])
all_exon <- fread(args[3])

#### intersect ####
bedtoolsr::bt.intersect(do_te,all_exon,wo=T,s=T) -> inter
# filter out transcript with dominant exons
!inter$V10 %in% do_te_per$transcript_id -> fil
inter[fil,] -> inter
# filter out overlap includes terminal exon backbone
inter$V8-inter$V2 >0 | inter$V9-inter$V3<0 -> fil
inter[fil,] -> inter

#### get clean dominant te ####
inter %>% distinct(V4) -> fil_gene
do_te[!do_te$gene_id %in% fil_gene$V4] -> do_te_fil
do_te_per[!do_te_per$gene_id %in% fil_gene$V4] -> do_te_per_fil

# extend end of do_te to 30 bp
rbind(filter(do_te_fil,strand=="+") %>% mutate(end=end+30),
      filter(do_te_fil,strand=="-") %>% mutate(start=start-30)) %>% mutate(start=case_when(start<0~0,start>=0~start)) -> do_te_fil_extend
# write
file_path=paste0(args[4],"/do_te_fil.bed")
fwrite_c(do_te_fil_extend,file_path)
file_path=paste0(args[4],"/do_te_per_fil.bed")
fwrite_c(do_te_per_fil,file_path)
