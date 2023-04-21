args=commandArgs(T)
# 1 input file
# 2 output path
source("/home/user/data2/lit/bin/lit_utils.R")
lib()
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/human/PAusage.fil.bed6+"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_fil.bed"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

# function
fwd_pro_dis_bin <- function(inter){
  # get proximal and distal end
  inter %>% group_by(V4) %>% summarise(min=min(V9),max=max(V9)) -> pro_dis_end
  # refine
  merge(te,pro_dis_end,by.x="gene_id",by.y = "V4") -> tmp
  data.frame(tmp$gene_id,tmp$chr,tmp$strand,tmp$start,tmp$min,tmp$min,tmp$max) -> res
  colnames(res) <- colnames(res) <- c("gene","chr","strand","pro_s","pro_e","dis_s","dis_e")
  return(res)
}
rvs_pro_dis_bin <- function(inter){
  # get proximal and distal end
  inter %>% group_by(V4) %>% summarise(min=min(V8),max=max(V8)) -> pro_dis_end
  # refine
  merge(te,pro_dis_end,by.x="gene_id",by.y = "V4") -> tmp
  data.frame(tmp$gene_id,tmp$chr,tmp$strand,tmp$max,tmp$end,tmp$min,tmp$max) -> res
  colnames(res) <- c("gene","chr","strand","pro_s","pro_e","dis_s","dis_e")
  return(res)
}

# file read in 
pas <- fread(args[1])
te <- fread(args[2])

# intersect
bedtoolsr::bt.intersect(te,pas,wa=T,wb=T,s=T) %>% filter(V4==V10)  -> inter
intersect(te$gene_id %>% unique(),pas$V4 %>% unique()) %>% length() ->n1
distinct(inter,V4) %>% nrow() ->n2

# pas in dominant te > 1
inter %>% count(V4) %>% filter(n>1) %>% .$V4 -> fil
inter[inter$V4 %in% fil,] -> inter_fil
distinct(inter_fil,V4) %>% nrow()->n3

# run
inter_fil %>% filter(V6=="+") -> fwd
inter_fil %>% filter(V6=="-") -> rvs
rbind(fwd_pro_dis_bin(fwd),rvs_pro_dis_bin(rvs)) -> pro_dis_bin
# filter on bin length
pro_dis_bin %>% mutate(pro_bin_len=pro_e-pro_s,dis_bin_len=dis_e-dis_s) %>% 
  filter(pro_bin_len>=100 & dis_bin_len>=100) %>% .$gene -> fil
pro_dis_bin[pro_dis_bin$gene %in% fil,] -> pro_dis_bin_fil
nrow(pro_dis_bin_fil) ->n4
inter_fil[inter_fil$V4 %in% fil,] -> inter_fil_fil

# write
filepath=paste0(args[3],"/pro_dis_bin.txt")
fwrite_c(pro_dis_bin_fil,filepath)
filepath=paste0(args[3],"/do_te_pas.txt")
fwrite_c(inter_fil_fil,filepath)
filepath=paste0(args[3],"/pas_map.sta.txt")
cat("Covered by iso-seq pas and identified to have dominant exons",n1,"\n",file = filepath,sep = '\t')
cat("Pas intersect with terminal exon backbone",n2,"\n",file = filepath,sep = '\t',append = T)
cat("Over one PAS after map",n3,"\n",file = filepath,sep = '\t',append = T)
cat("Cutoff on bin length",n4,"\n",file = filepath,sep = '\t',append = T)
