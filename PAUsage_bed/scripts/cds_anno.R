source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_per_fil.bed"
  args[3] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/pro_dis_bin.txt"
  args[4] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_pas.txt"
  args[5] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/do_te_fil.bed"
  args[6] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output/"
}

# function
type_based_on_annotation <- function(df){
  cds_s_e_fil_add_type[cds_s_e_fil_add_type$V12 %in% df$V4] %>% filter(type=="Coding but upstream of te") %>% .$V12 -> fil
  mutate(df,type=case_when(V4 %in% fil ~ "Coding but likely NMD",
                            TRUE ~ "Non-coding")) -> df
  df$general_type <- "Disrupt cds in terminal exon"
  return(df)

}
cnt_usage_rename_fwd <- function(df){
  do_te_pas[do_te_pas$V4 %in% df$V4,] %>% group_dt(V4,filter_dt(V9==min(V9))) %>% select(7:13) -> pro_cnt_usage
  merge(df,pro_cnt_usage,by.x = "V4",by.y = "V10") -> final
  colnames(final) <- c("gene_name","n_cds_end_pro","n_cds_end_dis","type","general_type","chr","start","end","cnt","strand","usage")
  return(final)
}
cnt_usage_rename_rvs <- function(df){
  do_te_pas[do_te_pas$V4 %in% df$V4,] %>% group_dt(V4,filter_dt(V8==max(V8))) %>% select(7:13) -> pro_cnt_usage
  merge(df,pro_cnt_usage,by.x = "V4",by.y = "V10") -> final
  colnames(final) <- c("gene_name","n_cds_end_pro","n_cds_end_dis","type","general_type","chr","start","end","cnt","strand","usage")
  return(final)
}
fwd <- function(strand){
  pro_dis_bin %>% filter(strand==strand) %>% select(chr,pro_s,pro_e,gene) -> pro_bed
  pro_bed$puppet = "."
  pro_bed$strand = strand
  pro_dis_bin %>% filter(strand==strand) %>% select(chr,pro_s,dis_e,gene) -> dis_bed
  dis_bed$puppet = "."
  dis_bed$strand = strand
  
  bedtoolsr::bt.intersect(pro_bed,cds_s_e_fil_coding,wa=T,wb=T,s=T) %>% 
    filter(V4==V11) %>% 
    filter(V3>=V9) %>% 
    distinct(V4,V9) %>% 
    group_by(V4) %>% 
    summarise(n=n()) -> pro
  
  bedtoolsr::bt.intersect(dis_bed,cds_s_e_fil_coding,wa=T,wb=T,s=T) %>% 
    filter(V4==V11) %>% 
    filter(V3>=V9) %>% 
    distinct(V4,V9) %>% 
    group_by(V4) %>% 
    summarise(n=n()) -> dis
  
  merge(pro,dis,"V4",all.y = T) %>% replace(is.na(.), 0)-> tmp
  
  filter(tmp,n.x==0) -> res
  
  res <- type_based_on_annotation(res)
  final <- cnt_usage_rename_fwd(res)
  
  return(final)
}
rvs <- function(strand){
  pro_dis_bin %>% filter(strand==strand) %>% select(chr,pro_s,pro_e,gene) -> pro_bed
  pro_bed$puppet = "."
  pro_bed$strand = strand
  pro_dis_bin %>% filter(strand==strand) %>% select(chr,dis_s,pro_e,gene) -> dis_bed
  dis_bed$puppet = "."
  dis_bed$strand = strand
  bedtoolsr::bt.intersect(pro_bed,cds_s_e_fil_coding,wa=T,wb=T,s=T) %>% 
    filter(V4==V11) %>% 
    filter(V2<=V8) %>% 
    distinct(V4,V8) %>% 
    group_by(V4) %>% 
    summarise(n=n()) -> pro
  
  bedtoolsr::bt.intersect(dis_bed,cds_s_e_fil_coding,wa=T,wb=T,s=T) %>% 
    filter(V4==V11) %>% 
    filter(V2<=V8) %>% 
    distinct(V4,V8) %>% 
    group_by(V4) %>% 
    summarise(n=n()) -> dis
  
  merge(pro,dis,"V4",all.y = T) %>% replace(is.na(.), 0)-> tmp
  
  filter(tmp,n.x==0) -> res
  
  res <- type_based_on_annotation(res)
  final <- cnt_usage_rename_rvs(res)
  
  return(final)
}

# read in
gpe <- fread(args[1])
transcript_fil <- fread(args[2])
pro_dis_bin <- fread(args[3])
do_te_pas <- fread(args[4])
do_te <- fread(args[5])

# filter cds annotation
gpe %>% select(V2,V6,V7,V1,V12,V3) -> cds_s_e
cds_s_e$V1 %in% transcript_fil$transcript_id -> fil

# [optional] classify transcripts with dominant terminal exons into coding, coding upstream of te and no-coding based on public annotation
cds_s_e[fil] -> cds_s_e_fil
merge(cds_s_e_fil,do_te,by.x="V12",by.y="gene_id") -> tmp
tmp %>% filter(V3=="+") %>% mutate(type=case_when(V6==V7 ~ "Non-coding",
                                                  V6!=V7 & V7>=start ~ "Coding",
                                                  V6!=V7 & V7<start ~ "Coding but upstream of te")) -> tmp_fwd
tmp %>% filter(V3=="-") %>% mutate(type=case_when(V6==V7 ~ "Non-coding",
                                                  V6!=V7 & V6<=end ~ "Coding",
                                                  V6!=V7 & V6>end ~ "Coding but upstream of te")) -> tmp_rvs

rbind(tmp_fwd,tmp_rvs) ->cds_s_e_fil_add_type

# cds annotation of transcripts with dominant terminal exons
cds_s_e_fil %>% filter(V6!=V7) -> cds_s_e_fil_coding

# run
rbind(fwd("+"),rvs("-")) -> res

# write
filepath <- paste0(args[6],"/disrupt_cds_pro_pa.txt")
fwrite_c(res,filepath)


