args <- commandArgs(T)
# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/terminal_exon/gencode.v41.basic.annotation.terminal_exon.genePredExt"
  args[2] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output"
}
# 1 terminal exon
# 2 output path
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

# function
dominant_te_fwd <- function(df){
  df %>% 
    unite("gene_id_pos",gene_id,start,remove = F,sep = ":") -> df_add_unique_id
  df_add_unique_id %>% count(gene_id_pos) %>% 
    separate(gene_id_pos,into = c("gene_id","pos"),sep = ":") %>% group_dt(gene_id,filter_dt(n==max(n))) %>%
    rename(n_isoform=n) %>% 
    add_count_dt(gene_id) %>% filter(n==1) %>% unite("gene_id_pos",gene_id,pos,remove = F,sep = ":") -> dominant_id
  merge(dominant_id[,"gene_id_pos"],df_add_unique_id,"gene_id_pos") -> res
  res %>% group_dt(gene_id,filter_dt(end==max(end))) %>% distinct(gene_id,end,.keep_all = T) -> te_backbone
  return(list(res,te_backbone))
}
dominant_te_rvs <- function(df){
  df %>% 
    unite("gene_id_pos",gene_id,end,remove = F,sep = ":") -> df_add_unique_id
  df_add_unique_id %>% count(gene_id_pos) %>% 
    separate(gene_id_pos,into = c("gene_id","pos"),sep = ":") %>% group_dt(gene_id,filter_dt(n==max(n))) %>%
    rename(n_isoform=n) %>% 
    add_count_dt(gene_id) %>% filter(n==1) %>% unite("gene_id_pos",gene_id,pos,remove = F,sep = ":") -> dominant_id
  merge(dominant_id[,"gene_id_pos"],df_add_unique_id,"gene_id_pos") -> res
  res %>% group_dt(gene_id,filter_dt(start==min(start))) %>% distinct(gene_id,start,.keep_all = T) -> te_backbone
  return(list(res,te_backbone))
}

# te 
te <- fread(args[1])
colnames(te) <- c("chr","strand","start","end","gene_id","transcript_id")
filter(te,strand=="+") -> fwd
filter(te,strand=="-") -> rvs

# run
dominant_te_fwd(fwd) -> fwd_do_te
dominant_te_rvs(rvs) -> rvs_do_te
rbind(fwd_do_te[[2]],rvs_do_te[[2]]) -> do_te_collapse
do_te_collapse %>% count(gene_id) %>% filter(n!=1) %>% select("gene_id") -> gene_id_on_both_strand
do_te_collapse$gene_id %in% gene_id_on_both_strand$gene_id -> fil
do_te_collapse[!fil] %>% select(chr,start,end,gene_id,gene_id_pos,strand) -> do_te_collapse_bed

# statistics
te$gene_id %>% unique() %>% length() -> Gene_number
do_te_collapse_bed$gene_id %>% unique() %>% length() -> Gene_with_dominant_te_number

# per dominant transcripts 
rbind(fwd_do_te[[1]],rvs_do_te[[1]]) -> do_te_per
do_te_per$gene_id %in% gene_id_on_both_strand$gene_id -> fil
do_te_per[!fil] %>% select(chr,start,end,gene_id,gene_id_pos,strand,transcript_id) -> do_te_per_bed

# write
file_path=paste0(args[2],"/number.txt")
cat("Gene_number",Gene_number,"\n",file = file_path,sep = '\t')
cat("Gene_with_dominant_te_number",Gene_with_dominant_te_number,"\n",file = file_path,sep = '\t',append = T)
file_path=paste0(args[2],"/do_te.bed")
fwrite_c(do_te_collapse_bed,file_path)
file_path=paste0(args[2],"/do_te_per.bed")
fwrite_c(do_te_per_bed,file_path)


