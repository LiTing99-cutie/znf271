
hr_pro <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.pro.RheMac8.bed6+",header = F)
hr_dis <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.dis.RheMac8.bed6+",header = F)
hm_pro <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.pro.Mm10.bed6+",header = F)
hm_dis <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.dis.Mm10.bed6+",header = F)
h_pro <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.pro.bed6+",header = F)
h_dis <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.dis.bed6+",header = F)

unmapped_r <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.ToRheMac8.unmapped.gene.lst",header = F)
unmapped_m <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/hg38.ToMm10.unmapped.gene.lst",header = F)
ortholog <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)

h_pro %>% group_by(V4) %>% summarise(cnt=sum(V11)) -> h_pro_cnt
h_dis %>% group_by(V4) %>% summarise(cnt=sum(V11)) -> h_dis_cnt
merge(h_pro_cnt,h_dis_cnt,by="V4",all = T) ->cnt_h

MergeAndCnt_r <- function(df){
  merge(df,ortholog,by.x = "V4",by.y = "V1") %>% filter(V10==V2.y) %>% group_by(V4) %>% summarise(cnt=sum(V11)) ->res
  return(res)
}
MergeAndCnt_m <- function(df){
  merge(df,ortholog,by.x = "V4",by.y = "V1") %>% filter(V10==V3.y) %>% group_by(V4) %>% summarise(cnt=sum(V11)) ->res
  return(res)
}
RmUnmap <- function(df,unmapped){
  df$V4 %in% unmapped$V1 -> inx
  df[!inx]->df_clean
  return(df_clean)
}
Ortho_Cnt <- function(pro,dis,MergeAndCnt,unmapped){
  RmUnmap(pro,unmapped) -> pro
  RmUnmap(dis,unmapped) -> dis
  MergeAndCnt(pro) ->pro_cnt
  MergeAndCnt(dis) ->dis_cnt
  merge(pro_cnt,dis_cnt,by="V4",all = T) ->cnt
  cnt[is.na(cnt)] <- 0
  return(cnt)
}
COMPARE <- function(cnt_1,cnt_2,compare){
  merge(cnt_1,cnt_2,by="V4") ->compare
  df <- data.frame()
  for(i in 1:nrow(compare)){
    compare[i,2:5] %>% unlist() %>% matrix(,nrow = 2) %>% fisher.test() %>% .$p.value -> p_value
    df <- rbind(df,p_value)
  }
  colnames(df) <- "p_value"
  mutate(df,p_adjust=p.adjust(df$p_value)) ->df
  cbind(compare,df) -> cp
  
  return(cp)
}

Ortho_Cnt(hr_pro,hr_dis,MergeAndCnt_r,unmapped_r) -> cnt_r
Ortho_Cnt(hm_pro,hm_dis,MergeAndCnt_m,unmapped_m) -> cnt_m

COMPARE(cnt_r,cnt_m,) -> cp_rm
COMPARE(cnt_h,cnt_m) -> cp_hm
COMPARE(cnt_h,cnt_r) -> cp_hr
  
l_1 <- filter(cp_hm,p_adjust<0.05) %>% select("V4")
l_2 <- filter(cp_hr,p_adjust<0.05) %>% select("V4")
l_3 <- filter(cp_rm,p_adjust<0.05) %>% select("V4")

union(l_1,l_2) %>% union(.,l_3) %>% rename(gene_name=V4) %>% 
  fwrite("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_region_PA_compare/analysis/evo_diff_pro_dis_region.1.lst")
