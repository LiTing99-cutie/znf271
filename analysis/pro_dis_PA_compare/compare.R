
library(data.table)
library(dplyr)
library(tidyfst)

func <- function(spe){
  path=paste0("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/",spe,"/PAusage.te.bed6+")
  PA <- fread(path)
  colnames(PA) <- c("chr","start","end","gene_name","cnt","strand")
  
  # coverage on this gene
  PA %>% group_dt(.,gene_name,summarise_dt(sum=sum(cnt))) -> sum
  PA %>% count(gene_name) -> number
  PA %>% count(gene_name,strand) %>% count(gene_name) %>% filter(n==1) -> same_strand
  
  # proximal PA coverage
  filter(PA,strand=="+") %>% group_dt(by = gene_name,filter_dt(start==min(start))) -> tmp.1 
  filter(PA,strand=="-") %>% group_dt(by = gene_name,filter_dt(start==max(start))) -> tmp.2
  rbind(tmp.1,tmp.2) ->tmp
  
  tmp %>% select(gene_name,cnt) %>% merge(sum,"gene_name") %>% 
    merge(number,"gene_name") %>% merge(same_strand,"gene_name") %>% filter(n.x>1) %>%
    mutate(dis=sum-cnt) %>% select(gene_name,cnt,dis) %>% rename(pro=cnt) -> final_2
  
  l <- list(final_2,number)
  
  return(l)
}

func("h") ->h
func("r") ->r
func("m") ->m



# merge
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hrm.ortholog.txt",header = F)

##### 3 species #####
#1. merge to compare proximal and distal PAS coverage
h[[1]] %>% merge(ortho,by.x = "gene_name",by.y = "V1") %>% merge(r[[1]],by.x = "V2",by.y = "gene_name") %>% 
  merge(m[[1]],by.x = "V3",by.y = "gene_name") %>% select(gene_name,pro.x,dis.x,pro.y,dis.y,pro,dis) ->compare

#2. merge to compare PAS number
h[[2]] %>% merge(ortho,by.x = "gene_name",by.y = "V1") %>% merge(r[[2]],by.x = "V2",by.y = "gene_name") %>% 
  merge(m[[2]],by.x = "V3",by.y = "gene_name") %>% select(gene_name,n.x,n.y,n) -> cp_pas_n

# test
df <- data.frame()
for(i in 1:nrow(compare)){
  compare[i,2:7] %>% unlist() %>% matrix(,nrow = 2) %>% chisq.test() %>% .$p.value -> p_value
  df <- rbind(df,p_value)
}
colnames(df) <- "p_value"
mutate(df,p_adjust=p.adjust(df$p_value)) ->df
cbind(compare,df) -> cp

lo_1 <- cp_pas_n$n.x==1 & cp_pas_n$n.y==1 & cp_pas_n$n.x==1
lo_2 <- cp_pas_n$n.x>1 & cp_pas_n$n.y>1 & cp_pas_n$n.x>1
lo_3 <- !lo_1 & !lo_2
cp_pas_n[lo_3] -> cp_pas_n_diff

fwrite(cp,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp.txt",col.names = T,sep = '\t')
fwrite(cp_pas_n_diff,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n_diff.txt",col.names = T,sep = '\t')
fwrite(cp_pas_n,file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/cp_pas_n.txt",col.names = T,sep = '\t')

##### 2 species #####

test_two <- function(s_1,s_2,c_1,c_2){
s_1 %>% merge(ortho,by.x = "gene_name",by.y = c_1) %>% merge(s_2,by.x = c_2,by.y = "gene_name") %>% select(gene_name,pro.x,dis.x,pro.y,dis.y) ->compare
df <- data.frame()
for(i in 1:nrow(compare)){
  compare[i,2:5] %>% unlist() %>% matrix(,nrow = 2) %>% fisher.test() %>% .$p.value -> p_value
  df <- rbind(df,p_value)
}
colnames(df) <- "p_value"
mutate(df,p_adjust=p.adjust(df$p_value)) ->df
cbind(compare,df) -> cp
mutate(cp,orien=dis.y/(dis.y+pro.y)-dis.x/(dis.x+pro.x)) -> cp
return(cp) 
}


test_two(h[[1]],r[[1]],"V1","V2") ->hr
test_two(h[[1]],m[[1]],"V1","V3") ->hm
test_two(r[[1]],m[[1]],"V2","V3") ->rm

# 2207
union(filter(hr,p_adjust<0.05)$gene_name,filter(hm,p_adjust<0.05)$gene_name) %>% union(.,filter(rm,p_adjust<0.05)$gene_name) -> union_hrm
# 1237
filter(cp,p_adjust<0.05)$gene_name -> spe_three_res
# 1054
intersect(spe_three_res,union_hrm) %>% length()

fwrite(data.frame(gene_name=union_hrm),file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/output/union_hrm.txt",col.names = T,sep = '\t')
