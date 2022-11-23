
library(stringr)
library(dplyr)
library(magrittr)

#----------------
metadata <- readRDS(file="/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.rds")
colnames(metadata) <- c("sample","Sex","Developmental_Stage","Age","Unit","organism")
write.table(metadata,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.txt",quote = FALSE,sep = '\t',row.names = FALSE)
saveRDS(metadata,"/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata.rds")

#----------------
period <- read.table("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period_description.awk.txt",sep = '\t')
colnames(period) <- c("Developmental_Stage","period")
add <- data.frame(Developmental_Stage=c("neonate","infant","toddler","school age child","adolescent","young adult","middle adult","elderly"),
                  period=rep("after birth")) 
rbind(period,add) -> period
period$period <- str_to_sentence(period$period)

period$period %<>% factor(.,level=c("Embryonic","Early fetal","Early mid-fetal","Late mid-fetal","Late fetal","After birth"))
period %<>% mutate(.,period.abbre=case_when(period=="Embryonic"~"Embryonic",
                                           period=="Early fetal" | period=="Early mid-fetal" | period=="Late fetal" | period=="Late mid-fetal"~"Fetal",
                                           period=="After birth"~"After birth"))
period$period.abbre %<>% factor(.,level=c("Embryonic","Fetal","After birth"))
period %<>% mutate(.,period.abbre.abbre=case_when(period.abbre=="Fetal"| period.abbre=="Embryonic" ~"Before birth",
                                           period.abbre=="After birth"~"After birth"))
period$period.abbre.abbre %<>% factor(.,level=c("Before birth","After birth"))
saveRDS(period,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period_description.rds")
write.table(period,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/period.txt",quote = FALSE,sep = '\t',row.names = FALSE)

#------------
merge(metadata,period,by="Developmental_Stage") -> metadata_period
merge(metadata,period,by="Developmental_Stage") %>% 
  write.table(.,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/metadata_period.txt",quote = FALSE,sep = '\t',row.names = FALSE)

#------------
frag_score <- read.table("/home/user/data2/lit/project/ZNF271/02-APA/analysis/fragmentation/brain/fragmentation.score.txt")
colnames(frag_score) <- c("sample","frag_score")
frag_score$sample %<>% gsub("E-MTAB-6814.","",.)
saveRDS(frag_score,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.rds")
write.table(frag_score,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/frag_score.txt",quote = FALSE,sep = '\t',row.names = FALSE)

#------------
libSize <- read.table("/home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapping/brain/uniq/bam_stat/library_size.txt")
colnames(libSize) <- c("sample","mapped_reads")
libSize$sample <- gsub("E-MTAB-6814.","",libSize$sample) %>% gsub(".sorted","",.)


  

merge(libSize,metadata_period,by= "sample") %>% group_by(period.abbre.abbre) %>% summarise(mapped_reads_s=sum(mapped_reads)/10^6)->tmp_1
colnames(tmp_1) <- c("period","mapped_reads_s")
merge(libSize,metadata_period,by= "sample") %>% group_by(period) %>% summarise(mapped_reads_s=sum(mapped_reads)/10^6) ->tmp_2
rbind(tmp_1[1,],tmp_2) %>% 
write.table(.,file = "/home/user/data/lit/project/ZNF271/data/rna-seq/brain/mapped_m_read_period.txt",quote = FALSE,sep = '\t',row.names = FALSE)
