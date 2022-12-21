args <- commandArgs(TRUE)
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyfst))

file=read.table(args[1])

tmp <- mutate(file,gene_id=str_split_fixed(file$V1,'_',2)[,1]) %>% group_dt(gene_id,filter_dt(V2==max(V2)))
prefix=gsub(".txt","",args[1])
write.table(tmp,quote = F,col.names = F,row.names = F,file = paste0(prefix,"_filter.txt"))