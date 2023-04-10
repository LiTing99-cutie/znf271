
setwd("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/typeII-III")
end <- fread("3end_validated_typeII_III.merge.txt")
iso_seq <- fread("pas_supported_sample/supported_sample.final.txt",fill = TRUE)

merge(end,iso_seq,by.x = "V5",by.y = "V1",all.y = T) -> tmp

union(end$V5,na.omit(iso_seq)$V1) -> g_h_c

tmp[tmp$V5 %in% g_h_c] ->tmp_1
