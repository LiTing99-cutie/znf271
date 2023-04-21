args=commandArgs(T)

# 1 file
# 2 cutoff of overall coverage
# 3 cutoff of single PAS PAU
# 4 cutoff of single PAS count
# 5 output path

if(is.na(args[1])){
args <- c()
args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/human/PAusage.bed6+"
args[2] <- 0
args[3] <- 0
args[4] <- 2
args[5] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/human/"
}
pas <- fread(args[1])


if(args[2]==0 & args[3]==0){
  pas %>% filter(V5>=as.numeric(args[4])) -> res
}else{
  pas %>% group_by(V4) %>% summarise(coverage=sum(V5)) %>% filter(coverage>=as.numeric(args[2])) -> fil_cov
  merge(fil_cov,pas) %>% filter(V7>as.numeric(args[3])) %>% select(V1,V2,V3,V4,V5,V6,V7) -> res
}

# PAS number and gene number
res %>% count(V4) -> gene_pas_n
cat("pa_isoforms_number",nrow(res),"\n",file = paste0(args[5],"/stat.txt"),sep = '\t')
cat("gene_number",nrow(gene_pas_n),"\n",file = paste0(args[5],"/stat.txt"),append = T,sep = '\t')

# PAS number on each gene histogram
gene_pas_n %<>% mutate(n_sta=case_when(n==1~"1",
                                       n==2~"2",
                                       n==3~"3",
                                       n>=4~"4+"))
table(gene_pas_n$n_sta) %>% data.frame() %>% ggplot(aes(x=Var1,y=Freq))+geom_col(width = 0.6)+
  theme_1()+
  xlab("PAS number")+ylab("Gene number")+
  geom_text(mapping = aes(label=Freq),size=4,vjust=-0.5) -> p
ggsave(paste0(args[5],"/pas_number.pdf"),width = 5,height = 5) 
p
dev.off()

# write to file
fwrite(res,file = paste0(args[5],"/PAusage.fil.bed6+"),sep='\t',col.names = F)
