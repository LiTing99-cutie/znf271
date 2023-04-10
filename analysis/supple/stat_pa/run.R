args <- commandArgs(T)
library("data.table")
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

# PAS number and gene number
len_cnt <- fread(args[1])
len_cnt %>% count(V1) -> gene_pas_n
cat("pa_isoforms_number",nrow(len_cnt),"\n",file = "./stat.txt",sep = '\t')
cat("gene_number",nrow(gene_pas_n),"\n",file = "./stat.txt",append = T,sep = '\t')

# PAS number on each gene histogram
gene_pas_n %<>% mutate(n_sta=case_when(n==1~"1",
                                      n==2~"2",
                                      n==3~"3",
                                      n>=4~"4+"))
table(gene_pas_n$n_sta) %>% data.frame() %>% ggplot(aes(x=Var1,y=Freq))+geom_col(width = 0.6)+
theme_1()+
xlab("PAS number")+ylab("Gene number")+
geom_text(mapping = aes(label=Freq),size=4,vjust=-0.5) -> p
ggsave("./pas_number.pdf",width = 5,height = 5) 
p
dev.off()

# PA length and distance between proximal PAS and distal PAS

len_cnt %>% ggplot(aes(x=V2))+geom_histogram(color="red")+theme_1()+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  xlab("PA length")+ylab("Number") -> p
ggsave("./pa_length.pdf",width = 5,height = 5)
p
dev.off()


len_cnt %>% group_by(V1) %>% summarise(dis=range(V2)) %>% filter(dis!=0) %>% 
ggplot(aes(x=dis))+geom_histogram(color="red")+theme_1()+
scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
xlab("Distance between proximal PAS and distal PAS")+ylab("Number") -> p
ggsave("./distance.pdf",width = 6,height = 5)
p
dev.off()
