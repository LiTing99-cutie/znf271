

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(pheatmap))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/h/cds_type_pa_count.txt") -> h

fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/m/cds_type_pa_count.txt") -> m

fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/r/cds_type_pa_count.txt") -> r

fread("/home/user/data2/lit/project/ZNF271/02-APA-1/evolution-diff/clean/homolog.add.txt",
      header = F,col.names = c("h_g_n","r_g_n","m_g_n")) ->homolog

merge(m,select(homolog,h_g_n,m_g_n),by.x="gene_name",by.y = "m_g_n") %>% mutate(gene_name=NULL) %>% rename(gene_name=h_g_n)->m_c_g_n
merge(r,select(homolog,h_g_n,r_g_n),by.x="gene_name",by.y = "r_g_n") %>% mutate(gene_name=NULL) %>% rename(gene_name=h_g_n)->r_c_g_n

m_w <-dcast(m_c_g_n,gene_name~cds_type,value.var="count_sum")
h_w <-dcast(h,gene_name~cds_type,value.var="count_sum")
r_w <-dcast(r_c_g_n,gene_name~cds_type,value.var="count_sum")


#### compare between human and rhesus ####
merge(r_w,h_w,"gene_name")-> h_r_com

h_r_com[is.na(h_r_com)] <- 0

p_value=data.frame()
for (i in 1:nrow(h_r_com)){
  h_r_com[i,2:5]  %>% unlist()%>% matrix(,nrow = 2) %>% fisher.test() %>% .$p.value %>% data.frame() ->tmp
  p_value=rbind(p_value,tmp)
}

colnames(p_value) <- "p_value"
cbind(h_r_com,p_value) -> res_hr

#### compare between mouse and rhesus ####
merge(r_w,m_w,"gene_name")-> r_m_com

r_m_com[is.na(r_m_com)] <- 0

p_value=data.frame()
for (i in 1:nrow(r_m_com)){
  r_m_com[i,2:5]  %>% unlist()%>% matrix(,nrow = 2) %>% fisher.test() %>% .$p.value %>% data.frame() ->tmp
  p_value=rbind(p_value,tmp)
}

colnames(p_value) <- "p_value"
cbind(r_m_com,p_value) -> res_rm
#### merge res_hr and res_rm ####
merge(res_hr,res_rm,by="gene_name")->final
final %<>%  select(-7,-8)
colnames(final) <- c("gene_name","alternative_r","primary_r","alternative_h","primary_h","p_value_hr","alternative_m","primary_m","p_value_rm")
print("Significantly differential between human and rhesus")
filter(final,p_value_hr<0.05)

#### draw a heatmap ####
d_pau_r <- final[,3]/(final[,2]+final[,3])
d_pau_h <- final[,5]/(final[,4]+final[,5])
d_pau_m <- final[,8]/(final[,7]+final[,8])
d_pau <- cbind(d_pau_h,d_pau_r,d_pau_m)
d_pau <- t(d_pau) # matrix
rownames(d_pau) <- c("Human","Rhesus","Mouse")
colnames(d_pau) <- final$gene_name


# colnames(d_pau) <- c("Human","Rhesus","Mouse")
# rownames(d_pau) <- final$gene_name
# d_pau <- t(d_pau)

# set breaks
bk=seq(0,1,by=0.01)
ph <- pheatmap(d_pau, color = colorRampPalette(c("white", "#bd6a9a"))(100),
         cluster_cols = FALSE,cluster_rows = FALSE,angle_col = 45,breaks = bk,fontsize = 8,show_rownames = T,show_colnames = T)
ggsave(ph,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/evo_heatmap.pdf",width = 5,height =3,units = "in")

ph_n <- pheatmap(d_pau, color = colorRampPalette(c("white", "#bd6a9a"))(100),
         cluster_cols = FALSE,cluster_rows = FALSE,angle_col = 45,breaks = bk,fontsize = 8,display_numbers = TRUE)
ggsave(ph_n,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/output/R/evo_heatmap_n.pdf",width = 5,height =3,units = "in")
