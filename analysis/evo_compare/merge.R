
# read pairwise results in 
hm <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hm_diff.h_g.txt")
hr <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hr_diff.h_g.txt")
rm <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/rm_diff.h_g.txt")

# merge 
ortho <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/pro_dis_PA_compare/ortholog/hr.txt",header = F)
colnames(ortho) <- c("h_n","r_n")
merge(hm,hr,by="gene_name") %>% merge(ortho,by.x = "gene_name",by.y = "h_n") %>% 
  merge(rm,by.x = "r_n",by.y = "gene_name") %>% mutate(r_n=NULL) ->m

# union
m$diff_union_1=if_else(m$diff_1.y=="Diff"|m$diff_1.x=="Diff"|m$diff_1=="Diff","Diff","NotDiff")
m$diff_union_2=if_else(m$diff_2.y=="Diff"|m$diff_2.x=="Diff"|m$diff_2=="Diff","Diff","NotDiff")

# intersect
m$diff_inter_1=if_else(m$diff_1.x=="Diff"&m$diff_1.y=="Diff"&m$diff_1=="Diff","Diff","NotDiff")
m$diff_inter_2=if_else(m$diff_2.y=="Diff"&m$diff_2.x=="Diff"&m$diff_2=="Diff","Diff","NotDiff")

# fwrite
fwrite(m,sep = '\t',file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/hrm_diff.txt")
