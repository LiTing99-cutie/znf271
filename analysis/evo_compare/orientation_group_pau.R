
# all positive or negative vs not (PAU)

# process devo data
res=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")
evo_1 <- fread(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_evo.txt")
evo_3 <- readRDS(file ="/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/skewness_solu_1_2_evo.rds")

merge(res,evo_1,"gene_name") -> m
mutate(m,abs_delta_pau=abs(delta_pau)) -> m

compare_means(abs_delta_pau~diff,m)
m %>% group_by(diff) %>% summarise(mean=mean(abs_delta_pau))

merge(res,evo_3[[1]],"gene_name") -> m
mutate(m,abs_delta_pau=abs(delta_pau)) -> m

compare_means(abs_delta_pau~diff,m)
m %>% group_by(diff) %>% summarise(mean=mean(abs_delta_pau))