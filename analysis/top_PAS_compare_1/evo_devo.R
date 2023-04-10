
source("/home/user/data2/lit/bin/lit_utils.R")
lib()

# read data in 
evo <- fread(file = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare_1/output/hr.delta_pau_sum.txt")

# process devo data
devo=readRDS("/home/user/data2/lit/project/ZNF271/02-APA-1/output/final_list/res.rm_outlier.rds")

# merge
merge(evo,devo,"gene_name") %>% filter(diff_expr=="NotDiff") %>% mutate(delta_pau_abs=abs(delta_pau)) ->df
merge(evo,devo,"gene_name") %>% mutate(delta_pau_abs=abs(delta_pau))-> df_1


# skewness range distribution
ggplot(evo,aes(x=delta_pau_sum_all))+
  geom_histogram()+
  theme_1()+
  xlab("Sum of delta pau")+ylab("Count") -> p1
ggplot(devo,aes(x=abs(delta_pau)))+
  geom_histogram()+
  theme_1()+
  xlab("Delta PAU")+ylab("Count") -> p3

p1+p3 -> p
ggsave(p,filename = "/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/top_PAS_compare_1/output/dis.pdf",width = 8,height = 5)

df_1 %>% group_by(diff_pau) %>% summarise(mean=mean(delta_pau_sum_all)) 
compare_means(delta_pau_sum_all~diff_pau,df_1)
df %>% group_by(diff_pau) %>% summarise(mean=mean(delta_pau_sum_all)) 
compare_means(delta_pau_sum_all~diff_pau,df)

cor.test(df_1$delta_pau_abs,df_1$delta_pau_sum_all)
cor.test(df$delta_pau_abs,df$delta_pau_sum_all)

df_1 %>% mutate(bin=cut_interval(df_1$delta_pau_abs,5)) ->df_binned
df_1 %>% mutate(bin=cut_number(df_1$delta_pau_abs,5)) ->df_binned_cut_number

compare_means(delta_pau_sum_all~bin,df_binned)
compare_means(delta_pau_sum_all~bin,df_binned_cut_number)
