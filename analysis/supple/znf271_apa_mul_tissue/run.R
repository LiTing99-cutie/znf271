

setwd("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/supple/znf271_apa_mul_tissue")
res <- fread("apa_m1.PDUI.0.NR_024566.ENSG00000257267.txt")

# Merge tissues of different part
res %<>%
  mutate(Tissue_c = case_when(
    grepl("Adipose", Tissue) ~ "Adipose",
    grepl("Artery", Tissue) ~ "Artery",
    grepl("Brain", Tissue) ~ "Brain",
    grepl("Breast", Tissue) ~ "Breast",
    grepl("Cervix", Tissue) ~ "Cervix",
    grepl("Colon", Tissue) ~ "Colon",
    grepl("Esophagus", Tissue) ~ "Esophagus",
    grepl("Heart", Tissue) ~ "Heart",
    grepl("Kidney", Tissue) ~ "Kidney",
    grepl("Muscle", Tissue) ~ "Muscle",
    grepl("Nerve", Tissue) ~ "Nerve",
    grepl("Skin", Tissue) ~ "Skin",
    grepl("Small Intestine", Tissue) ~ "Small Intestine",
    TRUE ~ Tissue
  ))

# Exclude cells
res %<>% filter(!grepl("Cells",Tissue_c))

# Reorder tissue according to median value of PDUI
res %>% group_by(Tissue_c) %>% summarise(median=median(PDUI)) %>% arrange(median) %>% .$"Tissue_c" -> factor
res$Tissue_c <- factor(res$Tissue_c,levels=factor)

# number of each group
res %>% count(Tissue_c) -> sample_cnts

# plot
p <- ggplot(res,aes(x=Tissue_c,y=PDUI))+geom_boxplot(fill=colorRampPalette(brewer.pal(9,"Set1"))(30),width=0.5)+theme_1()+xlab("Tissue")+
  stat_compare_means(label.y.npc = 0.1,label.x.npc = 1,size=5)+
  scale_x_discrete(labels=paste0(sample_cnts$Tissue_c," (n=",sample_cnts$n,")"))+
  coord_flip()
ggsave(p,filename = "merge.output.pdf",height = 10,width = 10)