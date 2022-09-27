rm(list=ls())

setwd("/Users/katherine/project/ZNF271/apa")

# library packages
library(data.table)
library(tidyverse)
library(magrittr)
library(dplyr)

# load data
output <- fread("/Users/katherine/project/ZNF271/apa/pau_results.txt")
metadata <- readRDS(file="metadata.rds")

# reorganize data

## output
select(output,-c(1:12)) %>% .[5,] -> output.clean
melt(output.clean) -> output.clean.long
output.clean.long %>% .[grep("PAU",.$variable,ignore.case = TRUE),] -> output.clean.long.PAU
output.clean.long.PAU$variable %<>% sub("E-MTAB-6814.","",.) %>% sub(".PAU","",.) 
mutate(output.clean.long.PAU,"PA_usage"=100-value) -> output.clean.long.PAU.distal

# merge output and metadata
merge(output.clean.long.PAU.distal,metadata,by.x="variable",by.y="Id") -> output.clean.long.PAU.distal.metadata

# plot

###### all the developmental stages ####
pdf(file = "brain.qapa.pdf",height = 4,width = 5)
ggplot(output.clean.long.PAU.distal.metadata,aes(Developmental_Stage,PA_usage,group=Developmental_Stage))+
  geom_boxplot(outlier.colour = NA,outlier.fill = NA)+theme_classic()+
  geom_jitter(size=1,color="gray40",alpha=0.6,position=position_jitter(width = 0.2,height = 0))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"),
        axis.title = element_text(size=8,face = "bold"))+
  labs(x="Developmental Stage",y = "PA usage",title = "APA usage of ZNF271 in developmental brain")+
  facet_grid(. ~ organism)
dev.off()




###### subset developmental stages correspond to before studies ####

output.clean.long.PAU.distal.metadata %>% 
  .[.$Developmental_Stage %in% c("4 pcw","5 pcw","11 pcw","16 pcw","20 pcw"),] -> 
  output.clean.long.PAU.distal.metadata.subset

pdf(file = "brain.qapa.subset_age.pdf",height = 4,width = 5)
ggplot(output.clean.long.PAU.distal.metadata.subset,aes(Developmental_Stage,PA_usage,group=Developmental_Stage))+
  geom_boxplot(outlier.colour = NA,outlier.fill = NA)+theme_classic()+
  geom_jitter(size=1,color="gray40",alpha=0.6,position=position_jitter(width = 0.2,height = 0))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"),
        axis.title = element_text(size=8,face = "bold"))+
  labs(x="Developmental Stage",y = "PA usage",title = "APA usage of ZNF271 in developmental brain")+
  facet_grid(. ~ organism)
dev.off()

