rm(list=ls())

setwd("/Users/katherine/project/ZNF271/apa")

# library R packages
library(data.table)
library(tidyverse)
library(magrittr)
library(dplyr)

# load data
output <- fread("/Users/katherine/Downloads/result.brain.txt")
metadata <- readRDS(file="metadata.rds")

# reorganize data

## output
colnames(output) <- c("Id","Down","Up")
output$Id <- sub("E-MTAB-6814.","",output$Id)
mutate(output,"PA_usage"=Down/Up) -> output.pa.usage

# sample number
if(F){
  metadata[grep("brain|cerebellum",metadata$`Id`,ignore.case = TRUE),] -> metadata.brain
  sample_number <- as.data.frame(table(metadata.brain$Developmental_Stage),col.names=c("stage","number")) 
  ggplot(sample_number,aes(Var1,Freq))+geom_bar(stat = "identity")+theme_classic()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
}

# merge output and metadata
merge(output.pa.usage,metadata,by="Id") -> output.pa.usage.metadata

###### all the developmental stages ####
pdf(file = "brain.custom.pdf",height = 4,width = 5)
ggplot(output.pa.usage.metadata,aes(Developmental_Stage,PA_usage,group=Developmental_Stage))+
  geom_boxplot(outlier.colour = NA,outlier.fill = NA)+theme_classic()+
  geom_jitter(size=1,color="gray40",alpha=0.6,position=position_jitter(width = 0.2,height = 0))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"),
        axis.title = element_text(size=8,face = "bold"))+
  labs(x="Developmental Stage",y = "PA usage",title = "APA usage of ZNF271 in developmental brain")+
  facet_grid(. ~ organism)
dev.off()




###### subset developmental stages correspond to before studies ####

output.pa.usage.metadata %>% 
  .[.$Developmental_Stage %in% c("4 pcw","5 pcw","11 pcw","16 pcw","20 pcw"),] -> 
  output.pa.usage.metadata.subset

pdf(file = "brain.custom.subset_age.pdf",height = 4,width = 5)
ggplot(output.pa.usage.metadata.subset,aes(Developmental_Stage,PA_usage,group=Developmental_Stage))+
  geom_boxplot(outlier.colour = NA,outlier.fill = NA)+theme_classic()+
  geom_jitter(size=1,color="gray40",alpha=0.6,position=position_jitter(width = 0.2,height = 0))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        plot.title = element_text(hjust = 0.5,size = 8,face = "bold"),
        axis.title = element_text(size=8,face = "bold"))+
  labs(x="Developmental Stage",y = "PA usage",title = "APA usage of ZNF271 in developmental brain")+
  facet_grid(. ~ organism)
dev.off()



