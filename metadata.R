rm(list=ls())

setwd("/Users/katherine/project/ZNF271/apa")

# library R packages
library(data.table)
library(tidyverse)
library(magrittr)
library(dplyr)

# reorganize data

metadata <- fread("/Users/katherine/Downloads/E-MTAB-6814.sdrf.txt")
metadata %<>% select("Source Name","Characteristics[sex]" ,
                     "Characteristics[developmental stage]",
                     "Characteristics[age]","Unit [time unit]",
                     "Characteristics[organism part]")
colnames(metadata) <- c("Id","Sex","Developmental_Stage","Age","Unit","organism")
metadata$Developmental_Stage %<>% gsub("week post conception","pcw",.)
metadata$Developmental_Stage <- factor(metadata$Developmental_Stage,
                                       levels = c("neonate","4 pcw","5 pcw","6 pcw","7 pcw","8 pcw","9 pcw","10 pcw",
                                                  "11 pcw","12 pcw","13 pcw","16 pcw","18 pcw","19 pcw","20 pcw","infant","toddler",
                                                 "school age child","adolescent","young adult","middle adult","elderly"))

saveRDS(metadata,file = "metadata.rds")
