source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args=commandArgs(T)

# default
if(is.na(args[1])){
  args <- c()
  args[1] <- "/home/user/data2/lit/project/ZNF271/02-APA-1/PAUsage_bed/output_m/fragmentation.score.txt"
}