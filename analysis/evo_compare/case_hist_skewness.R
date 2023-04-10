

h <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/h/h.cnt.bed")
m <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/m/m.cnt.bed")
r <- fread("/home/user/data2/lit/project/ZNF271/02-APA-1/analysis/evo_compare/r/r.cnt.bed")

# max length in different species
# 2135
# 3016
# 2721

GetHist <- function(df,gene){
  df %>% filter(V1==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"V2"],times=input[i,"V3"])) %>% unlist()->Terminal_exon_length
  }
  hist(Terminal_exon_length)
}

GetHist_modi <- function(df,gene){
  df %>% filter(V1==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"V2"],times=input[i,"V3"])) %>% unlist()->Terminal_exon_length
  }
  Terminal_exon_length <- c(1,3016,Terminal_exon_length)
  hist(Terminal_exon_length)
}

GetHist(h,"ZNF271P")
GetHist(r,"ZNF271")
GetHist(m,"Zfp35")

GetHist_modi(h,"ZNF271P")
GetHist_modi(r,"ZNF271")
GetHist_modi(m,"Zfp35")

GetSknewness <- function(df,gene){
  df %>% filter(V1==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"V2"],times=input[i,"V3"])) %>% unlist()->Terminal_exon_length
  }
  skewness(Terminal_exon_length)
}

GetSknewness_modi <- function(df,gene){
  df %>% filter(V1==gene) ->input
  arr <- c()
  for (i in 1:nrow(input)){
    arr <- c(arr,rep(input[i,"V2"],times=input[i,"V3"])) %>% unlist()->Terminal_exon_length
  }
  Terminal_exon_length <- c(1,3016,Terminal_exon_length)
  skewness(Terminal_exon_length)
}
GetSknewness(h,"ZNF271P")
GetSknewness(r,"ZNF271")
GetSknewness(m,"Zfp35")

GetSknewness_modi(h,"ZNF271P")
GetSknewness_modi(r,"ZNF271")
GetSknewness_modi(m,"Zfp35")