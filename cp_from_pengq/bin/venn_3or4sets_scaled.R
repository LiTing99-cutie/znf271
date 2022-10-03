#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script will creat a proportional venn diagram for 3/4 sets by using R packages eulerr and Vennerable.\n',file=stderr())
  cat('Notice: Need R version >= 3.3.0 for 4 sets\n',file=stderr())
  cat('Usage:venn_mutipleRegions.R -i=file1,file2,file3 -n=sample1,sample2,sample3 -o=venn.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\t\tFILE\tThe input list files separated by comma,each with one column.\n',file=stderr())
  cat('\t-n\t\tSTRING\tNames for each files separated by comma.\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput file name[venn.pdf]\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

out="venn.pdf"

if(length(args)>1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -i')
      }else{
        inFiles=arg.split[2]
      }
    }
    if(grepl('^-n=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -n')
      }else{
        names=arg.split[2]
      }
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -o')
      }else{
        out=arg.split[2]
      }
    }
  }
}else if(length(args)==0 || args[1]=="-h"){
  usage()
}

files=strsplit(inFiles,",")[[1]]
names=strsplit(names,",")[[1]]
library("Vennerable")
file1<-read.delim(file=files[1],header=F)
file2<-read.delim(file=files[2],header=F)
file3<-read.delim(file=files[3],header=F)
if(length(files)==4){
  library("eulerr")
  file4<-read.delim(file=files[4],header=F)
  list<-list(file1$V1,file2$V1,file3$V1,file4$V1)
  names(list)<-names
  overlap<-Venn(list)
  result<-names(overlap@IntersectionSets)
  #Claulate logical matrix
  set<-grepl("1",strsplit(result[2],"")[[1]])
  length<-length(overlap@IntersectionSets[2][[1]])
  matrix<-matrix(set,ncol=4,nrow=length,byrow=T)
  for(i in 3:length(result)){
    set<-grepl("1",strsplit(result[i],"")[[1]])
    length<-length(overlap@IntersectionSets[i][[1]])
    if(length>0){
      mat<-matrix(set,ncol=4,nrow=length,byrow=T)
      matrix<-rbind(matrix,mat)
    }
  } 
  colnames(matrix)<-names
  venn<-euler(matrix)
  pdf(out)
  plot(venn,counts=T,fill=c("red","green","blue","purple"))
}else{
  list<-list(file1$V1,file2$V1,file3$V1)
  names(list)<-names
  overlap<-Venn(list)
  pdf(out)
  plot(overlap, doWeights = TRUE, show=list(Faces=F))
}

dev.off()
