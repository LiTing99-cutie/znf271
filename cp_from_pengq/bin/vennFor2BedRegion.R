#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script find the overlap regions for two bed file and plot a venn diagram.(Please insure bedtools are in your current PATH)\n',file=stderr())
  cat('Usage:vennForBedRegion.R -b1=<bed1> -b2=<bed2> -c=[6] -n1="sample1" -n2="sample2" -o=venn.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-b1\t\tFILE\tThe first bed file\n',file=stderr())
  cat('\t-b2\t\tFILE\tThe second bed file\n',file=stderr())
  cat('\t-c\tINT\tThe total column number for the first bed file[6]\n',file=stderr())
  cat('\t-n1\tSTRING\tThe label for the first bed[sample1]\n',file=stderr())
  cat('\t-n2\tSTRING\tThe label for the second bed[sample2]\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput file name[venn.pdf]\n',file=stderr())
  q(save="no")
}

col=6
name1="sample1"
name2="sample2"
out="venn.pdf"

if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-b1=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -b1')
      }else{
        bed1=arg.split[2]
      }
     }
      if(grepl('^-b2=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
          stop('Please specify the value of -b2')
        }else{
          bed2=arg.split[2]
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
      if(grepl('^-c=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        col=as.numeric(arg.split[2])
      }
      if(grepl('^-n1=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        name1=arg.split[2]
      }
      if(grepl('^-n2=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        name2=arg.split[2]
      }
  	if(grepl('^-h', arg)){
	 usage()
	 }

  }
}else{
  usage()
}

system(paste("bedtools intersect -wo -a ",bed1," -b ",bed2," >",name1,"_",name2,".intersect.bed",sep=""))
intersect=paste(name1,"_",name2,".intersect.bed",sep="")

num1=as.numeric(system(paste("cat ",bed1,"|wc -l"),intern = T))
num2=as.numeric(system(paste("cat ",bed2,"|wc -l"),intern = T))
overlap=as.numeric(system(paste("cut -f1-",col," ",intersect,"|sort|uniq|wc -l",sep=""),intern=T))

library("VennDiagram")
A=seq(1,num1,1)
B=seq(-(num2-overlap-1),overlap,1)

T<-venn.diagram(list(A,B),category.names=c(name1,name2),file=NULL,fill=rainbow(2),alpha=c(0.5,0.5),cat.pos=6)
pdf(out)
grid.draw(T)
dev.off()


