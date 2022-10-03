#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:Run ChIPseeker package to plot heatmap for peaks distribution across genomic regions\n',file=stderr())
  cat('Usage:peakCoverage.heatmap.R -i=input.peaks -g=genomicRegions.bed -r=range -c=color -h=heatmap.pdf -p=aggPlot.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tPeak file in bed/narrowPeak/broadPeak format.\n',file=stderr())
  cat('\t-f\tFILE\tThe genomic regions in bed format, the strand information will be ignored\n',file=stderr())
  cat('\t-r\tSTRING\tThe upstream and downstream distance from the region center separated by comma[default:-3000,3000]\n',file=stderr())
  cat('\t-c\tSTRING\tThe color for heatmap[default:red]\n',file=stderr())
  cat('\t-ht\tFILE\tFile name for heatmap.[heatmap.pdf]\n',file=stderr())
  cat('\t-p\tFILE\tFile name for aggregate plot if define.[Optical]\n',file=stderr())
  cat('\t-o\tFILE\tFile name for heatmap matrix.[Optical]\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
heatmapFile="heatmap.pdf"
range="-3000,3000"
color="red"
strand="F"
if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        inFile=arg.split[2]
      }
    }else if(grepl('^-f=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the genomic regions -f')
      }else{
        regionFile=arg.split[2]
      }
    }else if(grepl('^-r=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the ploted range -s')
      }else{
        range=arg.split[2]
      }
    }else if(grepl('^-ht=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      heatmapFile=arg.split[2]
    }else if(grepl('^-p=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      aggFile=arg.split[2]
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      outFile=arg.split[2]
    }
  }
}

library(ChIPseeker)
library(GenomicRanges)
ranges<-as.numeric(strsplit(range,",")[[1]])
region<-read.delim(file=regionFile,header=F)
regionGR=region
for(i in 1:nrow(region)){
  regionL=region[i,3]-region[i,2]
  start=region[i,2]
  regionGR[i,2]<-start+floor((regionL)/2)+ranges[1]
  regionGR[i,3]<-start+1+floor((regionL)/2)+ranges[2]
}
regionGR<-makeGRangesFromDataFrame(regionGR,seqnames.field="V1",start.field="V2",end.field="V3",starts.in.df.are.0based=TRUE)
tagMatrix <- getTagMatrix(inFile, windows=regionGR)
pdf(file=heatmapFile,width=7,height=14)
tagHeatmap(tagMatrix, xlim=c(ranges[1],ranges[2]), color=color)
dev.off()
peak=read.delim(file=inFile,header=F)
if(exists("outFile")){
  write.table(tagMatrix,file=outFile,quote = F,sep="\t",row.names = F,col.names = F)
}
if(exists("aggFile")){
  pdf(file=aggFile)
  meanCol<-colSums(tagMatrix)/nrow(peak)
  plot(c(ranges[1]:(-1),1:ranges[2]),meanCol,xlab="Genomic Region (5'->3')", ylab = "Peak Count Fraction",type="l")
  dev.off()
}
