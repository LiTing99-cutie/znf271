#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script plot GC metrix chart from picard result.\n',file=stderr())
  cat('Usage:picard.GCbiasMetrix.plot.R -i=GCbias.metrix.tsv -o=GCbias.metrix.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe detailed metrixs from picard CollectGcBiasMetrics\n',file=stderr())
  cat('\t-y1\tsTRING\tThe min value for the normalized read coverage(optional)\n',file=stderr())
  cat('\t-y2\tSTRING\tThe max value for the normalized read coverage(Optional)\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput pdf file\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

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
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        out=arg.split[2]
      }
    }
  }
}
data<-read.delim(file=inFile,comment.char = "#",header=T,blank.lines.skip=TRUE)
y1min=min(data$NORMALIZED_COVERAGE)
y1max=max(data$NORMALIZED_COVERAGE)
for(i in 1:length(args)){
  arg=args[i]
  if(grepl('^-y1=',arg)){
    arg.split = strsplit(arg, '=', fixed = T)[[1]]
    y1min=as.numeric(arg.split[2])
  }
  if(grepl('^-y2=',arg)){
    arg.split = strsplit(arg, '=', fixed = T)[[1]]
    y1max=as.numeric(arg.split[2])
  }
}
pdf(file=out)
par(mar=c(4.1,4.1,4.1,4.1))
plot(data$GC,data$NORMALIZED_COVERAGE,xlab="%GC of 100bp windows",ylab="",ylim=c(y1min,y1max),col="#4169E1")
par(new=T)
plot(data$GC,data$WINDOWS/sum(as.numeric(data$WINDOWS)),yaxt="n",ylab="",xlab="",col="#FFAAAA",type="h",ylim=c(0,0.1))
axis(4,at=seq(0,0.1,length.out = 6))
mtext("Normalized coverage",2,col="#4169E1",adj=0.5,padj=-4)
mtext("Fraction of 100bp windows",4,col="#FFAAAA",adj=0.5,padj=4)
dev.off()


