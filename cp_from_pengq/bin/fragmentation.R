#!/bin/env Rscript
args <- commandArgs(TRUE)

usage = function(){
  cat('Usage: scriptName.R -p=output.pdf <fragmentation.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-p -pdf\tFILE\tThe output figure in pdf[fragmentation.pdf]\n',file=stderr())
  cat('\t-h\t\tShow help\n',file=stderr())
  q(save='no')
}

myPdf = 'fragmentation.pdf'

if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    if(arg == '-h'){
      usage()
    }
    if(grepl('^-p(df)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p')
      }else{
        myPdf=arg.split[2]
      }
    }
  }
}

pdf(myPdf, bg = 'lightyellow', width = 10)

data = read.delim(file('stdin'), header = F)
barplot(data[,2], names.arg = data[,1], main = 'Fragmentation', xlab = 'Interval', ylab = 'Fraction')
abline(h=0.04, lty = 2, lwd = 1)
