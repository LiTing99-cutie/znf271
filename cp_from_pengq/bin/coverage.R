#!/usr/bin/env Rscript
library(tools)
args <- commandArgs(TRUE)

usage = function(){
  cat('Usage: scriptName.R -p=output.pdf <cov_all.tsv >cov_all_cum.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-p -pdf\t\tFILE\tThe output figure in pdf[cov.pdf]\n',file=stderr())
  cat('\t-c -minCov\tINT\tThe minimal coverage to be taken as sufficient depth[10]\n',file=stderr())
  q(save='no')
}

myPdf = 'cov.pdf'
minCov = 10

if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-p(df)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p')
      }else{
        myPdf=arg.split[2]
      }
    }
    if(grepl('-c=', arg) || grepl('-minCov=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -c')
      }else{
        minCov=as.numeric(arg.split[2])
      }
    }
    if(arg == '-h'){
      usage()
    }
  }
}

input = read.delim(file('stdin'), header=F)
coverage = input[,1]
fre = input[,3]

pdf(myPdf, bg = 'lightyellow', width = 10)
cdf = rev(cumsum(rev(fre)))

plot(coverage, cdf, type = 'l', main = 'Target Coverage', xlab = 'Depth', ylab = 'Cumulative Fraction', las = 1, lwd = 2, tcl = -0.5)
lines(c(minCov, minCov), c(0, cdf[minCov+1]), lty = 2, lwd = 1)
grid (10, 10 , lwd = 1)
text(minCov, cdf[minCov+1], cdf[minCov+1], adj = c(-0.1, -0.5))
axis(1, at = c(minCov))
write.table(cbind(coverage, cdf), stdout(), row.names = F, col.names = F, sep = "\t")
