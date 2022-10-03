#!/usr/bin/env Rscript
library(tools)
args <- commandArgs(TRUE)

usage = function(){
  cat('Usage: scriptName.R -p1=mutRate.pdf <mutRate.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-p1 -pdf1\tFILE\tThe output figure in pdf (read1 when PE)[mutRateRead1.pdf]\n',file=stderr())
  cat('\t-p2 -pdf2\tFILE\tThe read2 figure in pdf[mutRateRead2.pdf]\n',file=stderr())
  cat('\t-p3 -pdf3\tFILE\tThe both reads figure in pdf[mutRateBoth.pdf]\n',file=stderr())
  q(save='no')
}

pdf1 = 'mutRateRead1.pdf'
pdf2 = 'mutRateRead2.pdf'
pdf3 = 'mutRateBoth.pdf'

if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-p(df)?1=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p1')
      }else{
        pdf1=arg.split[2]
      }
    }
    if(grepl('^-p(df)?2=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p2')
      }else{
        pdf2=arg.split[2]
      }
    }
    if(grepl('^-p(df)?3=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p3')
      }else{
        pdf3=arg.split[2]
      }
    }
    if(arg == '-h'){
      usage()
    }
  }
}

input = read.delim(file('stdin'), header=FALSE)

pos = input[,1]

if(ncol(input)==5){
  posMuteFrac = input[,5]
  pdf(pdf1, bg='lightyellow',width=15)
  barplot(posMuteFrac, space=0, col='lightblue', names.arg=pos, 
          main = 'Distribution of Mutation Rate', xlab='Position', ylab='Mutation Fraction', las=1)
  abline(h=0.01, lty=2)
}else{
  posMuteRead1Frac = input[,5]
  pdf(pdf1, bg='lightyellow',width=15)
  barplot(posMuteRead1Frac, space=0, col='lightblue', names.arg=pos, 
          main = 'Distribution of Mutation Rate on Read1', xlab='Position', ylab='Mutation Fraction', las=1)
  abline(h=0.01, lty=2)
  posMuteRead2Frac = input[,9]
  pdf(pdf2, bg='lightyellow',width=15)
  barplot(posMuteRead2Frac, space=0, col='lightblue', names.arg=pos, 
          main = 'Distribution of Mutation Rate on Read2', xlab='Position', ylab='Mutation Fraction', las=1)
  abline(h=0.01, lty=2)
  posMuteBothFrac = input[,10]
  pdf(pdf3, bg='lightyellow',width=15)
  barplot(posMuteBothFrac, space=0, col='lightblue', names.arg=pos, 
          main = 'Distribution of Mutation Rate on Both Reads', xlab='Position', ylab='Mutation Fraction', las=1)
  abline(h=0.01, lty=2)
}
