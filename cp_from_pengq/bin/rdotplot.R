################################################
#File Name: rdotplot.R
#Author: Wanqiu Ding       
#Mail: wanqiuding@163.com
#Created Time: Wed Aug  8 10:23:46 2018
################################################

#!/bin/env Rscript

args <- commandArgs()
scriptPath = strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName = basename(scriptPath)
scriptDir = dirname(scriptPath)
args = args[-(1:5)]

parseArg = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(arg.split[2])
        }
    }
}
parseArgNum = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(as.numeric(arg.split[2]))
        }
    }
}

parseArgAsNum = function(arg, pattern, msg){
    return(parseArgNum(arg, pattern, msg))
}

usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat("-p=outputName.pdf <input.rdotplot
Option:
    Common:
    -p|pdf      FILE    The output figure in pdf[figure.pdf]
    -w|width    INT     The figure width
    -m|main     STR     The main title
    -x|xlab     STR     The xlab
    -y|ylab     STR     The ylab
    -x1         INT     The xlim start
    -x2         INT     The xlim end
    -y1         INT     The ylim start
    -y2         INT     The ylim end
    -h|help             Show help

")
	q(save = 'no')
}

mainS=20
myPdf = 'rdotplot.pdf'

if(length(args) >= 1){
    for(i in 1:length(args)){
        arg = args[i]
		if(arg == '-h' || arg == '-help') usage()
        tmp = parseArg(arg, 'p(df)?', 'p'); if(!is.null(tmp)) myPdf = tmp
        tmp = parseArgAsNum(arg, 'w(idth)?', 'w'); if(!is.null(tmp)) width = tmp
        tmp = parseArgAsNum(arg, 'x1', 'x1'); if(!is.null(tmp)) x1 = tmp
        tmp = parseArgAsNum(arg, 'x2', 'x2'); if(!is.null(tmp)) x2 = tmp
		tmp = parseArgAsNum(arg, 'y1', 'y1'); if(!is.null(tmp)) y1 = tmp
		tmp = parseArgAsNum(arg, 'y2', 'y2'); if(!is.null(tmp)) y2 = tmp
        tmp = parseArg(arg, 'm(ain)?', 'm'); if(!is.null(tmp)) main = tmp
        tmp = parseArg(arg, 'x(lab)?', 'x'); if(!is.null(tmp)) xLab = tmp
        tmp = parseArg(arg, 'y(lab)?', 'y'); if(!is.null(tmp)) yLab = tmp
		}
	}
if(exists('width')){
    pdf(myPdf, width = width)
}else{
    pdf(myPdf)
}

data = read.table(file('stdin'), header = T)
xmin=min(data[,1])
xmax=max(data[,1])
ymin=min(data[,2])
ymax=max(data[,2])
mainTitle="dotplot"
xlabel=colnames(data)[1]
ylabel=colnames(data)[2]

if(exists('x1') && exists('x2')){
	xmin=x1
	xmax=x2
	}
if(exists('y1') && exists('y2')){
    ymin=y1
	ymax=y2
	}
if(exists('main')) mainTitle=main
if(exists('xLab')) xlabel=xLab
if(exists('yLab')) ylabel=yLab
plot(x=data[,1],y=data[,2],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlabel,ylab=ylabel,main=mainTitle,pch=20)
dev.off()
