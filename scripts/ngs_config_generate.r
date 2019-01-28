#!/usr/bin/env Rscript

# R scripts for generating config file for ngsplot
# Author @chuan-wang https://github.com/chuan-wang

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)
dir <-args[1]
datafiles <-args[-1]
# Check input args
if (length(args) < 2) {
  stop("Usage: ngs_config_generate.r [outdir] [input-1] [input-2]..", call.=FALSE)
}

r <- length(datafiles)
r <- r-1
table<-matrix(0, nrow=length(datafiles), ncol=3)
table<-as.data.frame(table)
for (i in 1:length(datafiles)){
   table[i,1]<-datafiles[i]
   table[i,2]<-(-1)
   table[i,3]<-paste('"', gsub(".filtered.sorted.bam", "", as.character(basename(datafiles[i]))), '"', sep="")
}
table<-table[order(table[,1]),]
write.table(table, file=paste0(dir,"/ngsplot_config"),sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
