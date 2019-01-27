#!/usr/bin/env Rscript

# R scripts for generating config file for ngsplot
# Author @chuan-wang https://github.com/chuan-wang

# Command line argument processing
datafiles <- commandArgs(trailingOnly=TRUE)

dir <-args[1]

# Check input args
if (length(datafiles) < 2) {
  stop("Usage: ngs_config_generate.r [outdir] [input-1] [input-2]..", call.=FALSE)
}

table<-matrix(0, nrow=(length(datafiles)-1), ncol=3)
table<-as.data.frame(table)
for (i in 2:length(datafiles)){
   table[i,1]<-datafiles[i]
   table[i,2]<-(-1)
   table[i,3]<-paste('"', gsub(".filtered.sorted.bam", "", as.character(datafiles[i])), '"', sep="")
}
table<-table[order(table[,1]),]
write.table(table, file=paste0(dir,"/ngsplot_config"),sep='\t', quote=FALSE, row.names=FALSE,
col.names=FALSE)
