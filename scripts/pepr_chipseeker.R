#!/usr/bin/env Rscript

#############
# Performs peak annotation and plotting of PePr produced peaks
# Requirements:
# - must be PePr peaks output
# - must be from a reference genome that has an available
#   - TxDb.XXXX.UCSC.XX.knownGene
#   - org.XX.eg.db
# - Databases above must already be downloaded
#############

print_options <- function(opt) {
  opt_df <- data.frame(matrix(unlist(opt), byrow=T),stringsAsFactors=FALSE)
  row.names(opt_df)<-names(opt)
  names(opt_df) <- NULL
  message('options:')
  opt_df
}

library(optparse)

#create parser object
option_list = list(
  make_option(c('-p', '--peaks'), type = 'character', default = NULL,
              help = 'Path to PePr peaks.'),
  make_option(c('-t', '--txdb'), type = 'character', default = NULL,
              help = 'Name of TxDb to load'),
  make_option(c('-a', '--annodb'), type = 'character', default = NULL,
              help = 'Name of AnnotationDb to load'),
  make_option(c('-o', '--outDir'), type = 'character', default = 'deseq_output',
              help = 'Output directory to write all generated subdirectories.')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#check for necessary files
if (is.null(opt$peaks)) {
  print_help(opt_parser)
  stop('Path to PePr called peaks must be supplied (--peaks).', call.=FALSE)
} else if (is.null(opt$txdb)) {
  print_help(opt_parser)
  stop('TxDb must be supplied (--txdb).', call.=FALSE)
} else if (is.null(opt$annodb)) {
  print_help(opt_parser)
  stop('AnnotationDb must be supplied (--annodb).', call.=FALSE)
} else if (is.null(opt$outDir)) {
  print_help(opt_parser)
  stop('Output directory to write all generated plots and files (--outDir).', call.=FALSE)
} else
  message('Inputs detected: peaks, txdb, annodb, outDir Continue...')

print_options(opt)

#load libs
message('loading libraries begins')
suppressPackageStartupMessages(library('ChIPseeker', character.only=TRUE))
suppressPackageStartupMessages(library('clusterProfiler', character.only=TRUE))
suppressPackageStartupMessages(library('data.table', character.only=TRUE))

#test options
#opt$peaks <- '/Users/cjsifuen/ActiveProjects/LSA_Denver_rdenver_CS3_cjsifuen_HI-2582/analysis_01_25/pepr/sharp/peaks/klf13__PePr_peaks_fixed_chr.bed'
#opt$txdb <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#opt$annodb <- 'org.Mm.eg.db'
#opt$outDir <- '/Users/cjsifuen/ActiveProjects/LSA_Denver_rdenver_CS3_cjsifuen_HI-2582/analysis_01_25/pepr/sharp/peaks/'

basefile <- basename(tools::file_path_sans_ext(opt$peaks))

#### need to account for possibility that an uninstalled database is named
if (as.character(opt$txdb) %in% rownames(installed.packages())) { 
  suppressPackageStartupMessages(library(opt$txdb, character.only = TRUE))
} else { #### stop
  stop(paste0(opt$txdb, ' package was not found. Please install TxDb database package for species of interest.'))
}

if (as.character(opt$annodb) %in% rownames(installed.packages())) { 
  suppressPackageStartupMessages(library(opt$annodb, character.only = TRUE))
} else { #### stop
  stop(paste0(opt$annodb, ' package was not found. Please install AnnotationDb database package for species of interest.'))
}

#####
# Get peaks and txdb
####

#get txdb
assign(x = 'txdb', value = get(opt$txdb))

#read in peaks
peak <- readPeakFile(peakfile = opt$peaks)


#####
# Coverage plots
#####

cat('Calculating and plotting coverage...', fill = T)
pdf(file = paste(opt$outDir,'CoverageByChrom.pdf', sep = '/'), onefile = TRUE, bg = 'white')
covplot(peak = peak, 
        weightCol = 'V5')
dev.off()

#####
# TSS related plots
#####

cat('Beginning TSS plots...', fill = T)
#heatmap peaks at tss regions
pdf(file = paste(opt$outDir,'PeaksHeatmapTSS.pdf', sep = '/'), onefile = TRUE)
peakHeatmap(peak = peak, 
            weightCol = 'V5', 
            TxDb = txdb, 
            upstream = 3000, 
            downstream = 3000,
            xlab = "Genomic Region (5'->3')", 
            color = 'red')
dev.off()

#profile peaks at tss regions
pdf(file = paste(opt$outDir,'PeaksAvgProfileTSS.pdf', sep = '/'), onefile = TRUE)
plotAvgProf2(peak = peak, 
             TxDb=txdb, 
             upstream=3000, 
             weightCol = 'V5',
             downstream=3000,
             conf = 0.95,
             xlab="Genomic Region (5'->3')",
             ylab = "Read Count Frequency")
dev.off()

#####
# Genebody related plots
#####

cat('Getting genebody regions...', fill = T)

#profile peaks at gene
gene <- getBioRegion(TxDb = txdb, upstream = 3000, downstream = 3000, by = 'gene')
geneTagMatrix <- getTagMatrix(peak = peak, 
                              weightCol = 'V5',
                              windows = gene)

cat('Beginning genebody plots...', fill = T)

#heatmap peaks at genebody regions
pdf(file = paste(opt$outDir,'PeaksHeatmapGenebody.pdf', sep = '/'), onefile = TRUE)
tagHeatmap(tagMatrix = geneTagMatrix, 
           color='red',
           xlab = "Genomic Region (5'->3')", 
           xlim = c(-3000,3000))
dev.off()

#profile peaks at genebody regions
pdf(file = paste(opt$outDir,'PeaksAvgProfileGenebody.pdf', sep = '/'), onefile = TRUE)
plotAvgProf(tagMatrix = geneTagMatrix, 
                                xlim=c(-3000, 3000), 
                                xlab="Genomic Region (5'->3')", 
                                ylab = "Read Count Frequency", 
                                conf = 0.95)
dev.off()


#####
# Exon related plots
#####

cat('Getting exon regions...', fill = T)
#profile peaks at exon
exon <- getBioRegion(TxDb = txdb, upstream = 1000, downstream = 1000, by = 'exon')
exonTagMatrix <- getTagMatrix(peak = peak,
                              weightCol = 'V5',
                              windows = exon)

cat('Beginning exon plots...', fill = T)
#heatmap peaks at exon regions
pdf(file = paste(opt$outDir,'PeaksHeatmapExon.pdf', sep = '/'), onefile = TRUE)
tagHeatmap(tagMatrix = exonTagMatrix, 
           xlab = "Genomic Region (5'->3')", 
           color='red', 
           xlim = c(-1000,1000))
dev.off()

#profile peaks at exon regions
pdf(file = paste(opt$outDir,'PeaksAvgProfileExon.pdf', sep = '/'), onefile = TRUE)
plotAvgProf(tagMatrix = exonTagMatrix, 
            xlim=c(-1000, 1000), 
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency", 
            conf = 0.95)
dev.off()

#####
# Intron related plots
#####

cat('Getting inton regions...', fill = T)
#profile peaks at intron
intron <- getBioRegion(TxDb = txdb, upstream = 1000, downstream = 1000, by = 'intron')
intronTagMatrix <- getTagMatrix(peak = peak, 
                               weightCol = 'V5',
                               windows = intron)

cat('Beginning intron plots...', fill = T)
#heatmap peaks at intron regions
pdf(file = paste(opt$outDir,'PeaksHeatmapIntron.pdf', sep = '/'), onefile = TRUE)
tagHeatmap(tagMatrix = intronTagMatrix, 
           xlab = "Genomic Region (5'->3')", 
           color='red', 
           xlim = c(-1000,1000))
dev.off()

#profile peaks at intron regions
pdf(file = paste(opt$outDir,'PeaksAvgProfileIntron.pdf', sep = '/'), onefile = TRUE)
plotAvgProf(tagMatrix = intronTagMatrix, 
            xlim=c(-1000, 1000), 
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency", 
            conf = 0.95)
dev.off()




#####
# Peak annotation and plots
#####

cat('Annotating...', fill = T)

#annotate peaks
peakAnno <- annotatePeak(peak = peak, 
                         tssRegion = c(-3000,3000), 
                         TxDb = txdb, 
                         annoDb = opt$annodb)
peakAnnoDf <- as.data.frame(x = peakAnno)

cat('Writing out annotation...', fill = T)
write.table(x = peakAnnoDf, 
            file = paste(opt$outDir, '/' ,basefile, '_annot.bed', sep = ''), 
            quote = F, 
            sep = '\t', 
            row.names = F, 
            col.names = T)


cat('Annotation plots...', fill = T)
pdf(file = paste(opt$outDir,'PeaksDistByFeature.pdf', sep = '/'), onefile = TRUE, bg = 'white', paper = 'USr')
plotAnnoPie(x = peakAnno)
plotAnnoBar(x = peakAnno)
vennpie(x = peakAnno, r = 0.08)
upsetplot(x = peakAnno)
upsetplot(x = peakAnno, vennpie = T)
dev.off()

pdf(file = paste(opt$outDir,'PeaksByDistanceToTSS.pdf', sep = '/'), onefile = TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

cat('Finished.', fill = T)
