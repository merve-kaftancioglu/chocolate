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
  make_option(c('-s', '--species'), type = 'character', default = NULL,
              help = 'Species.'),
  make_option(c('-o', '--outDir'), type = 'character', default = 'deseq_output',
              help = 'Output directory to write all generated subdirectories.')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#check for necessary files
if (is.null(opt$peaks)) {
  print_help(opt_parser)
  stop('Path to PePr called peaks must be supplied (--peaks).', call.=FALSE)
} else if (is.null(opt$species)) {
  print_help(opt_parser)
  stop('Species must be supplied (--species).', call.=FALSE)
} else if (is.null(opt$outDir)) {
  print_help(opt_parser)
  stop('Output directory to write all generated plots and files (--outDir).', call.=FALSE)
} else
  message('Inputs detected: peaks, species, outDir Continue...')

print_options(opt)

#load libs
message('loading libraries begins')
suppressPackageStartupMessages(library('ChIPseeker', character.only=TRUE))
suppressPackageStartupMessages(library('clusterProfiler', character.only=TRUE))
suppressPackageStartupMessages(library('data.table', character.only=TRUE))
suppressPackageStartupMessages(library('BiocParallel', character.only=TRUE))


#load database for organism
ParseSpeciesInfo <- function(speciesInfoFile, species) {
  df <- read.table(file=speciesInfoFile, header = TRUE, sep = '\t', stringsAsFactors = FALSE, as.is = TRUE, quote = '') #read in file
  kegg.db <- df$kegg[df$species == species] #get kegg database name for species
  go.db <- df$go[df$species == species] #get go database name for species
  react.db <- df$reactome[df$species == species] #get reactome database name for species
  return (list(kegg.db, go.db, react.db))
}

db.data <- ParseSpeciesInfo(speciesInfoFile = "~/Documents/Core/watermelon/functionalAnnotation/species_database_names.txt", species = species) 
go.db <- unlist(db.data[2]) #data[2] = go database

#### need to account for possibility that an unspecified species is named
if (as.character(go.db) %in% rownames(installed.packages())) { 
  suppressPackageStartupMessages(library(go.db, character.only = TRUE))
} else { #### specify stop
  stop(paste0(go.db, ' package was not found. Please install GO database package for species of interest.'))
}


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


#read in peaks
peak <- readPeakFile(peakfile = '/Users/cjsifuen/ActiveProjects/LSA_Denver_rdenver_CS3_cjsifuen_HI-2582/analysis_01_25/pepr/sharp/peaks/klf13__PePr_peaks_fixed_chr.bed')

#coverage plots
peak_covplot <- covplot(peak = peak, 
                 weightCol = 'V5')

#heatmap peaks at tss regions
#promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
#tagMatrix <- getTagMatrix(peak, windows = promoter)
#tagHeatmap(tagMatrix = tagMatrix, color='red', xlim = c(-3000,3000))
peak_tssHeatmap <- peakHeatmap(peak = peak, 
                               weightCol = 'V5', 
                               TxDb = txdb, 
                               upstream = 3000, 
                               downstream = 3000, 
                               color = 'red')

#profile peaks at tss regions
peak_tssProfile <- plotAvgProf2(peak = peak, 
                                TxDb=txdb, 
                                upstream=3000, 
                                weightCol = 'V5',
                                downstream=3000,
                                conf = 0.95,
                                xlab="Genomic Region (5'->3')", 
                                ylab = "Read Count Frequency")

#profile peaks at gene, intron, exon
gene <- getBioRegion(TxDb = txdb, upstream = 3000, downstream = 3000, by = 'gene')
exon <- getBioRegion(TxDb = txdb, upstream = 1000, downstream = 1000, by = 'exon')
intron <- getBioRegion(TxDb = txdb, upstream = 1000, downstream = 1000, by = 'intron')

geneTagMatrix <- getTagMatrix(peak = peak, 
                              weightCol = 'V5',
                              windows = gene)
exonTagMatrix <- getTagMatrix(peak = peak,
                              weightCol = 'V5',
                              windows = exon)
intronTagMatrix <- getTagMatrix(peak = peak, 
                               weightCol = 'V5',
                               windows = intron)
#profile peaks Heatmaps 
peak_geneHeatmap <- tagHeatmap(tagMatrix = geneTagMatrix, 
                               color='red', 
                               xlim = c(-3000,3000))

peak_exonHeatmap <- tagHeatmap(tagMatrix = exonTagMatrix, 
                               color='red', 
                               xlim = c(-1000,1000))

peak_intronHeatmap <- tagHeatmap(tagMatrix = intronTagMatrix, 
                               color='red', 
                               xlim = c(-1000,1000))

#profile peaks plot
peak_geneProfile <- plotAvgProf(tagMatrix = geneTagMatrix, 
                                xlim=c(-3000, 3000), 
                                xlab="Genomic Region (5'->3')", 
                                ylab = "Read Count Frequency", 
                                conf = 0.95)
peak_exonProfile <- plotAvgProf(tagMatrix = exonTagMatrix, 
                                xlim=c(-1000, 1000), 
                                xlab="Genomic Region (5'->3')", 
                                ylab = "Read Count Frequency", 
                                conf = 0.95)
peak_intronProfile <- plotAvgProf(tagMatrix = intronTagMatrix, 
                                  xlim=c(-1000, 1000), 
                                  xlab="Genomic Region (5'->3')", 
                                  ylab = "Read Count Frequency", 
                                  conf = 0.95)

#annotate peaks
peakAnno <- annotatePeak(peak = peak, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = 'org.Mm.eg.db')
peakAnnoDf <- as.data.frame(x = peakAnno)


#plot annotation pie
plotAnnoPie(x = peakAnno)
plotAnnoBar(x = peakAnno)
vennpie(x = peakAnno)
upsetplot(x = peakAnno)
upsetplot(x = peakAnno, vennpie = T)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
