#!/usr/bin/env Rscript

library(clusterProfiler)
library(ReactomePA)
library(data.table)



speciesInfoFile <- "~/Documents/Core/watermelon/functionalAnnotation/species_database_names.txt"
#opt$species <- 'Human'
species <- opt$species

MSigDBInfoFile <- "~/Documents/Core/watermelon/functionalAnnotation/mbsigdb_names.txt"

#function to give database information for msigdb
ParseMSigDBInfo <- function(MSigDBInfoFile) {
  df <- read.table(file=MSigDBInfoFile, header = TRUE, sep = '\t', stringsAsFactors = FALSE, as.is = TRUE, quote = '') #read in file
  GMT_FILES <- df$name
  names(GMT_FILES) <- df$file
  return (GMT_FILES)
}

#function to give species database based on species name and either go or kegg
ParseSpeciesInfo <- function(speciesInfoFile, species) {
  df <- read.table(file=speciesInfoFile, header = TRUE, sep = '\t', stringsAsFactors = FALSE, as.is = TRUE, quote = '') #read in file
  kegg.db <- df$kegg[df$species == species] #get kegg database name for species
  go.db <- df$go[df$species == species] #get go database name for species
  react.db <- df$reactome[df$species == species] #get reactome database name for species
  return (list(kegg.db, go.db, react.db))
}

#writer
ResultWriter <- function(result.non.readable, file.suffix, ont.type=NULL, column.of.interest=NULL, simplify.result=FALSE, convert.id=FALSE, go.db=NULL,kegg.db=NULL, key.type) {
  if (!is.na(result.non.readable[1]$ID)){#if over-represented terms are present
    
    if (simplify.result == TRUE) { #simplify = true; convert, simplify, merge, write
      result.readable <- setReadable(x = result.non.readable, OrgDb = go.db, keytype = key.type) 
      result.non.readable.simplified <- simplify(result.non.readable, cutoff=0.7, by = 'p.adjust', select_fun=min)
      result.readable.simplified <- simplify(result.readable, cutoff=0.7, by = 'p.adjust', select_fun=min)
      result.merged <- merge(x = as.data.frame(result.non.readable), y = as.data.frame(result.readable)[,c(1,column.of.interest)], by = 'ID', all.x = TRUE)
      result.simplified.merged <- merge(x = as.data.frame(result.non.readable.simplified), y = as.data.frame(result.readable.simplified)[,c(1,column.of.interest)], by = 'ID', all.x = TRUE) 
      write.table(x = result.merged, file = paste0(out.dir, '/',base.file, file.suffix, ont.type, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
      write.table(x = result.simplified.merged, file = paste0(out.dir, '/',base.file, file.suffix,'simplified_' , ont.type, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
      
    } else if (simplify.result == FALSE & convert.id == TRUE){ #simplify = false & convert = true; convert, merge, write
      result.readable <- setReadable(x = result.non.readable, OrgDb = go.db, keytype = key.type)
      result.merged <- merge(x = as.data.frame(result.non.readable), y = as.data.frame(result.readable)[,c(1,column.of.interest)], by = 'ID', all.x = TRUE) 
      write.table(x = result.merged, file = paste0(out.dir, '/',base.file, file.suffix, ont.type, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
      
    } else {
      write.table(x = as.data.frame(result.non.readable), file = paste0(out.dir, '/',base.file, file.suffix, ont.type, '.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
    }
    
  } else {
    message('Skipping output for ', base.file, file.suffix, ont.type)
  }
}

#general KEGG function to call specific functions
KeggAnalysis <- function(filtered.geneset, full.geneset, db.name, pvalue, analysis.type) {
  if (analysis.type == 'enrichKEGG' | analysis.type == 'enrichMKEGG') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          organism      = db.name,
                                          keyType       = 'kegg',
                                          pAdjustMethod = 'BH',
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2)
    return (analysis.result)
  } else if (analysis.type == 'gseKEGG' | analysis.type == 'gseMKEGG') {
    analysis.result <- get(analysis.type)(geneList     = full.geneset,
                                          organism     = db.name,
                                          keyType      = 'kegg',
                                          nPerm        = 1000,
                                          minGSSize    = 10,
                                          maxGSSize    = 300,
                                          pvalueCutoff = pvalue,
                                          verbose      = FALSE)
      return (analysis.result)
  } else {
    stop('Analysis type must be "enrichKEGG", "enrichMKEGG", "gseKEGG", or "gseMKEGG".', call.=FALSE)
  }
} 


#general GO function to call specific functions
GoAnalysis <- function(filtered.geneset, full.geneset, db.name, ont.type, pvalue, analysis.type){
  if (analysis.type == 'groupGO') {
    analysis.result <- get(analysis.type)(gene     = names(filtered.geneset),
                                          OrgDb    = db.name,
                                          ont      = ont.type,
                                          keytype  = 'ENTREZID',
                                          level    = 3,
                                          readable = FALSE)
    return (analysis.result)
  } else if (analysis.type == 'enrichGO') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          keytype       = 'ENTREZID',
                                          OrgDb         = db.name,
                                          ont           = ont.type,
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          pAdjustMethod = 'BH', 
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2,
                                          readable      = FALSE)
    return (analysis.result) 
  } else if (analysis.type == 'gseGO'){
    analysis.result <- get(analysis.type)(geneList      = full.geneset,
                                          OrgDb         = db.name,
                                          ont           = ont.type,
                                          keytype       = 'ENTREZID',
                                          nPerm         = 1000,
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          pAdjustMethod = 'BH',
                                          pvalueCutoff  = pvalue,
                                          verbose       = FALSE)
    return (analysis.result)
  } else {
    stop('Analysis type must be "groupGO", "enrichGO", or "gseGO".', call.=FALSE)
  }
} 

##general DO function to call specific functions
DoAnalysis <- function(filtered.geneset, full.geneset, db.name, pvalue, analysis.type){
  if (analysis.type == 'enrichDO') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          ont           = 'DO',
                                          pAdjustMethod = 'BH',
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2,
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          readable      = FALSE)
    return (analysis.result)
  } else if (analysis.type == 'enrichNCG' | analysis.type == 'enrichDGN') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          pAdjustMethod = 'BH',
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2,
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          readable      = FALSE)
    return (analysis.result)
  } else if (analysis.type == 'gseDO'){
    analysis.result <- get(analysis.type)(geneList      = full.geneset,
                                          nPerm         = 1000,
                                          pAdjustMethod = 'BH',
                                          pvalueCutoff  = pvalue,
                                          minGSSize     = 10,
                                          maxGSSize     = 300)
    return (analysis.result)
  } else if (analysis.type == 'gseNCG' | analysis.type == 'gseDGN') {
    analysis.result <- get(analysis.type)(geneList      = full.geneset,
                                          nPerm         = 1000,
                                          pAdjustMethod = 'BH',
                                          pvalueCutoff  = pvalue,
                                          minGSSize     = 10,
                                          maxGSSize     = 300)
    return(analysis.result)
  } else {
    stop('Analysis type must be "enrichDO", "enrichNCG", "enrichDGN", "gseDO", "gseNCG", or "gseDGN".', call.=FALSE)
  }
}

#general MBSigDB function to call specific functions
MbSigAnalysis <- function(filtered.geneset, full.geneset, db.name, pvalue, gmt.name, analysis.type) {
  if (analysis.type == 'enricher') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          TERM2GENE     = gmt.name,
                                          pAdjustMethod = 'BH',
                                          minGSSize     = 5,
                                          maxGSSize     = 300,
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2)
    return (analysis.result)
  } else if (analysis.type == 'GSEA') {
    analysis.result <- get(analysis.type)(geneList = full.geneset,
                                          TERM2GENE = gmt.name,
                                          pvalueCutoff = pvalue, 
                                          nPerm = 1000, 
                                          pAdjustMethod = 'BH',
                                          minGSSize = 5, 
                                          maxGSSize = 300)
      return (analysis.result)
    } else {
      stop('Analysis type must be "enricher" or "GSEA".', call.=FALSE)
  }
} 

#general Reactome function to call specific functions
ReactomeAnalysis <- function(filtered.geneset, full.geneset, species, pvalue, analysis.type) {
  if (analysis.type == 'enrichPathway') {
    analysis.result <- get(analysis.type)(gene          = names(filtered.geneset),
                                          universe      = names(full.geneset),
                                          organism      = species,
                                          pAdjustMethod = 'BH',
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          pvalueCutoff  = pvalue,
                                          qvalueCutoff  = 0.2,
                                          readable      = FALSE)
    return (analysis.result)
  } else if (analysis.type == 'gsePathway') {
    analysis.result <- get(analysis.type)(geneList      = full.geneset,
                                          organism      = species,
                                          nPerm         = 1000,
                                          minGSSize     = 10,
                                          maxGSSize     = 300,
                                          pvalueCutoff  = pvalue,
                                          pAdjustMethod = 'BH',
                                          verbose       = FALSE)
    return (analysis.result)
  } else {
    stop('Analysis type must be "enrichPathway", or "gsePathway".', call.=FALSE)
  }
} 






