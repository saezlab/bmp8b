#!/usr/bin/Rscript

#
# # (c) Dénes Türei EMBL-EBI 2016-2017
# # turei.denes@gmail.com
# # thanks to Aurelien Dugourd
#

library(ggplot2)
library(gskb)
library(piano)
library(BioNet)
library(igraph)
library(reshape)
library(data.table)
library(pheatmap)
library(limma)
library(grid)
library(gridExtra)
library(GSEABase)

require(gplots)
require(RColorBrewer)
require(data.table)

standard <- 'none'
filespath <- ''
datapath <- 'src/'
gmtpath <- 'piano/'
setsfpost <- 'v5.2.symbols.gmt'
mapfile  <- 'uniprot_genesymbol.tab'

geneset.df <- function(geneset){
    return(as.data.frame(cbind(geneIds(geneset), setName(geneset))))
}

load_genesets <- function(setspref){
    
    genesets <- list()
    
    for(setpref in setspref){
        genesets <- c(genesets,
                        unlist(
                            getGmt(con = paste(gmtpath, setpref, '.', setsfpost, sep = ''))
                        )
                    )
    }
    
    gene_to_term <- NULL
    
    gene_to_term <- rbindlist(lapply(genesets, geneset.df))
    
    setnames(gene_to_term, c('genesymbol', 'term'))
    
    unip_to_gsym <- read.table(mapfile, header = FALSE, sep = '\t')
    colnames(unip_to_gsym) <- c('uniprot', 'genesymbol')
    gene_to_term <- merge(gene_to_term, unip_to_gsym, all = FALSE, by = 'genesymbol')
    
    geneSet <- loadGSC(gene_to_term[,c(3, 2), with = FALSE])
    
    return(geneSet)
    
}

run_piano <- function(infile, geneSet, signed = FALSE){
    
    data <- read.table(infile, header = TRUE, sep = '\t')
    
    rownames(data) <- data$uniprot
    
    if(signed){
        logfc <- data['logfc'] * data['effect']
    }else{
        logfc <- data['logfc']
    }
    
    pval  <- data['pval']
    tval  <- data['tval']
    
    maxgssz <- dim(data)[1] / 3
    
    gsaRes1 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'mean', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes2 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'median', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes3 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'sum', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes4 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'maxmean', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes5 <- runGSA(geneLevelStats = pval, directions = logfc,
                      gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'reporter', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes6 <- runGSA(geneLevelStats = pval, directions = logfc,
                      gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'tailStrength', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes7 <- runGSA(geneLevelStats = pval, directions = logfc,
                      gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'wilcoxon', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    gsaRes8 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr',
                      geneSetStat = 'page', ncpus = 4,
                      gsSizeLim = c(1, maxgssz))
    
    resList <- list(gsaRes1, gsaRes2, gsaRes3,
                    gsaRes4, gsaRes5, gsaRes6,
                    gsaRes7, gsaRes8)
    names(resList) <- c('mean', 'median', 'sum',
                        'maxmean', 'reporter',
                        'tailStrength', 'wilcoxon',
                        'page')
    
    return(resList)
}

run_piano_series <- function(treatments, geneSet, signed){
    
    result <- list()
    
    strsign <- ifelse(signed, 'sign_', '')
    
    for(tr in treatments){
        
        cat(paste(':::: Working on', tr, '\n'))
        
        infile <- paste(filespath, 'fctop_uniqp_', strsign, tr, '_', standard, '.csv', sep = '')
        
        result[[tr]] <- run_piano(infile, geneSet)
        
    }
    
    return(result)
    
}
