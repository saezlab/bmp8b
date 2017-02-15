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
require(snowfall)

source('heatmap.2.R')

standard <- 'none'
filespath <- ''
datapath <- 'src/'
gmtpath <- 'piano/'
setsfpost <- 'v5.2.symbols.gmt'
mapfile  <- 'uniprot_genesymbol.tab'
filepref <- 'fctop_uniqp_'
filepost <- '.csv'
plotpref <- 'pathway_heatmap_'

my_palette <- colorRampPalette(c('#FFFFFF', '#990000'))(n = 299)

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
        
        cat(paste('\n:::: Working on', tr, '\n\n'))
        
        infile <- paste(filespath, filepref, strsign, tr, '_',
                        standard, filepost, sep = '')
        
        result[[tr]] <- run_piano(infile, geneSet, signed = signed)
        
    }
    
    return(result)
    
}

merge_args <- function(defaults, final){
    
    for(arg in names(defaults)){
        
        if(is.null(final[[arg]]) & !(arg %in% names(final))){
            
            final[[arg]] <- defaults[[arg]]
            
        }
        
    }
    
    return(final)
    
}

piano_plot <- function(result, tr, plotname, cols,
                       dir_col = NULL,
                       consensus_args = list(),
                       heatmap_args = list(),
                       cairo_args = list(),
                       maintitle = ''){
    
    maintitle <- ifelse(is.null(dir_col), paste('Geneset ranks:', tr), maintitle)
    
    consensus_defaults <- list(
        cutoff = 100,
        method = 'median',
        ncharLabel = 60,
        cellnote = 'consensusScore',
        cex = 1,
        colorkey = FALSE,
        plot = FALSE
    )
    
    consensus_args <- merge_args(consensus_defaults, consensus_args)
    
    if(!is.null(dir_col)){
        cutoff <- consensus_args[['cutoff']]
        cutoff <- ifelse(cutoff, cutoff, 50)
        consensus_args[['cutoff']] <- 1000000
    }
    
    if(!is.null(dir_col)){
        ch <- data.frame()
        for(tr0 in names(result)){
            ch0 <- do.call('consensusHeatmap', c(list(resList = result[[tr0]]),
                                                 consensus_args))
            ch0 <- ch0$rankMat[,c(dir_col)]
            ch <- transform(merge(ch, ch0, by = 0, all = TRUE),
                            row.names = Row.names, Row.names = NULL)
        }
        fullmat <- as.matrix(ch)
        colnames(fullmat) <- names(result)
        plotname2 <- as.character(dir_col)
        mat <- fullmat[apply(fullmat, MARGIN = 1, function(x) any(x < cutoff)), ]
    }else{
        
        ch <- do.call('consensusHeatmap', c(list(resList = result[[tr]]),
                                            consensus_args))
        mat <- ch$rankMat[,cols]
        plotname2 <- tr
    }
    
    rownames(mat) <- sapply(rownames(mat), function(x){substr(x, 1, 40)})
    
    heatmap_defaults <- list(
        x = mat,
        Colv = "NA",
        dendrogram = 'row',
        density.info = 'none',
        key = FALSE,
        col = my_palette,
        margins = c(12,30),
        notecol = '#000000',
        cellnote = mat,
        sepwidth = c(0.01/dim(mat)[1], 0.01/dim(mat)[2]),
        colsep = seq(dim(mat)[2]),
        rowsep = seq(dim(mat)[1]),
        sepcolor = 'white',
        trace = 'none',
        #lmat = matrix(c(2, 1), nrow = 1),
        lhei = c(0.03, 1.0),
        xlab = 'Treatments',
        ylab = 'Genesets',
        main = maintitle,
        cexCol = 1.8,
        cexRow = 1.2,
        cex.lab = 2.4,
        adjCol = c(NA, 0.5)
    )
    
    cairo_defaults <- list(
        filename = paste(plotpref, plotname, '_', plotname2, '.pdf', sep = ''),
        onefile = TRUE,
        width = 9,
        height = 36,
        pointsize = 12
    )
    
    cat(paste(':::: Plotting to `', cairo_defaults$filename, '`.\n'))
    
    
    heatmap_args <- merge_args(heatmap_defaults, heatmap_args)
    cairo_args   <- merge_args(cairo_defaults,   cairo_args)
    
    do.call('cairo_pdf', c(cairo_args))
    
        #par(mgp=c(0, 0, 0))
        par(oma = c(2, 2, 2, 2))
        do.call('heatmap.2', c(heatmap_args))
    
    dev.off()
    
    return(list(full = fullmat, shown = mat))
    
}

piano_plot_series <- function(result, plotname,
                              cols = c(1, 2, 3, 4, 5),
                              dir_col = NULL,
                              consensus_args = list(),
                              heatmap_args = list(),
                              cairo_args = list(),
                              maintitle = '',
                              return_result = FALSE){
    
    if(is.null(dir_col)){
        
        result <- list()
        
        for(tr in names(result)){
            
            result[[tr]] <- (
                piano_plot(result, tr, plotname, cols,
                        consensus_args, heatmap_args, cairo_args)
            )
            
            if(return_result){
                return(result)
            }
        }
    }else{
        result <- (
            piano_plot(result, names(result)[0], plotname, cols, dir_col,
                    consensus_args, heatmap_args, cairo_args,
                    maintitle)
        )
        
    }
    
    if(return_result){
        return(result)
    }
    
}
