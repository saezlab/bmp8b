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

standard <- 'none'
treatments <- c('NE', 'BMP8b', 'BMP8b_NE')
datapath <- 'src/'
gmtpath <- 'piano/'
setsfpost <- 'v5.2.symbols.gmt'
setspref <- c('c2.cp', 'c2.cp.kegg', 'c5.bp', 'h.all', 'c2.cp.reactome', 'c2.cp.biocarta')
setspref <- c('c2.cp', 'c5.bp', 'h.all')
mapfile  <- 'uniprot_genesymbol.tab'
result <- list()
filespath <- ''

###

result <- list()

genesets <- list()
for(setpref in setspref){
    genesets <- c(genesets,
                    unlist(
                        getGmt(con = paste(gmtpath, 'c2.cp.', setsfpost, sep = ''))
                    )
                )
}

gene_to_term <- NULL

for(geneset in genesets){
    
    temp1 <- geneIds(geneset)
    temp2 <- setName(geneset)
    temp3 <- as.data.frame(cbind(temp1, rep(temp2, length(temp1))))
    names(temp3) <- c('genesymbol', 'term')
    gene_to_term <- rbind(gene_to_term, temp3)
    
}

unip_to_gsym <- read.table(mapfile, header = FALSE, sep = '\t')
colnames(unip_to_gsym) <- c('uniprot', 'genesymbol')
gene_to_term <- merge(gene_to_term, unip_to_gsym, all = FALSE)

geneSet <- loadGSC(gene_to_term[,c(3, 2)])

for(tr in treatments){
    
    cat(paste('Working on', tr, setpref, '\n'))
    
    infile <- paste(filespath, 'fctop_uniqp_', tr, '_', standard, '.csv', sep = '')
    data <- read.table(infile, header = TRUE, sep = '\t')
    
    rownames(data) <- data$uniprot
    
    logfc <- data['logfc']
    pval  <- data['pval']
    tval  <- data['tval']
    
    gsaRes1 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'mean')
    gsaRes2 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'median')
    gsaRes3 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'sum')
    gsaRes4 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'maxmean')
    gsaRes5 <- runGSA(geneLevelStats = pval, directions = logfc, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'reporter')
    gsaRes6 <- runGSA(geneLevelStats = pval, directions = logfc, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'tailStrength')
    gsaRes7 <- runGSA(geneLevelStats = pval, directions = logfc, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'wilcoxon')
    gsaRes8 <- runGSA(tval, gsc = geneSet, adjMethod = 'fdr', geneSetStat = 'page')
    
    resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6, gsaRes7, gsaRes8)
    names(resList) <- c('mean', 'median', 'sum', 'maxmean', 'reporter', 'tailStrength', 'wilcoxon', 'page')
    
    result[[tr]] <- resList
    
}

for(tr in treatments){
    
    cairo_pdf(filename = paste('pathways_heatmap_', 'all_', tr, '.pdf', sep = ''),
                onefile = TRUE, width = 12, height = 24, pointsize = 12)
    
        ch <- consensusHeatmap(result[[tr]], cutoff = 100, method = 'median',
                                ncharLabel = 60, cellnote = 'consensusScore',
                                cex = 1, colorkey = FALSE)
    
    dev.off()
    
    consensus_cp <- ch$pMat
    
    write.csv(consensus_cp, paste('consensus_cp_', 'all_', tr, '.csv', sep = ''))
    
}

save.image(file = 'moregenesets2.Rdata')

####

cairo_pdf(filename = paste('pathways_heatmap_', setpref, '_', tr, '.pdf', sep = ''),
          onefile = TRUE, width = 12, height = 24, pointsize = 12)

ch <- consensusHeatmap(result[[setpref]][[tr]], cutoff = 100, method = 'median',
                        ncharLabel = 60, cellnote = 'consensusScore',
                        cex = 1, colorkey = FALSE, plot = FALSE)

dev.off()

my_palette <- colorRampPalette(c("red", "black"))(n = 299)

cairo_pdf(filename = paste('pathways_heatmap_g_', setpref, '_', tr, '.pdf', sep = ''),
          onefile = TRUE, width = 12, height = 24, pointsize = 12)

par(mgp=c(0, 0, 0))
heatmap.2(ch$rankMat,
          Colv = "NA",
          dendrogram = 'row',
          density.info = 'none',
          key = FALSE,
          col = my_palette,
          margins =c(15,15),
          notecol = '#CCCCCC',
          cellnote = ch$rankMat,
          sepwidth = c(0.01/dim(ch$rankMat)[1], 0.01/dim(ch$rankMat)[2]),
          colsep = seq(dim(ch$rankMat)[2]),
          rowsep = seq(dim(ch$rankMat)[1]),
          sepcolor = 'white',
          trace = 'none')

dev.off()
