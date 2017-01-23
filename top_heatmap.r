#!/usr/bin/Rscript

#
# # Dénes Türei EMBL-EBI 2017
#

require(ggplot2)
require(reshape2)
require(ggdendro)

inFile <- 'fctop_none.csv'
fctop <- read.table(inFile, header = TRUE, sep = '\t')

fctop <- fctop[, names(fctop) != 'resnum']

fctop.m <- melt(fctop)

fctop.m$label <- factor(fctop$label, levels = unique(as.character(fctop$label)))

fctop.m['logfc'] <- sign(fctop.m$value) * log2(abs(fctop.m$value)) / max(log2(abs(fctop.m$value)))

p <- ggplot(fctop.m, aes(variable, label)) +
    geom_tile(aes(fill = logfc)) +
    scale_fill_gradient2(low = 'orangered3', mid = 'white', high = 'steelblue') +
    xlab('Treatment') +
    ylab('Phosphorylation sites') +
    theme(
        axis.text.y = element_text(size = 4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank()
    )

ggsave('gg_fcdiff_top_heatmap.pdf', device = cairo_pdf, width = 4, height = 48)

# Only top 50:

fctop50 <- head(fctop, 111)
rnames <- fctop50$psite
rownames(fctop50) <- rnames
mfctop50 <- as.matrix(fctop50[,6:8])
rownames(mfctop50) <- rnames
psites_dendro <- as.dendrogram(hclust(dist(mfctop50)))
ordr <- order.dendrogram(psites_dendro)
ofctop50 <- fctop50[ordr,]
onames <- attr(ofctop50, 'dimnames')

fctop.m <- melt(fctop50)

fctop.m$psite <- factor(fctop50$psite, levels = unique(as.character(fctop50$psite)))

fctop.m['logfc'] <- sign(fctop.m$value) * log2(abs(fctop.m$value)) / max(log2(abs(fctop.m$value)))


p <- ggplot(fctop.m, aes(variable, psite)) +
geom_tile(aes(fill = logfc)) +
scale_fill_gradient2(low = 'orangered3', mid = 'white', high = 'steelblue') +
xlab('Treatment') +
ylab('Phosphorylation sites') +
theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank()
    )

ggsave('gg_fcdiff_top50_heatmap.pdf', device = cairo_pdf, width = 6, height = 14)

### top heatmap with dendrogram

