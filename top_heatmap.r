#!/usr/bin/Rscript

#
# # Dénes Türei EMBL-EBI 2017
#

require(ggplot2)

inFile <- 'fctop_none.csv'
fctop <- read.table(inFile, header = TRUE, sep = '\t')

fctop <- fctop[, names(fctop) != 'resnum']

fctop.m <- melt(fctop)

fctop.m$psite <- factor(fctop$psite, levels = unique(as.character(fctop$psite)))

fctop.m['logfc'] <- sign(fctop.m$value) * log2(abs(fctop.m$value)) / max(log2(abs(fctop.m$value)))

p <- ggplot(fctop.m, aes(variable, psite)) +
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
