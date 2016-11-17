#!/usr/bin/Rscript

library(ggplot2)
library(RColorBrewer)
library(GGally)
library(reshape2)

cat("\t:: Reading data.\n")

getPalette = colorRampPalette(brewer.pal(12, "Paired"))

data <- read.csv('bmp8.csv', sep = '\t', header = TRUE)

mdata <- melt(data, measure.vars = c('norm', 'cv', 'signal'), variable.name = 'var')

pdata <- melt(data, measure.vars = c('pratio', 'pratio_actin', 'pratio_gapdh', 'pratio_pkc'), variable.name = 'pratio')
fdata <- melt(data, measure.vars = c('fc', 'fc_actin', 'fc_gapdh', 'fc_pkc'), variable.name = 'fc')
ndata <- melt(data, measure.vars = c('norm', 'norm_actin', 'norm_gapdh', 'norm_pkc'), variable.name = 'norm')

p <- ggplot(data, aes(x = phos, y = log10(signal))) +
geom_violin() +
facet_wrap(~group) +
scale_x_discrete(labels = c('p' = 'P', 'np' = 'non-P')) +
xlab('Phosphorylated [no/yes]') +
ylab('Signal [log10]') +
ggtitle('Raw signal of spots by phosphorylation state and treatment') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_phos_signal_viol.pdf', device = cairo_pdf, width = 12, height = 12)

###

p <- ggpairs(data[,c('numof_kin', 'degree', 'group', 'signal', 'cv', 'norm', 'phos', 'pratio', 'fc')], mapping = aes(colour = group, alpha = 0.4))

ggsave('gg_ggpairs.pdf', device = cairo_pdf, width = 18, height = 18)

###

p <- ggplot(data, aes(x = group, y = fc)) +
geom_boxplot() +
xlab('Treatment') +
ylab('Fold change') +
ggtitle('Fold change vs. treatment') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_fc_box.pdf', device = cairo_pdf, width = 12, height = 12)

###

p <- ggplot(fdata, aes(x = group, y = value)) +
geom_boxplot() +
facet_wrap(~fc) +
# scale_x_discrete(labels = c('p' = 'P', 'np' = 'non-P')) +
xlab('Treatment') +
ylab('Fold change') +
ggtitle('Fold change vs. treatment\nvs. normalization by standards') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_fc_std_box.pdf', device = cairo_pdf, width = 6, height = 12)

###

p <- ggplot(pdata, aes(x = group, y = value)) +
geom_boxplot() +
facet_wrap(~pratio) +
# scale_x_discrete(labels = c('p' = 'P', 'np' = 'non-P')) +
xlab('Treatment') +
ylab('Phosphorylation [relative]') +
ggtitle('Phosphorylation state vs. treatment\nvs. normalization by standards') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_phos_std_box.pdf', device = cairo_pdf, width = 6, height = 12)

###

p <- ggplot(pdata, aes(x = group, y = log10(value))) +
geom_boxplot() +
facet_wrap(~pratio) +
# scale_x_discrete(labels = c('p' = 'P', 'np' = 'non-P')) +
xlab('Treatment') +
ylab('Phosphorylation [log]') +
ggtitle('Phosphorylation state vs. treatment\nvs. normalization by standards') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_log-phos_std_box.pdf', device = cairo_pdf, width = 6, height = 12)

###

p <- ggplot(ndata, aes(x = group, y = log(value))) +
geom_boxplot() +
facet_wrap(~norm) +
# scale_x_discrete(labels = c('p' = 'P', 'np' = 'non-P')) +
xlab('Treatment') +
ylab('Normalized value [log]') +
ggtitle('Value vs. treatment\nvs. normalization by standards') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

dcast(data, uniprot + resaa + resnum + group ~ )

ggsave('gg_group_log-value_std_box.pdf', device = cairo_pdf, width = 6, height = 12)


###
fc <- dcast(data, uniprot + name + gsymbol + resaa + resnum ~ group, value.var = 'fc_gapdh', fun.aggregate = max)

ggplot(fc, aes(x = BMP8b_NE, y = BMP8b)) +
geom_point(data = fc, aes(x = BMP8b_NE, y = NE), color = 'gray') +
geom_point() +
xlab('BMP8b & NE') +
ylab('BMP8b or NE') +
ggtitle('Fold change') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_group_fc_point.pdf', device = cairo_pdf, width = 6, height = 6)

###

fc <- dcast(data, uniprot + name + gsymbol + resaa + resnum ~ group, value.var = 'fc_pkc', fun.aggregate = max)

fc.pc <- prcomp(fc[c('BMP8b', 'NE', 'BMP8b_NE')], center = TRUE, scale. = TRUE)
ggbiplot(fc.pc, obs.scale = 1, var.scale = 1, circle = TRUE) +
ggtitle('PCA of fold changes at treatments')
ggsave('gg_fc_pca_biplot_c.pdf', device = cairo_pdf, width = 6, height = 6)

###

pr <- dcast(data, uniprot + name + gsymbol + resaa + resnum ~ group, value.var = 'pratio_actin', fun.aggregate = max)

pr.pc <- prcomp(pr[c('BMP8b', 'NE', 'BMP8b_NE')], center = TRUE, scale. = TRUE)
ggbiplot(pr.pc, obs.scale = 1, var.scale = 1, circle = TRUE) +
ggtitle('PCA of ratios of phosphorylation at treatments')
ggsave('gg_phos_pca_biplot_c.pdf', device = cairo_pdf, width = 6, height = 6)

###

mcor <- t(apply(cor(pr[,6:9]), 2, rev))
mmcor <- melt(mcor)
ggplot(mmcor, aes(Var1, Var2)) +
geom_tile(aes(fill = value), colour = "white") +
geom_text(aes(label = round(value, 3)), size = 7) +
scale_fill_gradient(low = "white", high = "steelblue", guide = guide_legend(title = 'Corr. coeff.')) +
xlab('Treatments') +
ylab('Treatments') +
ggtitle('Correlations between phosphorylation patterns at treatments') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 21))

ggsave('gg_phos_corr_tile.pdf', device = cairo_pdf, width = 12, height = 10)

###

mcor <- t(apply(cor(fc[,6:9]), 2, rev))
mmcor <- melt(mcor)
ggplot(mmcor, aes(Var1, Var2)) +
geom_tile(aes(fill = value), colour = "white") +
geom_text(aes(label = round(value, 3)), size = 7) +
scale_fill_gradient(low = "white", high = "steelblue", guide = guide_legend(title = 'Corr. coeff.')) +
xlab('Treatments') +
ylab('Treatments') +
ggtitle('Correlations between fold changes at treatments') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 21))

ggsave('gg_fc_corr_tile.pdf', device = cairo_pdf, width = 12, height = 10)

###

tfc <- t(fc[,6:9])
colnames(tfc) <- paste(fc$gsymbol, '-', fc$resaa, fc$resnum, sep = '')
mcor <- t(apply(cor(tfc), 2, rev))
mmcor <- melt(mcor)

###

hc <- hclust(dist(t(tfc)), 'ward.D2')
hcdata <- dendro_data(hc, type = 'rectangle')

ggplot() +
geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
geom_text(data=label(hcdata), aes(x=x, y=y, label=label, hjust=0), size=2) +
coord_flip() +
scale_y_reverse(expand=c(0.2, 0)) +
xlab('Phosphorylation sites') +
ylab('Fold changes (distance)') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.y = element_blank())

ggsave('gg_fc_dendro.pdf', device = cairo_pdf, width = 10, height = 49)

###

tpr <- t(pr[,6:9])
colnames(tpr) <- paste(pr$gsymbol, '-', pr$resaa, pr$resnum, sep = '')
mcor <- t(apply(cor(tpr), 2, rev))
mmcor <- melt(mcor)

###

hc <- hclust(dist(t(tpr)), 'ward.D2')
hcdata <- dendro_data(hc, type = 'rectangle')

ggplot() +
geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
geom_text(data=label(hcdata), aes(x=x, y=y, label=label, hjust=0), size=2) +
coord_flip() + 
scale_y_reverse(expand=c(0.2, 0)) +
xlab('Phosphorylation sites') +
ylab('Phosphorylation patterns (distance)') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.y = element_blank())

ggsave('gg_psite_dendro.pdf', device = cairo_pdf, width = 10, height = 49)
