#!/usr/bin/Rscript

#
# # (c) Dénes Türei EMBL-EBI 2016-2017
# # turei.denes@gmail.com
# # thanks to Aurelien Dugourd
#

treatments <- c('NE', 'BMP8b', 'BMP8b_NE')
setspref <- c('c2.cp', 'c2.cp.kegg', 'c5.bp', 'h.all', 'c2.cp.reactome', 'c2.cp.biocarta')
setspref_all <- c('c2.cp', 'c5.bp', 'h.all')
setspref_all_nogo <- c('c2.cp', 'h.all')

### loading genesets
source('piano.incl.r')
geneSet_go   <- load_genesets(setspref_all)
geneSet_nogo <- load_genesets(setspref_all_nogo)

### running piano
s_nogo_result <- run_piano_series(treatments, geneSet_nogo, TRUE)
s_go_result   <- run_piano_series(treatments, geneSet_go,   TRUE)
u_nogo_result <- run_piano_series(treatments, geneSet_nogo, FALSE)
u_go_result   <- run_piano_series(treatments, geneSet_go,   FALSE)

### saving results
save.image(file = 'piano_all_together.Rdata')

### plotting heatmaps
source('piano.incl.r')
piano_plot_series(u_nogo_result, plotname = 'nosign_noGO_undir', dir_col = 3,
                  maintitle = '\t\tGeneset ranks\n\t\t(non-directed)',
                  consensus_args = list(cutoff = 50))

piano_plot_series(u_nogo_result, plotname = 'nosign_noGO_down2', dir_col = 2,
                  maintitle = '\t\t\tGeneset ranks\n\t\t\t(less phosphorylated)',
                  consensus_args = list(cutoff = 50))
piano_plot_series(u_nogo_result, plotname = 'nosign_noGO_up2', dir_col = 4,
                  maintitle = '\t\t\tGeneset ranks\n\t\t\t(more phosphorylated)',
                  consensus_args = list(cutoff = 50))

piano_plot_series(u_nogo_result, plotname = 'nosign_noGO_down', dir_col = 1,
                  maintitle = '\t\t\tGeneset ranks\n\t\t\t(less phosphorylated)',
                  consensus_args = list(cutoff = 50))
piano_plot_series(u_nogo_result, plotname = 'nosign_noGO_up', dir_col = 5,
                  maintitle = '\t\t\tGeneset ranks\n\t\t\t(more phosphorylated)',
                  consensus_args = list(cutoff = 50))


piano_plot_series(s_nogo_result, plotname = 'sign_noGO_down', dir_col = 1,
                  maintitle = '\t\tGeneset ranks\n\t\t(downregulated)',
                  consensus_args = list(cutoff = 50))
piano_plot_series(s_nogo_result, plotname = 'sign_noGO_up', dir_col = 5,
                  maintitle = '\t\tGeneset ranks\n\t\t(upregulated)',
                  consensus_args = list(cutoff = 50))

piano_plot_series(s_nogo_result, plotname = 'sign_noGO_undir', dir_col = 3,
                  maintitle = '\t\tGeneset ranks\n\t\t(non-directed)',
                  consensus_args = list(cutoff = 50))

piano_plot_series(s_nogo_result, plotname = 'sign_noGO_down2', dir_col = 2,
                  maintitle = '\t\tGeneset ranks\n\t\t(downregulated)',
                  consensus_args = list(cutoff = 50),
                  cairo_args = list(height = 20),
                  heatmap_args = list(lhei = c(0.05, 1.0)))
