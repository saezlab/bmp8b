#!/usr/bin/Rscript

#
# # Dénes Türei EMBL-EBI 2017
#

require(ggplot2)
require(reshape2)
require(ggdendro)
require(grid)

### setting up constants
shortfuncnames <- list(
    angiogenesis = 'A',
    inflammation = 'I',
    neurogenesis = 'N'
    )

shortfuncnames[['lipid metabolism']] <- 'L'
shortfuncnames[['survival, cell cycle']] <- 'C'

# setting up a blank theme
theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour = NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), 'cm')
    )


### top heatmap with dendrogram

get_fctop <- function(fname = 'fctop_none.csv'){
    fctop <- read.table(fname, header = TRUE, sep = '\t')
    
    fctop <- fctop[, names(fctop) != 'resnum']
    
    return(fctop)
}

get_func <- function(fname = 'functional.csv'){
    func <- read.table('functional.csv', sep = '\t', header = TRUE)
    
    newlevels <- NULL
    for(fn in levels(func$func)){
        newlevels <- c(newlevels, shortfuncnames[[fn]])
    }
    
    levels(func$func) <- newlevels
    
    return(func)
}

plot_dendrogram <- function(dendrodata){
    p <- ggplot(segment(dendrodata)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.1) +
        #scale_y_continuous(limits = c(0, 100)) +
        scale_x_continuous(
            #limits = c(0, 100)
        ) +
        coord_flip(
            #xlim = c(0, 100
        ) +
        theme_none
    
    return(p)
}

plot_functional <- function(data){
    p <- ggplot(data, aes(y = label, fill = value, x = category)) +
        facet_grid(. ~ func) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = 'white', high = '#333333') +
        xlab('Functional annotations') +
        ylab('') +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
            axis.text.y = element_blank(),
            legend.position = 'none',
            axis.ticks.y = element_blank(),
            panel.border = element_rect(color = '#777777', fill = NA, size = .1)
        )
    
    return(p)
}

plot_heatmap <- function(data){
    p <- ggplot(data, aes(x = variable, y = psite)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(high = '#990000', mid = 'white', low = '#333333',
                             limits = c(-2, 2),
                             guide = guide_colorbar(
                                 title = 'Fold change (log2)',
                                 direction = 'horizontal',
                                 label.position = 'top'
                             )
        ) +
        theme(legend.position = 'left') +
        xlab('Treatment') +
        ylab('Phosphorylation site') +
        #guides(fill = guide_legend(title = 'Fold change (log2)')) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
            panel.background = element_blank(),
            legend.position = 'bottom'
        )
    
    return(p)
}

get_pdfname <- function(wfunc, signs, top){
    return(
        paste(
            'gg_heatmap_dendro',
            ifelse(wfunc, '_func', ''),
            ifelse(signs, '_effect-signs', ''),
            ifelse(top, sprintf('_%d', top), '_all'),
            '.pdf',
            sep = ''
        )
    )
}

combined_plot <- function(fctop,
                          func,
                          top = 100,
                          signs = TRUE,
                          hheatmap = .9584,
                          wheatmap = .38,
                          xheatmap = .2,
                          yheatmap = .48,
                          hannotat = .9154,
                          wannotat = .3,
                          xannotat = .502,
                          yannotat = .5194,
                          hdendrog = .9297,
                          wdendrog = .2,
                          xdendrog = .7435,
                          ydendrog = .5252,
                          hpaper   = 14,
                          wpaper   = 8){
    
    if(signs){
        # here we drop all items where the effect sign is unknown
        signedfctop <- fctop[fctop$effect != 0,]
        # and set the signs of the fold changes based on the effect signs
        signedfctop[,9:11] <-  signedfctop[,9:11] * signedfctop$effect
        if(top & typeof(top) == 'double'){
            headfctop <- head(signedfctop, top)
        }else{
            headfctop <- signedfctop
        }
    }else{
        # here we use all rows and the fold changes with their original signs
        if(top & typeof(top) == 'double'){
            headfctop <- head(fctop, top)
        }else{
            headfctop <- fctop
        }
    }
    rnames <- headfctop$label
    rownames(headfctop) <- rnames
    mheadfctop <- as.matrix(headfctop[,9:11])
    rownames(mheadfctop) <- rnames
    
    #return(mheadfctop)
    psites_dendro <- as.dendrogram(hclust(dist(mheadfctop / apply(abs(mheadfctop), 1, max))))
    ordr <- order.dendrogram(psites_dendro)
    oheadfctop <- mheadfctop[ordr,]
    onames <- attr(oheadfctop, 'dimnames')
    dfheadfctop <- as.data.frame(oheadfctop)
    colnames(dfheadfctop) <- onames[[2]]
    dfheadfctop$psite <- onames[[1]]
    dfheadfctop$psite <- with(dfheadfctop, factor(psite, levels = psite, ordered = TRUE))
    mltdfheadfctop <- melt(dfheadfctop, id.vars = 'psite')
    dendrodatapsite <- dendro_data(psites_dendro)
    
    # the functional annotations
    # selecting only those for the top
    headfunc <- func[func$label %in% rownames(dfheadfctop),]
    # reaordering factor levels
    headfunc$label <- factor(headfunc$label, levels = onames[[1]])

    # plotting the heatmap
    p1 <- plot_heatmap(mltdfheadfctop)
    # plotting the dendrogram
    p3 <- plot_dendrogram(dendrodatapsite)
    # plotting functional annotations
    p4 <- plot_functional(headfunc)
    
    #ggsave('gg_func_table_dendro-order.pdf', device = cairo_pdf, width = 5, height = 12)
    
    # combined plot of only the dendrogram and fold changes heatmap
    fname <- get_pdfname(FALSE, signs, top)
    
    cat(sprintf('\t:: Plotting to %s\n', fname))
    cat(sprintf('\t:: Values at NFKB1-S927: %s\n', paste(as.character(headfctop[headfctop$psite == 'NFKB1_S927',]))))
    
    cairo_pdf(filename = fname, width = 8, height = 14)
        grid.newpage()
        print(p1, vp = viewport(0.8, 1.0, x = .4, y = .5))
        print(p3, vp = viewport(0.2, 1.076, x = 0.88, y = 0.498))
    dev.off()
    
    # combined plot of all 3 plots
    fname <- get_pdfname(TRUE, signs, top)
    cat(sprintf('\t:: Plotting to %s\n', fname))
    cairo_pdf(filename = get_pdfname(TRUE, signs, top),
              width = wpaper,
              height = hpaper)
        grid.newpage()
        print(p4, vp = viewport(wannotat, hannotat, x = xannotat, y = yannotat))
        print(p1, vp = viewport(wheatmap, hheatmap, x = xheatmap, y = yheatmap))
        print(p3, vp = viewport(wdendrog, hdendrog, x = xdendrog, y = ydendrog))
    dev.off()
    
#     p1g <- ggplotGrob(p1)
#     p4g <- ggplotGrob(p4)
#     p3g <- ggplotGrob(p3)
#     
#     plots <- list(p1, p4, p3)
#     heights <- list()
#     grobs <- list()
#     for(i in 1:length(plots)){
#         grobs[[i]] <- ggplot_gtable(ggplot_build(plots[[i]]))
#         heights[[i]] <- grobs[[i]]$height[1:3]
#     }
#     maxheight <- do.call(grid::unit.pmax, heights)
#     for (i in 1:length(grobs)){
#         grobs[[i]]$heights[1:3] <- as.list(maxheight)
#     }
    
    # solution with egg
#     cairo_pdf(filename = paste('gg_heatmap_dendro_func',
#                                ifelse(signs, '_effect-signs', ''),
#                                '_egg',
#                                '.pdf',
#                                sep = ''),
#               width = wpaper,
#               height = hpaper)
#         ggarrange(plots = list(p1, p4, p3), ncol = 3, widths = c(.38, .3, .2))
#     dev.off()
    
    # solution with cowplot
    #plot_grid(plotlist = grobs, rel_widths = c(.38, .3, .2), labels = 'AUTO', ncol = 3)
    # solution with gridExtra::grid.arrange
    # grid.arrange(grobs = grobs, ncol = 3, widths = c(.38, .3, .2))
#     ggsave(paste('gg_heatmap_dendro_func',
#                  ifelse(signs, '_effect-signs', ''),
#                  '_cowplot',
#                  '.pdf',
#                 sep = ''),
#            device = cairo_pdf,
#            width = wpaper,
#            height = hpaper)
    
    # end of block: using signs or not
}
