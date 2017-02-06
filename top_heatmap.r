#!/usr/bin/Rscript

#
# # Dénes Türei EMBL-EBI 2017
#

require(ggplot2)
require(reshape2)
require(ggdendro)
require(grid)
#require(cowplot)
#require(egg)

source('combined.plot.incl.r')

fctop <- get_fctop()
func <- get_func()

# with effect signs
combined_plot(fctop = fctop,
              func = func,
              top = 100,
              hheatmap = .9584,
              wheatmap = .38,
              xheatmap = .2,
              yheatmap = .48,
              hannotat = .9154,
              wannotat = .3,
              xannotat = .502,
              yannotat = .5194,
              hdendrog = .92,
              wdendrog = .2,
              xdendrog = .7435,
              ydendrog = .5295,
              hpaper   = 14,
              wpaper   = 8)

# without considering effect signs
combined_plot(fctop = fctop,
              func = func,
              top = 100,
              signs = FALSE,
              hheatmap = .9584,
              wheatmap = .38,
              xheatmap = .2,
              yheatmap = .48,
              hannotat = .9154,
              wannotat = .3,
              xannotat = .502,
              yannotat = .5194,
              hdendrog = .92,
              wdendrog = .2,
              xdendrog = .7435,
              ydendrog = .5295,
              hpaper   = 14,
              wpaper   = 8)

# all with effect signs
combined_plot(fctop = fctop,
                func = func,
                top = FALSE,
                hheatmap = .931,
                wheatmap = .38,
                xheatmap = .2,
                yheatmap = .5049,
                hannotat = .9154,
                wannotat = .3,
                xannotat = .502,
                yannotat = .5194,
                hdendrog = .975,
                wdendrog = .2,
                xdendrog = .7435,
                ydendrog = .523,
                hpaper   = 38,
                wpaper   = 8)

# all not using effect signs
combined_plot(fctop = fctop,
              func = func,
              top = FALSE,
              signs = FALSE,
              hheatmap = .93125,
              wheatmap = .38,
              xheatmap = .2,
              yheatmap = .5116,
              hannotat = .9227,
              wannotat = .3,
              xannotat = .502,
              yannotat = .5194,
              hdendrog = .9978,
              wdendrog = .2,
              xdendrog = .7435,
              ydendrog = .5215,
              hpaper   = 70,
              wpaper   = 8)


