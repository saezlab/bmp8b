#!/usr/bin/env Rscript

# Denes Turei EMBL 2018
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(readr)
require(purrr)
require(furrr)
require(tidyr)
require(rlang)
require(scales)

rwfile    <- 'rw_phos_tf.tsv'
rwfile    <- 'rw_signalink.tsv'
rwfile    <- 'rw_ca1.tsv'
rwfile    <- 'rw_s3.tsv'
rwfile    <- 'rw_s2.tsv'
rwfile    <- 'rw_spike.tsv'
rwfile    <- 'rw_signor.tsv'
spfile    <- 'sp_phos_tf.tsv'
phosfile  <- 'phos.tsv'
tffile    <- 'tfreg.tsv'
expfile   <- 'exp.tsv'
tfactfile <- 'tfact.tsv'
upgsfile  <- 'up_gs.tsv'


dexp  <- suppressMessages(read_tsv(expfile))
dtf   <- suppressMessages(read_tsv(tffile))
dphos <- suppressMessages(read_tsv(phosfile))
dsp   <- suppressMessages(read_tsv(spfile))
tfact <- suppressMessages(read_tsv(tfactfile))
upgs  <- suppressMessages(read_tsv(upgsfile))


rw.read_data <- function(rwfile, w = '5w'){
    
    excl <- c('B2m', 'Betaactin', 'BK', 'Bmp5', 'Erbb3', 'F480', 'Th', 'Tnfa')
    
    # joining expression and tf-target tables
    d <- dexp %>%
        filter(!(gs_mouse %in% excl)) %>%
        select(
            gs_mouse, weeks, tissue, tiss, gen,
            pval_t, pval_w, fc_mean, fc_med,
            target = up_human
        ) %>%
        group_by(weeks, target, tissue) %>%
        mutate(
            min_pval_w = min(pval_w)
        ) %>%
        ungroup() %>%
        filter(!pval_w > min_pval_w) %>%
        select(-min_pval_w) %>%
        left_join(
            dtf %>%
                select(tf = uptf, target = uptg, sign, level, curated),
            by = c('target')
        ) %>%
        filter(!is.na(level)) %>%
        group_by(target) %>%
        mutate(
            tfpertarget = n()
        ) %>%
        ungroup() %>%
        left_join(
            tfact %>%
                select(tf = uniprot, nes, reg_pval, reg_fdr),
            by = c('tf')
        ) %>%
        filter(weeks == w)
    
    # adding random walk data
    d <- d %>%
        left_join(
            drw   <- suppressMessages(read_tsv(rwfile)),
            by = c('tf')
        ) %>%
        filter(dist > 0)
    
    # adding shortest path data
    d <- d %>%
        left_join(dsp, by = c('tf', 'phos'))
    
    # adding phosphoassay data
    d <- d %>%
        left_join(
            dphos %>% select(phos = uniprot, gs_phos = genesymbol, effect, fc),
            by = c('phos')
        ) %>%
        filter(!is.na(fc))
    
    return(d)
    
}


rw.tfpertarget <- function(d){
    
    # plotting tf per target
    p <- ggplot(
            d %>% group_by(target) %>% summarize_all(first),
            aes(tfpertarget)
        ) +
        geom_density() +
        xlab('Number of TFs for one target') +
        ylab('Density') +
        theme_linedraw() +
        theme(
            text = element_text(family = 'DINPro')
        )
    
    ggsave('tf-per-target.pdf', device = cairo_pdf, width = 4, height = 3)
    
}

rw.tfactivity_plot <- function(d){
    
    p <- ggplot(
            d %>%
                filter(!is.na(nes)) %>%
                group_by(tf, target) %>%
                summarize_all(first),
            aes(x = nes, y = fc_med)
        ) +
        geom_point() +
        xlab('TF normalized enrichment score') +
        ylab('Target exp. fold change') +
        theme_linedraw() +
        theme(text = element_text(family = 'DINPro'))
    
    ggsave('tfact_expfc.pdf', device = cairo_pdf, width = 4, height = 4)
    
}


rw.altered <- function(d){
    
    (
        # grouping by altered proteins and tfs
        d %>%
        mutate(
            phosalt = abs(fc) >= 1.5,
            tfalt   = pval_t <= .05,
            phostf  = ifelse(
                phosalt,
                ifelse(
                    tfalt,
                    'P(+) --> TF(+)',
                    'P(+) --> TF(-)'
                ),
                ifelse(
                    tfalt,
                    'P(-) --> TF(+)',
                    'P(-) --> TF(-)'
                )
            ),
            alt = phosalt & tfalt
        ) %>%
        {`if`(
            !('dist_o' %in% names(d)),
            mutate(., dist_o = dist),
            .
        )} %>%
        group_by(phosalt) %>%
        mutate(dist = dist_o / mean(dist_o)) %>%
        ungroup() %>%
        mutate(
            tfscore = -log(pval_w) * dist,
            phscore = abs(fc) * dist
        )
    )
    
}


rw.nes_phos <- function(d){
    
    d0 <- d %>%
        filter(!is.na(nes)) %>%
        group_by(tf, phos) %>%
        summarize_all(first) %>%
        ungroup() %>%
        arrange(nes) %>%
        mutate(nes = factor(nes, levels = unique(nes), ordered = TRUE))
    
    p <- ggplot(d0, aes(x = nes, y = log(dist), color = phosalt)) +
        geom_boxplot(outlier.size = .5, outlier.shape = 18, lwd = .3) +
        scale_color_manual(
            guide = guide_legend(
                title = 'Signaling proteins\nphosphorylation state'),
            values = c(
                `TRUE`  = '#CC0000',
                `FALSE` = '#000000'
            ),
            labels = c(
                `TRUE`  = 'Altered',
                `FALSE` = 'Non-altered'
            )
        ) +
        xlab('TF normalized enrichment score') +
        ylab('Proximity (log)') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINpro'),
            panel.grid.major = element_line(color = '#CCCCCC'),
            panel.grid.minor = element_line(color = '#CCCCCC'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        )
    
    ggsave('nes-phos-dist.pdf', device = cairo_pdf, width = 10, height = 5)
    
}

rw.param_discovery <- function(){
    
    rwfiles <- c(
        'rw_phos_tf.tsv',
        'rw_signalink.tsv',
        'rw_s3.tsv',
        'rw_s2.tsv',
        'rw_spike.tsv',
        'rw_signor.tsv'
    )
    
    for(rwfile in rwfiles){
        
        rw.param_discovery0(rwfile, reg_pval)
        rw.param_discovery0(rwfile, nes)
        rw.param_discovery0(rwfile, pval_t)
        rw.param_discovery0(rwfile, pval_w)
        rw.param_discovery0(rwfile, fc_mean)
        rw.param_discovery0(rwfile, fc_med)
        rw.param_discovery0(rwfile, fc)
        
    }
    
}

rw.param_discovery0 <- function(rwfile, param){
    
    param <- enquo(param)
    tparam <- quo_text(param)
    
    pdfname <- sprintf('%s_%s_param.pdf', rwfile, tparam)
    
    d <- rw.read_data(rwfile) %>%
        rw.altered() %>%
        filter(!is.na(!!tparam))
    
    qparam <- quantile(d[[tparam]], seq(.1, 1, .1), na.rm = TRUE)
    
    d <- d %>%
        mutate(
            bin = unlist(purrr::map(
                .[[tparam]],
                function(val){
                    i <- min(which(qparam >= val))
                    return(sprintf('#%g (%g)', i, `if`(i <= length(qparam), qparam[i], NA)))
                }
            ))
        )
    
    p <- ggplot(d, aes(x = bin, y = log(dist), color = phosalt)) +
        geom_boxplot(outlier.size = .5, outlier.shape = 18, lwd = .3) +
        scale_color_manual(
            guide = guide_legend(
                title = 'Signaling proteins\nphosphorylation state'),
            values = c(
                `TRUE`  = '#CC0000',
                `FALSE` = '#000000'
            ),
            labels = c(
                `TRUE`  = 'Altered',
                `FALSE` = 'Non-altered'
            )
        ) +
        xlab('Parameter bins') +
        ylab('Proximity (log)') +
        ggtitle(sprintf('%s, %s', rwfile, tparam)) +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINpro'),
            panel.grid.major = element_line(color = '#CCCCCC'),
            panel.grid.minor = element_line(color = '#CCCCCC'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        )
    
    ggsave(pdfname, device = cairo_pdf, width = 21, height = 5)
    
    print(rwfile)
    print(tparam)
    print(dim(d))
    
}


p <- ggplot(d, aes(x = phostf, y = dist)) +
    geom_boxplot(outlier.size = .5, outlier.shape = 18, lwd = .3) +
    facet_grid(. ~ weeks) +
    theme_linedraw() +
    ylab('Proximity by random walks with return') +
    xlab('Altered vs. non-altered TFs and signaling proteins') +
    theme(
        text = element_text(family = 'DINPro'),
        panel.grid.major = element_line(color = '#CCCCCC'),
        panel.grid.minor = element_line(color = '#CCCCCC'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave('dist-boxplot-phos-tf.pdf', device = cairo_pdf, width = 7, height = 4)

d <- d %>%
    mutate(
        phostfscore = -log(pval_w) * dist * abs(fc)
    )

wilcox.test(
    (d %>% filter( phosalt))$dist %>% c,
    (d %>% filter(!phosalt))$dist %>% c
)
# W = 1.5416e+11, p-value < 2.2e-16

medphosalt <- median((d %>% filter( phosalt))$dist %>% c)
medphosnal <- median((d %>% filter(!phosalt))$dist %>% c)
(medphosalt - medphosnal) / medphosnal * 100 # 18.1 %

wilcox.test(
    (d %>% filter( tfalt))$dist %>% c,
    (d %>% filter(!tfalt))$dist %>% c
)
# W = 3.9882e+11, p-value < 2.2e-16

medtfalt <- median((d %>% filter( tfalt))$dist %>% c)
medtfnal <- median((d %>% filter(!tfalt))$dist %>% c)
(medtfalt - medtfnal) / medtfnal * 100 # -8.709282 %

# are the altered phosphosites closer to the altered TFs
# compared to the non altered phosphosites?
wilcox.test(
    (d %>% filter( tfalt &  phosalt))$splen %>% c,
    (d %>% filter( tfalt & !phosalt))$splen %>% c
)
# W = 7725900000, p-value < 2.2e-16
mp <- median((d %>% filter( tfalt &  phosalt))$dist %>% c)
mn <- median((d %>% filter( tfalt & !phosalt))$dist %>% c)
(mp - mn) / mn * 100 # 19.11203
# yes, their median proximity is 19.1 % higher

# are the altered TFs closer to the altered phosphosites
# compared to the non altered TFs?
wilcox.test(
    (d %>% filter( !tfalt & phosalt))$dist %>% c,
    (d %>% filter(  tfalt & phosalt))$dist %>% c
)
# W = 1959600000, p-value = 1.631e-05
mp <- median((d %>% filter(  tfalt & phosalt))$dist %>% c)
mn <- median((d %>% filter( !tfalt & phosalt))$dist %>% c)
(mp - mn) / mn * 100 # -8.222247 % lower
# no, their median proximity is -8.2 % lower

# are the altered phosphosites closer to the non altered TFs
# compared to the non altered phosphosites?
wilcox.test(
    (d %>% filter( !tfalt &  phosalt))$dist %>% c,
    (d %>% filter( !tfalt & !phosalt))$dist %>% c
)
# W = 9.286e+10, p-value < 2.2e-16
mp <- median((d %>% filter( !tfalt &  phosalt))$dist %>% c)
mn <- median((d %>% filter( !tfalt & !phosalt))$dist %>% c)
(mp - mn) / mn * 100 # 18.39658
# yes, their median proximity is 18.4 % higher

# are the altered TFs closer to the non altered phosphosites
# compared to the non altered TFs?
wilcox.test(
    (d %>% filter( !tfalt & !phosalt))$dist %>% c,
    (d %>% filter(  tfalt & !phosalt))$dist %>% c
)
# W = 3.5772e+11, p-value < 2.2e-16
mp <- median((d %>% filter(  tfalt & !phosalt))$dist %>% c)
mn <- median((d %>% filter( !tfalt & !phosalt))$dist %>% c)
(mp - mn) / mn * 100 # -8.773514 % lower
# no, their median proximity is -8.8 % lower

pos <- (d %>% filter(  tfalt & phosalt))$dist %>% c
neg <- (d %>% filter(  !tfalt | !phosalt))$dist %>% c
mp <- mean(pos)
mn <- mean(neg)
(mp - mn) / mn * 100

mp <- mean((d %>% filter(  tfalt & phosalt))$splen %>% c)
mn <- mean((d %>% filter( !tfalt | !phosalt))$splen %>% c)
(mp - mn) / mn * 100

ranks <- data.frame(
    x = seq(.1, 1, .01),
    alt  = rev(log(quantile(
        (d %>% filter(  tfalt &  phosalt))$dist %>% c,
        seq(.1, 1, .01)
    ))),
    nalt = rev(log(quantile(
        (d %>% filter( !tfalt | !phosalt))$dist %>% c,
        seq(.1, 1, .01)
    )))
) %>%
gather(alt, logprox, -x)

ranks <- data.frame(
    x = seq(.1, 1, .01),
    alt  = rev(quantile(
        (d %>% filter(  tfalt &  phosalt))$splen %>% c,
        seq(.1, 1, .01)
    )),
    nalt = rev(quantile(
        (d %>% filter(  !tfalt | !phosalt))$splen %>% c,
        seq(.1, 1, .01)
    ))
) %>%
gather(alt, logprox, -x)

p <- ggplot(ranks, aes(x = x, y = logprox, color = alt)) +
    geom_line() +
    scale_color_manual(
        guide = guide_legend(
            title = 'Proximity of signaling\nproteins and TFs (log)'),
        values = c(
            `alt`  = '#CC0000',
            `nalt` = '#000000'
        ),
        labels = c(
            `alt`  = 'Altered phosphosites\nto altered TFs',
            `nalt` = 'All others'
        )
    ) +
    xlab('Ranks') +
    ylab('Proximity (log)') +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINpro'),
        panel.grid.major = element_line(color = '#CCCCCC'),
        panel.grid.minor = element_line(color = '#CCCCCC')
    )

ggsave('prox-omnipath5-fc1.5p-p05.pdf', device = cairo_pdf, width = 4.2, height = 2.7)

p <- ggplot(d, aes(y = log(dist_o), x = gs_mouse, color = phosalt)) +
    geom_boxplot(outlier.size = .5, outlier.shape = 18, lwd = .3) +
    scale_color_manual(
        guide = guide_legend(
            title = 'Proximity to altered\nsignaling proteins'),
        values = c(
            `TRUE`  = '#CC0000',
            `FALSE` = '#000000'
        ),
        labels = c(
            `TRUE`  = 'Altered',
            `FALSE` = 'Non-altered'
        )
    ) +
    xlab('Downstream expressed proteins') +
    ylab('Proximity (log)') +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINpro'),
        panel.grid.major = element_line(color = '#CCCCCC'),
        panel.grid.minor = element_line(color = '#CCCCCC'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave('prox-proteins-omnipath5-fc1.5p-p05.pdf', device = cairo_pdf, width = 21, height = 2.7)


ggplot(d, aes(x = splen, color = tfalt & phosalt)) +
    geom_density()
    

# NRG4 TF activity

nrg4tfact <- tfact %>%
    rename(gstf = genesymbol) %>%
    left_join(dtf, by = c('gstf')) %>%
    filter(gstg == 'NRG4')

# NRG4 phospho
nrg4phos <- d %>%
    filter(gs_mouse == 'Nrg4') %>%
    arrange(desc(abs(fc))) %>%
    group_by(phos, tf) %>%
    summarize_all(first) %>%
    ungroup() %>%
    arrange(desc(dist_o)) %>%
    group_by(phos) %>%
    summarize_all(first) %>%
    ungroup() %>%
    arrange(desc(dist_o)) %>%
    select(gs_mouse, nes, phos, dist_o, fc, gs_phos, tf) %>%
    left_join(upgs %>% rename(tf = uniprot), by = c('tf')) %>%
    rename(gs_tf = genesymbol) %>%
    head(100) %>%
    mutate(gs_phos = factor(gs_phos, levels = unique(gs_phos), ordered = TRUE))

p <- ggplot(nrg4phos, aes(y = fc, x = gs_phos)) +
    geom_col(fill = '#000000') +
    geom_line(mapping = aes(x = seq(1, 100), y = log(dist_o) / 18 + 2.5), color = '#CC0000') +
    xlab('Signaling proteins') +
    ylab('Phosphorylation fold change') +
    scale_y_continuous(limits = c(1, NA), oob = rescale_none) +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)
    )

ggsave('nrg4-phos-top100.pdf', device = cairo_pdf, width = 15, height = 4)

nrg4phos30 <- nrg4phos %>% head(30)


p <- ggplot(nrg4phos30, aes(y = fc, x = gs_phos)) +
    geom_col(fill = '#000000') +
    geom_line(mapping = aes(x = seq(1, 30), y = log(dist_o) / 18 + 2.5), color = '#CC0000') +
    xlab('Signaling proteins') +
    ylab('Phosphorylation fold change') +
    scale_y_continuous(limits = c(1, NA), oob = rescale_none) +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave('nrg4-phos-top30.pdf', device = cairo_pdf, width = 8, height = 4)

# VEGFA TF activity

vegfatfact <- tfact %>%
    rename(gstf = genesymbol) %>%
    left_join(dtf, by = c('gstf')) %>%
    filter(gstg == 'VEGFA') %>%
    filter(reg_pval < .2) %>%
    arrange(reg_pval) %>%
    mutate(gstf = factor(gstf, levels = unique(gstf), ordered = TRUE))

p <- ggplot(vegfatfact, aes(y = nes, x = gstf)) +
    geom_bar(fill = '#000000', stat = 'identity') +
    xlab('Transcription factors (TFs)') +
    ylab('TF activity\n(normalized enrichment score)') +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave('vegfa-tf-top.pdf', device = cairo_pdf, width = 3, height = 4)

# VEGFA phospho
vegfaphos <- d %>%
    filter(gs_mouse == 'Vegfa166') %>%
    arrange(desc(abs(fc))) %>%
    group_by(phos, tf) %>%
    summarize_all(first) %>%
    ungroup() %>%
    arrange(desc(dist_o)) %>%
    group_by(phos) %>%
    summarize_all(first) %>%
    ungroup() %>%
    arrange(desc(dist_o)) %>%
    select(gs_mouse, nes, phos, dist_o, fc, gs_phos, tf) %>%
    left_join(upgs %>% rename(tf = uniprot), by = c('tf')) %>%
    rename(gs_tf = genesymbol) %>%
    head(100) %>%
    mutate(gs_phos = factor(gs_phos, levels = unique(gs_phos), ordered = TRUE))

p <- ggplot(vegfaphos, aes(y = fc, x = gs_phos)) +
    geom_col(fill = '#000000') +
    geom_line(mapping = aes(x = seq(1, 100), y = log(dist_o) / 18 + 2.5), color = '#CC0000') +
    xlab('Signaling proteins') +
    ylab('Phosphorylation fold change') +
    scale_y_continuous(limits = c(1, NA), oob = rescale_none) +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)
    )

ggsave('vegfa-phos-top100.pdf', device = cairo_pdf, width = 15, height = 4)

vegfaphos30 <- vegfaphos %>% head(30)


p <- ggplot(vegfaphos30, aes(y = fc, x = gs_phos)) +
    geom_col(fill = '#000000') +
    geom_line(mapping = aes(x = seq(1, 30), y = log(dist_o) / 18 + 2.5), color = '#CC0000') +
    xlab('Signaling proteins') +
    ylab('Phosphorylation fold change') +
    scale_y_continuous(limits = c(1, NA), oob = rescale_none) +
    theme_minimal() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave('vegfa-phos-top30.pdf', device = cairo_pdf, width = 8, height = 4)
