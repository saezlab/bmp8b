#!/usr/bin/env Rscript

# Denes Turei EMBL 2018
# turei.denes@gmail.com

require(dplyr)
require(readr)
require(tidyr)
require(viper)


fc2fc <- function(fc){
    
    ifelse(
        fc >= 0.0,
        fc,
        abs(1 / fc)
    )
    
}


expfile  <- 'qpcr.tsv'
exp2file <- 'exp2.tsv'
upgsfile <- 'up_gs.tsv'

excl <- c('B2m', 'Betaactin', 'BK', 'Bmp5', 'Erbb3', 'F480', 'Th', 'Tnfa')

exp  <- suppressMessages(read_tsv(expfile))

upgs <- suppressMessages(read_tsv(upgsfile))

d <- suppressMessages(read_tsv(exp2file)) %>%
    left_join(
        suppressMessages(read_tsv(expfile)) %>%
            select(gs_mouse, uniprot = up_human) %>%
            group_by(gs_mouse, uniprot) %>%
            summarize_all(first) %>%
            ungroup(),
        by = c('gs_mouse')
    ) %>%
    left_join(
        upgs,
        by = c('uniprot')
    ) %>%
    filter(!is.na(genesymbol))

dv <- d %>%
    arrange(weeks, tiss, gen, genesymbol, excelcol, excelrow) %>%
    group_by(weeks, tiss, gen, genesymbol) %>%
    mutate(sid = 1:n()) %>%
    ungroup() %>%
    mutate(
        sample = sprintf('%s_%s_%s_%i', weeks, tiss, gen, sid)
    ) %>%
    group_by(sample) %>%
    mutate(
        zscore1 = (raw - mean(raw)) / sd(raw)
    ) %>%
    ungroup() %>%
    select(sample, genesymbol, zscore1) %>%
    spread(sample, zscore1)

rn <- dv$genesymbol

dv <- dv %>% select(-genesymbol)

expm <- (dv %>% as.data.frame())[,which(colSums(is.na(dv)) < dim(dv)[1] / 3)]
rownames(expm) <- rn
expm <- expm[which(rowSums(!is.na(expm)) > 0),]

tfact <- viper(
    eset = expm,
    regulon = viper_regulon,
    nes = TRUE,
    method = 'none',
    minsize = 4,
    eset.filter = FALSE
)

dm <- d %>%
    filter(tissue == 'TG BAT' & weeks == '5w') %>%
    select(genesymbol, pval_w, fc_med) %>%
    group_by(genesymbol) %>%
    summarize_all(first) %>%
    ungroup() %>%
    mutate(fc = fc2fc(fc_med)) %>%
    mutate(logfc = log2(fc)) %>%
    mutate(sig = qnorm(pval_w / 2, lower.tail = FALSE) * sign(logfc)) %>%
    arrange(desc(sig)) %>%
    filter(!is.na(sig))

sig <- dm$sig
names(sig) <- dm$genesymbol

mrs <- msviper(
    ges = sig,
    regulon = viper_regulon,
    minsize = 4,
    ges.filter = FALSE
)

tfact2 <- as_tibble(data.frame(
        regulon = names(mrs$es$nes),
        size = mrs$es$size[ names(mrs$es$nes) ], 
        nes = mrs$es$nes, 
        reg_pval = mrs$es$p.value, 
        reg_fdr = p.adjust(mrs$es$p.value, method = 'fdr')
    )) %>%
    mutate(
        genesymbol = gsub('_.*', '', regulon)
    ) %>%
    left_join(
        upgs,
        by = c('genesymbol')
    ) %>%
    write_tsv('tfact.tsv')

