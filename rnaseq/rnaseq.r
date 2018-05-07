#!/usr/bin/env Rscript

# Dénes Türei EMBL 2018
# turei.denes@gmail.com

require(readr)
require(dplyr)
require(ggplot2)

example <- suppressMessages(read_tsv(
        'GSM2648244_GAR1056.counts.tsv',
        col_names = c('name', 'count')
    )) %>%
    left_join(
        suppressMessages(read_tsv('mouse_biotypes.tsv')) %>%
        group_by(name, biotype) %>%
        summarize_all(first) %>%
        ungroup(),
        by = c('name')
    )

example %>%
    group_by(biotype) %>%
    mutate(cnt = n()) %>%
    summarize_all(first) %>%
    select(biotype, cnt) %>%
    write_tsv('biotypes.tsv')

p <- ggplot(example, aes(log2(count + 1))) +
    geom_density() +
    #scale_x_log10() +
    ylab('Frequency') +
    xlab('Count (log2)') +
    theme_minimal() +
    theme(text = element_text(family = 'DINPro'))

ggsave('counts_density.pdf', device = cairo_pdf, width = 5, height = 4)

