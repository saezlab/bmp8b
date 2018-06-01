#!/usr/bin/env Rscript

# Denes Turei EMBL 2018
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(readr)
require(purrr)
require(furrr)
require(tidyr)

fc <- function(val, ctrl){
    
    ifelse(
        val == ctrl,
        0.0,
        `if`(
            ctrl < val,
            val / ctrl,
            -ctrl / val
        )
    )
    
}

fc2fc <- function(fc){
    
    ifelse(
        fc >= 0.0,
        fc,
        abs(1 / fc)
    )
    
}

fc_ci <- function(val, ctrl, aggr = mean){
    
    if(length(val) < 2 | length(ctrl) < 2){
        return(c(NA, NA))
    }
    
    if(length(val) > 9 | length(ctrl) > 9){
        return(fc_ci2(val, ctrl, aggr = aggr))
    }
    
    mctrl <- aggr(ctrl)
    
    pool <- c(apply(combn(val, m = ceiling(length(val) / 2)), 1, function(v1){
        av1 <- aggr(v1)
        apply(combn(ctrl, m = ceiling(length(ctrl) / 2)), 1, function(v2){
            av1 - aggr(v2)
        })
    }))
    
    vallo <- mctrl - sd(pool) * 1.645
    valhi <- mctrl + sd(pool) * 1.645
    
    return(c(fc(vallo, mctrl), fc(valhi, mctrl)))
    
}

fc_ci2 <- function(val, ctrl, aggr = mean){
    
    if(length(val) < 2 | length(ctrl) < 2){
        return(c(NA, NA))
    }
    
    mctrl <- aggr(ctrl)
    
    pool <- sapply(seq(100), function(i){
        aggr(base::sample(val,  ceiling(length(val)  / 2))) -
        aggr(base::sample(ctrl, ceiling(length(ctrl) / 2)))
    })
    
    vallo <- mctrl - sd(pool) * 1.645
    valhi <- mctrl + sd(pool) * 1.645
    
    return(c(fc(vallo, mctrl), fc(valhi, mctrl)))
    
}

infile <- 'qpcr.tsv'

qpcr.read_data <- function(){
    
    (
        suppressMessages(read_tsv(infile)) %>%
        mutate(tissue = recode(tissue, `WT SCW` = 'WT ScW')) %>%
        group_by(gs_mouse, weeks) %>%
        mutate(
            expmean = mean(raw),
            expsd = sd(raw),
            zscore = (raw - expmean) / expsd,
            other = gsub('TG|KO', 'WT', tissue)
        ) %>%
        ungroup() %>%
        group_by(gs_mouse, weeks, tissue, excelcol) %>%
        mutate(
            zmean = mean(zscore)
        ) %>%
        ungroup()
    )
    
}

d <- qpcr.read_data()

d %>% write_tsv('exp2.tsv')

conds <- d %>%
    filter(
        !startsWith(tissue, 'WT')
    ) %>%
    group_by(gs_mouse, tissue, weeks, other) %>%
    nest() %>%
    left_join(
        d %>%
            filter(startsWith(tissue, 'WT')) %>%
            group_by(gs_mouse, tissue, weeks, other) %>%
            nest(),
        by = c('gs_mouse', 'other', 'weeks'),
        suffix = c('_tr', '_co')
    ) %>%
    filter(!is.na(tissue_co)) %>%
    select(-tissue_co) %>%
    rename(tissue = tissue_tr) %>%
    mutate(
        pval_t = unlist(
            future_map2_dbl(
                data_co,
                data_tr,
                function(d1, d2){
                    `if`(
                        dim(d1)[1] > 1 & dim(d2)[1] > 1,
                        t.test(d1$raw, d2$raw)$p.value,
                        NA
                    )
                },
                .progress = TRUE
            )
        ),
        pval_w = unlist(
            future_map2_dbl(
                data_co,
                data_tr,
                function(d1, d2){
                    `if`(
                        dim(d1)[1] > 1 & dim(d2)[1] > 1,
                        wilcox.test(d1$raw, d2$raw)$p.value,
                        NA
                    )
                },
                .progress = TRUE
            )
        ),
        fc_med = unlist(
            future_map2_dbl(
                data_co,
                data_tr,
                function(d1, d2){
                    fc(median(d2$raw), median(d1$raw))
                }
            )
        ),
        fc_mean = unlist(
            future_map2_dbl(
                data_co,
                data_tr,
                function(d1, d2){
                    fc(mean(d2$raw), mean(d1$raw))
                }
            )
        ),
        fc_mean_ci = future_map2(
            data_co,
            data_tr,
            function(d1, d2){
                fc_ci(d2$raw, d1$raw, aggr = mean)
            }
        ),
        fc_med_ci = future_map2(
            data_co,
            data_tr,
            function(d1, d2){
                fc_ci(d2$raw, d1$raw, aggr = median)
            }
        ),
        fc_mean_lo = unlist(
            future_map_dbl(fc_mean_ci, function(x){x[1]})
        ),
        fc_mean_hi = unlist(
            future_map_dbl(fc_mean_ci, function(x){x[2]})
        ),
        fc_med_lo = unlist(
            future_map_dbl(fc_med_ci, function(x){x[1]})
        ),
        fc_med_hi = unlist(
            future_map_dbl(fc_med_ci, function(x){x[2]})
        ),
        medco = unlist(
            future_map_dbl(data_co, function(d){median(d$raw)})
        ),
        meanco = unlist(
            future_map_dbl(data_co, function(d){mean(d$raw)})
        ),
        medtr = unlist(
            future_map_dbl(data_tr, function(d){median(d$raw)})
        ),
        meantr = unlist(
            future_map_dbl(data_tr, function(d){mean(d$raw)})
        ),
        meantr_lo = meanco * fc2fc(fc_mean_lo),
        meantr_hi = meanco * fc2fc(fc_mean_hi),
        medtr_lo  = medco  * fc2fc(fc_med_lo),
        medtr_hi  = medco  * fc2fc(fc_med_hi),
        mean_out_ci = (
                meanco > meanco * fc2fc(fc_mean_hi) &
                meanco > meanco * fc2fc(fc_mean_lo)
            ) | (
                meanco < meanco * fc2fc(fc_mean_hi) &
                meanco < meanco * fc2fc(fc_mean_lo)
            ),
        med_out_ci = (
                medco > medco * fc2fc(fc_med_hi) &
                medco > medco * fc2fc(fc_med_lo)
            ) | (
                medco < medco * fc2fc(fc_med_hi) &
                medco < medco * fc2fc(fc_med_lo)
            )
    ) %>%
    select(-data_tr, -data_co, -fc_mean_ci, -fc_med_ci)

# tg <- (d %>% filter(tissue == 'TG BAT' & gs_mouse == 'Tnfa' & weeks == '4w'))$raw
# wt <- (d %>% filter(tissue == 'WT BAT' & gs_mouse == 'Tnfa' & weeks == '4w'))$raw


d <- d %>%
    mutate(
        tiss = gsub('.* ', '', tissue),
        gen  = gsub(' .*', '', tissue),
        gen  = factor(gen, levels = c('WT', 'KO', 'TG'), ordered = TRUE)
    ) %>%
    left_join(conds, by = c('gs_mouse', 'tissue', 'weeks'))

p <- ggplot(d, aes(x = weeks, y = zscore, color = gen)) +
    geom_boxplot(outlier.size = .5, outlier.shape = 18, lwd = .3) +
    geom_point(
        data = d,
        mapping = aes(
            y = 5,
            x = weeks,
            color = gen,
            shape = pval_t < 0.1
        ),
        size = 5,
        position = position_dodge(width = .5),
        show.legend = FALSE
    ) +
    scale_shape_manual(
        guide = guide_legend(title = 'Significant\ndifference\nto WT'),
        values = c(
            `TRUE`  = '*',
            `FALSE` = ''
        ),
        labels = c(
            `TRUE`  = 'T-test\np < 0.1',
            `FALSE` = ''
        )
    ) +
    scale_color_manual(
        guide = guide_legend(title = 'Bmp8b\ngenotype'),
        values = c(
            WT = '#44AA99',
            KO = '#999933',
            TG = '#CC6677'
        )
    ) +
    facet_grid(tiss ~ gs_mouse) +
    xlab('Age and Bmp8b genotype') +
    ylab('RNA level expression') +
    ylim(NA, 6) +
    theme_linedraw() +
    theme(
        text = element_text(family = 'DINPro'),
        panel.grid.major = element_line(color = '#CCCCCC'),
        panel.grid.minor = element_line(color = '#CCCCCC')
    )

ggsave('qpcr_zscore_ttests.pdf', device = cairo_pdf, width = 64, height = 7, limitsize = FALSE)

d %>%
    filter(tissue == 'TG BAT') %>%
    select(-other.x, -other.y) %>%
    group_by(gs_mouse, weeks, tissue) %>%
    summarize_all(first) %>%
    write_tsv('bat.qpcr.tests.tsv')
