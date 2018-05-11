#!/usr/bin/env Rscript

# Dénes Türei EMBL 2018
# turei.denes@gmail.com

# generic
require(readr)
require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)

# bioconductor
require(limma)
require(Biobase)
require(HTSFilter)
require(edgeR)
require(org.Mm.eg.db)


rnaseq.main <- function(
        geo_id = 'GSE99412',
        biotype_file = 'mouse_biotypes.tsv',
        biotypes = c('protein_coding'),
        use_samples = c('BAT_WT', 'BAT_UCP'),
        design_formula = ~ 0 + data$samples$sample_name,
        etpair = c('BAT_WT', 'BAT_UCP'),
        do_lowexp_filter = TRUE,
        do_boxplots = TRUE,
        do_enrichment = TRUE
    ){
    
    data <- rnaseq.read_counts(
        geo_id,
        biotype_file = biotype_file,
        biotypes = biotypes,
        use_samples = use_samples
    ) %>%
    {`if`(
        do_lowexp_filter,
        rnaseq.plot_counts_density(.) %>%
        rnaseq.get_design_matrix(form = design_formula) %>%
        rnaseq.get_dgelist() %>%
        rnaseq.cpm('dgelist', 'raw') %>%
        rnaseq.plot_raw() %>%
        rnaseq.remove_50('dgelist', 'c0', '') %>%
        rnaseq.remove_50('dgelist', 'cpm1', '') %>%
        rnaseq.htsf_rm50(var = 'c0') %>%
        rnaseq.htsf_rm50(var = 'cpm1') %>%
        rnaseq.htsf_raw() %>%
        rnaseq.htsf_after(htsfvar = 'htsf_rm50c0') %>%
        rnaseq.htsf_after(htsfvar = 'htsf_rm50cpm1') %>%
        rnaseq.htsf_after(htsfvar = 'htsf_raw') %>%
        rnaseq.remove_50(
            'dgelist.htsf_raw', 'c0', 'after.htsf.raw.') %>%
        rnaseq.remove_50(
            'dgelist.htsf_raw', 'cpm1', 'after.htsf.raw.') %>%
        rnaseq.remove_50(
            'dgelist.htsf_rm50c0', 'c0', 'after.htsf.rm50c0.') %>%
        rnaseq.remove_50(
            'dgelist.htsf_rm50c0', 'cpm1', 'after.htsf.rm50c0.') %>%
        rnaseq.remove_50(
            'dgelist.htsf_rm50cpm1', 'c0', 'after.htsf.rm50cpm1-') %>%
        rnaseq.remove_50(
            'dgelist.htsf_rm50cpm1', 'cpm1', 'after.htsf.rm50cpm1.') %>%
        rnaseq.select_lowexp_removed('htsf_rm50cpm1'),
        .
    )} %>%
    {`if`(
        do_boxplots,
        rnaseq.samples_boxplot_logcounts(.) %>%
        rnaseq.samples_boxplot_logcpmnormalized() %>%
        rnaseq.estimate_disp() %>%
        rnaseq.glm_fit() %>%
        rnaseq.plot_mds() %>%
        rnaseq.plot_bcv() %>%
        rnaseq.exact_test(pair = etpair) %>%
        rnaseq.plot_md(),
        .
    )} %>%
    {`if`(
        do_enrichment,
        rnaseq.map_entrez(var = 'et') %>%
        rnaseq.goenrich() %>%
        rnaseq.keggenrich(),
        .
    )}
    #rnaseq.hts_filter2() %>%
    #rnaseq.hts_filter3()
    
    return(data)
    
}


rnaseq.read_counts <- function(
        geo_id,
        only_samples = NULL,
        exclude_samples = NULL,
        fname_end = 'counts.tsv.gz',
        download_method = 'wget',
        datadir = 'data',
        biotype_file = NULL,
        biotypes = NULL,
        use_samples = NULL,
        sample_filter = sample_name
    ){
    
    sample_filter <- enquo(sample_filter)
    
    samples_label <- `if`(
        is.null(use_samples),
        'all',
        paste0(use_samples, collapse = '-')
    )
    
    # check for data dir, create if missing
    if(!dir.exists(datadir)){dir.create('data')}
    
    # check for GEO series file, download and extract TAR if missing
    tarname <- file.path(datadir, sprintf('%s_RAW.tar', geo_id))
    if(!file.exists(tarname)){
        geo_url <- sprintf(
            'https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file',
            geo_id
        )
        download.file(
            geo_url,
            destfile = tarname,
            method = download_method
        )
        untar(tarname, exdir = datadir, compressed = FALSE, tar = 'internal')
    }
    
    # list files in data dir
    files <- Map(
        # prepend data dir path for each file
        function(fname){
            file.path(datadir, fname)
        },
        Filter(
            # only files with specified ending
            function(fname){
                # GEO ID of sample
                sample <- gsub('.*/', '', gsub('_.*', '', fname))
                
                endsWith(fname, fname_end) & (
                    is.null(only_samples) | sample %in% only_samples
                ) & (
                    is.null(exclude_samples) | !(sample %in% exclude_samples)
                )
            },
            list.files(path = datadir)
        )
    )
    
    # check for series matrix, download if missing
    series_matrix_txt <- sprintf('%s_series_matrix.txt.gz', geo_id)
    series_matrix_name <- file.path(datadir, series_matrix_txt)
    if(!file.exists(series_matrix_name)){
        series_matrix_url <- sprintf(
            'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/%s/matrix/%s',
            geo_id,
            series_matrix_txt
        )
        download.file(
            series_matrix_url,
            destfile = series_matrix_name,
            method = download_method
        )
    }
    series_matrix <- rnaseq.read_series_matrix(series_matrix_name)
    
    message(' > Reading read counts data')
    
    # building main data frame
    counts_data <-
        # data frame with all read counts
        bind_rows(
            lapply(
                files,
                function(fname){
                    # read read count files
                    suppressMessages(
                        read_tsv(fname, col_names = c('name', 'count'))
                    ) %>%
                    # add sample ID
                    mutate(sample = gsub('.*/', '', gsub('_.*', '', fname)))
                }
            )
        )
    
    # create a data frame with sample names and IDs only
    samples <- counts_data %>%
        select(sample) %>%
        group_by(sample) %>%
        summarize_all(first) %>%
        # add annotations from the series matrix,
        # primarily sample names
        left_join(
            series_matrix$samples %>%
                mutate(Sample_title = gsub(' \\[.*', '', Sample_title)) %>%
                select(
                    sample = Sample_geo_accession,
                    sample_name = Sample_title
                ),
            by = c('sample')
        ) %>%
        mutate(
            tissue  = gsub('_.*', '', sample_name),
            ko_gene = gsub('.*_', '', sample_name)
        ) %>%
        mutate(
            tissue_type = ifelse(
                tissue %in% c('BAT', 'iWAT'),
                'adipose',
                ifelse(
                    tissue == 'Hypothalamus',
                    'nervous',
                    'muscle'
                )
            )
        )
    
    counts_data <- counts_data %>%
        # remove duplicates
        group_by(name, sample) %>%
        summarize_all(first) %>%
        ungroup() %>%
        # add biotypes data if available
        {`if`(
            !is.null(biotype_file),
            left_join(
                .,
                rnaseq.read_biotypes(biotype_file),
                by = c('name')
            ) %>%
            # filter for biotypes if required
            {`if`(
                !is.null(biotypes),
                filter(., biotype %in% biotypes),
                .
            )},
            # add a dummy column isntead of biotype
            # for smooth further processing
            mutate(., biotype = 'Unknown')
        )} %>%
        # keep only certain samples
        {`if`(
            !is.null(use_samples),
            left_join(., samples, by = c('sample')) %>%
            filter(!!sample_filter %in% use_samples) %>%
            select(-sample_name, -tissue, -ko_gene, -tissue_type),
            .
        )} %>%
        # transform to wide data frame
        spread(sample, count)
    
    samples <- samples %>%
        {`if`(
            !is.null(use_samples),
            filter(., !!sample_filter %in% use_samples),
            .
        )}
    
    return(
        list(
            samples = samples,
            data = counts_data,
            series = series_matrix,
            geo_id = geo_id,
            samples_label = samples_label,
            prefix = sprintf('%s_%s', geo_id, samples_label)
        )
    )
    
}


rnaseq.read_series_matrix <- function(fname){
    
    con <- file(fname, 'r')
    lns <- readLines(con)
    close.connection(con)
    
    data     <- list()
    sampledf <- list()
    
    for(l in lns){
        
        if(nchar(l)){
            lcon <- textConnection(l, open = 'r')
            fields <- as.character(
                read.delim(
                    lcon,
                    header = FALSE,
                    sep = '\t',
                    stringsAsFactors = FALSE
                )[1,]
            )
            if(length(fields) > 1){
                key <- gsub('^!', '', as.character(fields[1]))
                values <- fields[2:length(fields)]
                if(length(values) == 1){
                    data[[key]] <- values
                }else{
                    sampledf[[key]] <- values
                }
            }
        }
        
        data[['samples']] <- as.data.frame(sampledf)
        
    }
    
    return(data)
    
}


rnaseq.read_biotypes <- function(fname){
    
    (
        suppressMessages(read_tsv(fname)) %>%
        group_by(name, biotype) %>%
        summarize_all(first) %>%
        ungroup()
    )
    
}


rnaseq.export_biotype_stats <- function(data, outfile = 'biotypes.tsv'){
    
    invisible(return(
        
        data$data %>%
        group_by(biotype) %>%
        mutate(cnt = n()) %>%
        summarize_all(first) %>%
        select(biotype, cnt) %>%
        write_tsv(outfile)
        
    ))
    
}


rnaseq.plot_counts_density <- function(data){
    
    pdfname <- sprintf('%s_counts_density.pdf', data$prefix)
    
    d <- rnaseq.gather_data(data)
    
    message(sprintf(' > Plotting counts density into `%s`', pdfname))
    
    p <- ggplot(d, aes(log2(count + 1))) +
        geom_density() +
        #scale_x_log10() +
        ylab('Frequency') +
        xlab('Count (log2)') +
        theme_minimal() +
        theme(text = element_text(family = 'DINPro'))
    
    ggsave(pdfname, device = cairo_pdf, width = 5, height = 4)
    
    return(data)
    
}


rnaseq.gather_data <- function(data){
    
    message(' > Gather')
    
    (
        data$data %>%
        gather(key = sample, value = count, -name, -biotype)
    )
    
}


rnaseq.get_design_matrix <- function(
        data,
        form = ~ 0 + data$samples$sample_name
    ){
    
    message(' > Creating design matrix')
    
    data$design_matrix <- model.matrix(form)
    
    return(data)
    
}


rnaseq.get_dgelist <- function(data){
    
    message(' > Creating DGElist object')
    
    nm <- data$data$name
    
    countdf <- data$data %>%
        # just to be safe having no groups
        ungroup() %>%
        select(-biotype, -name) %>%
        as.data.frame()
    
    rownames(countdf) <- nm
    
    data$dgelist <- DGEList(
        counts = countdf,
        group  = data$samples$sample_name
    )
    
    data$dgelist <- calcNormFactors(data$dgelist)
    
    return(data)
    
}


rnaseq.cpm <- function(data, dgelistvar, var){
    
    message(sprintf(' > Calculating CPMs and TMMs for `%s`', var))
    
    data[[sprintf('cpm.log2.%s', var)]]     <-
        cpm(data[[dgelistvar]], normalized.lib.sizes = FALSE, log = TRUE)
    data[[sprintf('cpm.log2.tmm.%s', var)]] <-
        cpm(data[[dgelistvar]], normalized.lib.sizes = TRUE,  log = TRUE)
    data[[sprintf('cpm.tmm.%s', var)]]      <-
        cpm(data[[dgelistvar]], normalized.lib.sizes = TRUE,  log = FALSE)
    
    return(data)
    
}


rnaseq.plot_raw <- function(data){
    
    pdfname <- sprintf('%s_raw.pdf', data$prefix)
    
    message(sprintf(' > Plotting raw data into `%s`', pdfname))
    
    cairo_pdf(pdfname, width = 10, height = 5)
        
        par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
        
        data$voom <- voom(data$dgelist, data$design_matrix, plot = TRUE)
        
        rnaseq.plot_histogram(
            data$dgelist$counts,
            'Raw data'
        )
        
    dev.off()
    
    return(data)
    
}


rnaseq.remove_50 <- function(data, dgelistvar, label, var){
    
    pdfname <- sprintf('%s_%s_rm50%s.pdf', data$prefix, var, label)
    above <- `if`(label == 'cpm1', 'CPM 1', 'count 0')
    
    message(
        sprintf(' > Remove 50 above %s; plotting into `%s`', above, pdfname)
    )
    
    dgelist <- data[[dgelistvar]]
    rmvar   <- sprintf('%srm50%s', var, label)
    # transcripts detected in more than the 50% of the samples
    keep <- `if`(
        label == 'cpm1',
        rowSums(cpm(dgelist) > 1) >= ncol(dgelist) / 2,
        rowSums(dgelist$counts > 0) >= ncol(dgelist) / 2
    )
    data[[rmvar]] <- dgelist[keep, , keep.lib.sizes = TRUE]
    
    data <- rnaseq.cpm(data, rmvar, rmvar)
    
    cairo_pdf(pdfname, width = 10, height = 5)
        
        par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
        
        data[[sprintf('voom.%s', rmvar)]] <- voom(
            data[[rmvar]],
            data$design_matrix,
            plot = TRUE
        )
        
        rnaseq.plot_histogram(
            data[[rmvar]]$counts,
            main = sprintf(
                '%s in < 50%% of samples removed',
                above
            )
        )
        
    dev.off()
    
    return(data)
    
}


rnaseq.htsf_rm50 <- function(data, var){
    
    varname <- ifelse(var == 'c0', 'count', 'CPM')
    pdfname <- sprintf('%s_remove_50_%s_HTSF.pdf', data$prefix, var)
    taskname <- sprintf(
        'HTSF on %s %i > 50%%', varname, `if`(var == 'c0', 0, 1)
    )
    htsfvar <- sprintf('htsf_rm50%s', var)
    counts <- data[[sprintf('rm50%s', var)]]$counts
    main <- sprintf(
        'HTSF on %s > %i > 50%%',
        varname,
        `if`(var == 'c0', 0, 1)
    )
    
    data <- rnaseq.htsf_generic(data, counts, htsfvar, pdfname, taskname, main)
    
    return(data)
    
}


rnaseq.htsf_raw <- function(data){
    
    pdfname <- sprintf('%s_raw_HTSF.pdf', data$prefix)
    taskname <- 'HTSF on raw counts'
    htsfvar <- 'htsf_raw'
    counts <- data$dgelist$counts
    main <- 'HTSF on raw counts'
    
    data <- rnaseq.htsf_generic(data, counts, htsfvar, pdfname, taskname, main)
    
    return(data)
    
}


rnaseq.htsf_after <- function(data, htsfvar){
    
    pdfname <- sprintf('%s_%s_after_HTSF.pdf', data$prefix, htsfvar)
    taskname <- sprintf('HTSF after HTSF on `%s`', htsfvar)
    htsfvar_after <- sprintf('%s_after', htsfvar)
    counts <- data[[htsfvar]]
    main <- taskname
    
    data <- rnaseq.htsf_generic(
        data, counts, htsfvar_after, pdfname, taskname, main
    )
    
    return(data)
    
}


rnaseq.htsf_generic <- function(
        data, counts, htsfvar, pdfname, taskname, main
    ){
    
    message(sprintf(' > %s; plotting into `%s`', taskname, pdfname))
    
    cairo_pdf(pdfname, width = 15, height = 5)
        
        par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
        
        data[[htsfvar]] <- HTSFilter(
            counts,
            data$samples$sample_name
        )$filteredData
        
        rnaseq.plot_histogram(
            data[[htsfvar]],
            main = main,
            logg = log2
        )
        
        dgelistvar <- sprintf('dgelist.%s', htsfvar)
        data[[dgelistvar]] <- DGEList(
            data[[htsfvar]], group = data$samples$sample_name
        )
        data[[dgelistvar]] <- calcNormFactors(
            data[[dgelistvar]]
        )
        data[[sprintf('voom.%s', htsfvar)]] <- voom(
            data[[dgelistvar]],
            data$design_matrix,
            plot=TRUE
        )
        data <- rnaseq.cpm(data, dgelistvar, htsfvar)
        
    dev.off()
    
    return(data)
    
}


rnaseq.plot_histogram <- function(x, main = '', logg = log){
    
    hist(
        logg(x + 1),
        col = 'grey',
        breaks = 25,
        main = main,
        xlab = 'Log(counts+1)'
    )
    
}


rnaseq.select_lowexp_removed <- function(data, label){
    
    # for GSE99412 I decided to use label = 'dgelist.htsf_rm50cpm1'
    
    data$rmlowexp <- data[[sprintf('dgelist.%s', label)]]
    
    return(data)
    
}


rnaseq.samples_boxplot_logcounts <- function(data){
    
    return(
        rnaseq.samples_boxplot(data, data$rmlowexp$counts, 'logcounts')
    )
    
}


rnaseq.samples_boxplot_logcpmnormalized <- function(data){
    
    return(
        rnaseq.samples_boxplot(
            data,
            cpm(data$rmlowexp, normalized.lib.sizes=TRUE, log = FALSE),
            'logcpmnormalized'
        )
    )
    
}


rnaseq.samples_boxplot <- function(data, counts, label = '0'){
    
    pdfname <- sprintf('%s_samples_boxplot_%s.pdf', data$prefix, label)
    
    message(sprintf(' > Plotting samples boxplot into `%s`', pdfname))
    
    d <- gather(
            as.data.frame(counts),
            sample, count
        ) %>%
        left_join(data$samples, by = c('sample'))
    
    p <- ggplot(d, aes(y = count, x = sample, color = sample_name)) +
        geom_boxplot() +
        scale_y_log10() +
        scale_color_discrete(guide = guide_legend(title = 'Condition')) +
        xlab('Samples') +
        ylab('Read counts (log)') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text.x = element_text(
                size = 7,
                angle = 90,
                vjust = 0.5
            )
        )
    
    ggsave(pdfname, device = cairo_pdf, width = 10, height = 5)
    
    return(invisible(data))
    
}


rnaseq.estimate_disp <- function(data, label = 'rmlowexp'){
    
    message(
        sprintf(' > Estimating dispersion on expression list `%s`', label)
    )
    
    data[[label]] <- estimateDisp(data[[label]], design = data$design_matrix)
    
    return(data)
    
}


rnaseq.glm_fit <- function(data, label = 'rmlowexp'){
    
    message(
        sprintf(' > Fitting GLM on expression list `%s`', label)
    )
    
    glmfitlabel <- sprintf('glm.fit_%s', label)
    ftestlabel  <- sprintf('ftest_%s', label)
    
    data[[glmfitlabel]] <- glmQLFit(data[[label]], data$desig_matrix)
    #data[[ftestlabel]]  <- glmQLFTest(data[[glmfitlabel]], coef = 2)
    
    return(data)
    
}


rnaseq.plot_mds <- function(data, label = 'rmlowexp', tissues = NULL){
    
    tissue_lab <- `if`(
        is.null(tissues),
        'all',
        paste(tissues, sep = '-', collapse = '-')
    )
    
    pdfname <- sprintf('%s_mds_%s_%s.pdf', data$prefix, label, tissue_lab)
    
    message(sprintf(
        ' > Plotting MDS of expression list `%s`; plotting into `%s`',
        label, pdfname
    ))
    
    mdslab <- sprintf('mds_%s', label)
    data[[mdslab]] <- plotMDS(data[[label]], plot = FALSE)
    
    p <- ggplot(
            data.frame(
                sample = data$samples$sample,
                mdsx = data[[mdslab]]$x,
                mdsy = data[[mdslab]]$y
            ) %>%
            left_join(data$samples, by = c('sample')) %>%
            {`if`(
                is.null(tissues),
                .,
                filter(., gsub('_.*', '', sample_name) %in% tissues)
            )},
            aes(x = mdsx, y = mdsy, color = sample_name)
        ) +
        geom_point() +
        scale_color_brewer(
            guide = guide_legend(title = 'Conditions'),
            palette = 'Set1'
        ) +
        xlab('MDS1') +
        ylab('MDS2') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro')
        )
    
    ggsave(pdfname, device = cairo_pdf, width = 6, height = 6)
    
    return(data)
    
}


rnaseq.plot_bcv <- function(data){
    
    pdfname <- sprintf('%s_bcv.pdf', data$prefix)
    
    cairo_pdf(pdfname, width = 5, height = 5)
        
        plotBCV(data$rmlowexp)
        
    dev.off()
    
    return(data)
    
}


rnaseq.exact_test <- function(data, pair = 1:2){
    
    message(' > Performing exact tests')
    
    data$et <- exactTest(data$rmlowexp, pair = pair)
    
    return(data)
    
}


rnaseq.plot_md <- function(data){
    
    message(' > Plotting MD')
    
    plotMD(data$et)
    
    return(data)
    
}


rnaseq.map_entrez <- function(data, var, sub = 'table'){
    
    message(' > Mapping to Entrez')
    
    entrezvar <- sprintf('%s_entrez', var)
    
    entrez <- unlist(map(
        mget(
            rownames(data[[var]][[sub]]),
            org.Mm.eg.db::org.Mm.egSYMBOL2EG,
            ifnotfound = NA
        ),
        `[`,
        1
    ))
    
    data[[entrezvar]] <- data[[var]]
    data[[entrezvar]][[sub]] <-
        data[[entrezvar]][[sub]][which(!is.na(entrez)),]
    rownames(data[[entrezvar]][[sub]]) <- entrez[which(!is.na(entrez))]
    
    return(data)
    
}


rnaseq.goenrich <- funciton(data, var = 'et_entrez'){
    
    govar <- sprintf('go_%s', var)
    
    message(' > GO enrichment analysis')
    
    data[[govar]] <- goana(data[[var]], species = 'Mm')
    
    return(data)
    
}


rnaseq.keggenrich <- funciton(data, var = 'et_entrez'){
    
    keggvar <- sprintf('kegg_%s', var)
    
    message(' > KEGG pathway enrichment analysis')
    
    data[[keggvar]] <- keggaa(data[[var]], species = 'Mm')
    
    return(data)
    
}


rnaseq.hts_filter3 <-function(data){
    
    message(' > HTS filter with GLM fit')
    
    data$dgelist <- estimateDisp(data$dgelist, design = data$design_matrix)
    data$fit <- glmFit(data$dgelist, design = data$design_matrix)
    data$lrt <- glmLRT(data$fit, coef = 2)
    
    cairo_pdf('hts3.pdf', width = 6, height = 6)
        
        par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
        
        data$dgelist.hts3 <- HTSFilter(
            data$lrt,
            DGEGLM = data$fit,
            s.len=25,
            plot = TRUE
        )$filteredData
        
    dev.off()
    
    return(data)
    
}
