#!/usr/bin/env Rscript

# Dénes Türei EMBL 2018
# turei.denes@gmail.com

# generic
require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)

# bioconductor
require(limma)
require(Biobase)
require(HTSFilter)
require(edgeR)


rnaseq.main <- function(
        geo_id = 'GSE99412',
        biotype_file = 'mouse_biotypes.tsv',
        biotypes = c('protein_coding')
    ){
    
    data <- rnaseq.read_counts(
        geo_id,
        biotype_file = biotype_file,
        biotypes = biotypes
    )
    
    rnaseq.plot_counts_density(data)
    
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
        biotypes = NULL
    ){
    
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
        )
    
    counts_data <- counts_data %>%
        # remove duplicates
        group_by(name, sample) %>%
        summarize_all(first) %>%
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
            .
        )} %>%
        # transform to wide data frame
        spread(sample, count)
    
    return(
        list(
            samples = samples,
            data = counts_data,
            series = series_matrix
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
    
    d <- rnaseq.gather_data(data)
    
    p <- ggplot(d, aes(log2(count + 1))) +
        geom_density() +
        #scale_x_log10() +
        ylab('Frequency') +
        xlab('Count (log2)') +
        theme_minimal() +
        theme(text = element_text(family = 'DINPro'))
    
    ggsave('counts_density.pdf', device = cairo_pdf, width = 5, height = 4)
    
}


rnaseq.gather_data <- function(data){
    
    (
        data$data %>%
        gather(key = sample, value = count, -name, -biotype)
    )
    
}
