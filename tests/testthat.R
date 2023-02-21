library(testthat)
suppressPackageStartupMessages(library(oncoEnrichR))


load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))

# log4r_logger <- log4r::logger(
#   threshold = "INFO",
#   appenders = log4r::console_appender(oncoEnrichR:::log4r_layout))

myc_data <- read.csv(system.file("extdata","myc_data.csv",
                                 package = "oncoEnrichR"),
                     stringsAsFactors = F) %>%
  dplyr::inner_join(oedb$genedb$all, by = "symbol",
                    multiple = "all") %>%
  dplyr::filter(!is.na(entrezgene))

bg_set <-
  oedb[['genedb']][['all']] %>%
  dplyr::filter(.data$gene_biotype == "protein-coding") %>%
  dplyr::filter(!is.na(.data$entrezgene)) %>%
  dplyr::distinct()

background_sample <- bg_set %>%
  dplyr::sample_n(250)

background_sample_entrez <- as.character(background_sample$entrezgene)
background_full_entrez <- as.character(bg_set$entrezgene)

testthat::test_check("oncoEnrichR")
