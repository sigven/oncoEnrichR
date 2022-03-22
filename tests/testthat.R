library(testthat)
suppressPackageStartupMessages(library(oncoEnrichR))


load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))

log4r_logger <- log4r::logger(
  threshold = "INFO",
  appenders = log4r::console_appender(oncoEnrichR:::log4r_layout))

myc_data <- read.csv(system.file("extdata","myc_data.csv",
                                 package = "oncoEnrichR"),
                     stringsAsFactors = F) %>%
  dplyr::inner_join(oedb$genedb$all, by = "symbol") %>%
  dplyr::filter(!is.na(entrezgene))

background_sample_set <-
  oedb$genedb$all %>%
  dplyr::filter(gene_biotype == "protein-coding") %>%
  dplyr::sample_n(150)

testthat::test_check("oncoEnrichR")
