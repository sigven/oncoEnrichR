tcga_oncoplot_genes <-
  function(qgenes,
           qsource = "symbol",
           cstrata = "site",
           genedb = "NULL",
           site = "Breast"){

  oncoEnrichR:::log4r_info(paste0("TCGA: generating oncoplot, tissue =  ", site))
  stopifnot(!is.null(genedb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == 'entrezgene'){
    query_genes_df <- data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(
      genedb, query_genes_df, by = "entrezgene") %>%
      dplyr::distinct()
  }else{
    query_genes_df <- dplyr::inner_join(
      genedb, query_genes_df, by = "symbol") %>%
      dplyr::distinct()
  }

  top_mutated_genes <- oncoEnrichR::tcgadb[['aberration']] %>%
    dplyr::filter(variant_type == "snv_indel" &
                    primary_site == site) %>%
    dplyr::filter(clinical_strata == cstrata) %>%
    dplyr::inner_join(dplyr::select(query_genes_df, symbol),
                      by = c("symbol")) %>%
    dplyr::select(symbol,
                  variant_type,
                  primary_site,
                  percent_mutated,
                  samples_mutated,
                  tot_samples,
                  percentile) %>%
    dplyr::rename(cohort_size = tot_samples) %>%
    dplyr::arrange(desc(percent_mutated)) %>%
    dplyr::distinct() %>%
    head(50)

  #n_omitted <- nrow(query_genes_df) - nrow(tcga_gene_stats)
  oncoEnrichR:::log4r_info(
    paste0("Choosing genes for oncoplot - highest SNV/InDel frequency in TCGA cohort - ", site))

  return(top_mutated_genes)

}


tcga_aberration_matrix <- function(qgenes,
                                 qsource = "symbol",
                                 cstrata = "site",
                                 vtype = "cna_ampl",
                                 genedb = NULL,
                                 percentile = FALSE){

  oncoEnrichR:::log4r_info(paste0("TCGA: generating gene aberration matrix, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == 'entrezgene'){
    query_genes_df <-
      data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "entrezgene") %>%
      dplyr::distinct()
  }else{
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "symbol") %>%
      dplyr::distinct()
  }

  title <- 'SNVs/InDels - TCGA'
  color <- 'steelblue'
  plotly_colors <- "Blues"
  if(vtype == 'cna_ampl'){
    plotly_colors <- "YlOrRd"
    title <- 'Copy number amplifications (sCNA) - TCGA'
    color <- 'darkgreen'
  }
  if(vtype == 'cna_homdel'){
    plotly_colors <- "YlGn"
    title <- 'Homozygous deletions (sCNA) - TCGA'
    color <- 'firebrick'
  }

  tcga_gene_stats <- oncoEnrichR::tcgadb[['aberration']] %>%
    dplyr::filter(clinical_strata == cstrata) %>%
    dplyr::filter(primary_site != "Other/Unknown") %>%
    dplyr::inner_join(
      dplyr::select(query_genes_df, symbol),
      by=c("symbol"))


  gene_candidates_init <- data.frame()
  tcga_ttypes <- sort(unique(tcga_gene_stats$primary_site))
  for(i in 1:length(sort(unique(tcga_gene_stats$symbol)))){
    init <- data.frame('primary_site' <- tcga_ttypes,
                       'primary_diagnosis_very_simplified' = NA,
                       'symbol' = sort(unique(tcga_gene_stats$symbol))[i],
                       'clinical_strata' = cstrata,
                       'percent_mutated' = 0,
                       'percentile' = 0,
                       'variant_type' = vtype,
                       'decile' = 0,
                       stringsAsFactors = F)
    colnames(init) <- c('primary_site',
                        'primary_diagnosis_very_simplified',
                        'symbol',
                        'clinical_strata',
                        'percent_mutated',
                        'percentile',
                        'variant_type',
                        'decile')
    gene_candidates_init <- rbind(gene_candidates_init, init)
    i <- i + 1
  }
  gene_candidates_init <- gene_candidates_init %>%
    dplyr::filter(primary_site != 'Other/Unknown' &
                    primary_site != "Pancancer")

  gene_aberrations <- tcga_gene_stats %>%
    dplyr::filter(variant_type == vtype) %>%
    dplyr::filter(primary_site != "Pancancer" &
                    primary_site != "Other/Unknown")


  site_stats_zero <- tcga_gene_stats %>%
    dplyr::select(primary_site,tot_samples) %>%
    dplyr::distinct() %>%
    dplyr::mutate(samples_mutated = 0)

  pancan_order <- tcga_gene_stats %>%
    dplyr::filter(variant_type == vtype & primary_site == "Pancancer") %>%
    dplyr::mutate(pancancer_percent_mutated = percent_mutated) %>%
    dplyr::mutate(pancancer_percentile = percentile) %>%
    dplyr::select(symbol, pancancer_percent_mutated, pancancer_percentile)

  zero_frequency_genes <-
    dplyr::anti_join(gene_candidates_init, gene_aberrations,
                     by = c("symbol", "primary_site", "variant_type")) %>%
    dplyr::left_join(site_stats_zero,by=c("primary_site"))

  gene_aberrations <-
    dplyr::left_join(dplyr::bind_rows(gene_aberrations, zero_frequency_genes),
                     pancan_order, by = c("symbol")) %>%
    dplyr::mutate(
      pancancer_percent_mutated =
        dplyr::if_else(is.na(pancancer_percent_mutated),
                       as.numeric(0),
                       as.numeric(pancancer_percent_mutated)))


  gene_aberrations <- gene_aberrations %>%
    dplyr::arrange(pancancer_percent_mutated, primary_site)

  top_mutated <- gene_aberrations %>%
    dplyr::arrange(desc(pancancer_percent_mutated)) %>%
    dplyr::select(symbol) %>%
    dplyr::distinct() %>%
    head(75)

  gene_aberrations_top <- gene_aberrations %>%
    dplyr::inner_join(top_mutated, by = "symbol") %>%
    dplyr::mutate(symbol = factor(symbol, unique(symbol)))


  gene_aberration_top_mat <- as.data.frame(
    gene_aberrations_top %>%
    dplyr::select(symbol, primary_site, percent_mutated) %>%
    tidyr::pivot_wider(names_from = primary_site,
                       values_from = percent_mutated)
  )
  rownames(gene_aberration_top_mat) <-
    gene_aberration_top_mat$symbol
  gene_aberration_top_mat$symbol <- NULL
  gene_aberration_top_mat <- as.matrix(gene_aberration_top_mat)

  return(gene_aberration_top_mat)

}

tcga_aberration_table <- function(qgenes,
                                  qsource = "entrezgene",
                                  genedb = NULL,
                                  vtype = "snv_indel"){

  oncoEnrichR:::log4r_info(paste0("TCGA: collecting gene aberration data table, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == "entrezgene"){
    query_genes_df <- data.frame(entrezgene = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "entrezgene") %>% dplyr::distinct()
  }else{
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "symbol") %>% dplyr::distinct()
  }

  aberration_data <- oncoEnrichR::tcgadb[['aberration']] %>%
    dplyr::filter(clinical_strata == "site_diagnosis" &
                    variant_type == vtype &
                    primary_site != "Pancancer") %>%
    dplyr::select(symbol,
                  primary_site,
                  primary_diagnosis_very_simplified,
                  variant_type,
                  samples_mutated,
                  tot_samples,
                  percent_mutated,
                  percentile) %>%
    dplyr::rename(primary_diagnosis =
                    primary_diagnosis_very_simplified,
                  cohort_size = tot_samples) %>%
    dplyr::filter(!stringr::str_detect(primary_diagnosis,"^Other")) %>%
    dplyr::left_join(dplyr::select(genedb, symbol, entrezgene),by=c("symbol")) %>%
    dplyr::inner_join(dplyr::select(query_genes_df, entrezgene),by=c("entrezgene")) %>%
    dplyr::distinct() %>%
    dplyr::mutate(gene = paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                                entrezgene,"' target='_blank'>",symbol,"</a>")) %>%
    dplyr::select(-c(entrezgene,symbol)) %>%
    dplyr::select(gene, variant_type,
                  primary_site,
                  primary_diagnosis,
                  percent_mutated,
                  samples_mutated,
                  cohort_size,
                  percentile,
                  dplyr::everything())

  return(aberration_data)
}

tcga_co_expression <- function(qgenes, qsource = "symbol", genedb = NULL){

  oncoEnrichR:::log4r_info("TCGA: collecting co-expression data (strong negative and positive correlations)")
  stopifnot(!is.null(genedb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == "entrezgene"){
    query_genes_df <- data.frame(entrezgene = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(dplyr::select(genedb, entrezgene, symbol),
                                        query_genes_df, by = "entrezgene") %>%
      dplyr::distinct()
  }else{
    query_genes_df <- dplyr::inner_join(dplyr::select(genedb, entrezgene, symbol),
                                        query_genes_df, by = "symbol") %>%
      dplyr::distinct()
  }

  coexp_target_1 <- oncoEnrichR::tcgadb[['co_expression']] %>%
    dplyr::mutate(corrtype = dplyr::if_else(
      r < 0,
      "Negative",
      as.character("Positive"))) %>%
    dplyr::mutate(correlation = dplyr::case_when(
      r < -0.6 & r >= -0.8 ~"Strong negative",
      r < -0.8 ~ "Very strong negative",
      r >= 0.6 & r < 0.8 ~"Strong positive",
      r >= 0.8 ~ "Very strong positive",
      TRUE ~ as.character(NA))) %>%
    dplyr::select(symbol, symbol_partner,
                  corrtype, correlation, r, p_value, tumor) %>%
    dplyr::left_join(query_genes_df, by = c("symbol" = "symbol")) %>%
    dplyr::filter(!is.na(entrezgene))

  coexp_target_tcga <- oncoEnrichR::tcgadb[['co_expression']] %>%
    dplyr::mutate(corrtype = dplyr::if_else(
      r < 0,
      "Negative",
      as.character("Positive"))) %>%
    dplyr::mutate(correlation = dplyr::case_when(
      r < -0.6 & r >= -0.8 ~"Strong negative",
      r < -0.8 ~ "Very strong negative",
      r >= 0.6 & r < 0.8 ~"Strong positive",
      r >= 0.8 ~ "Very strong positive",
      TRUE ~ as.character(NA))) %>%
    dplyr::select(symbol, symbol_partner,
                  corrtype, correlation, r, p_value, tumor) %>%
    dplyr::left_join(query_genes_df, by = c("symbol_partner" = "symbol")) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::mutate(tmp = symbol) %>%
    dplyr::mutate(symbol = symbol_partner) %>%
    dplyr::mutate(symbol_partner = tmp) %>%
    dplyr::select(-tmp) %>%
    dplyr::bind_rows(coexp_target_1) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dplyr::select(genedb, name,
                                   oncogene, cancer_driver,
                                   tumor_suppressor,
                                   symbol, ot_tractability_compound),
                     by = c("symbol_partner" = "symbol")) %>%
    dplyr::rename(target_gene = symbol,
                  partner_gene = symbol_partner,
                  partner_genename = name,
                  target_tractability = ot_tractability_compound) %>%
    dplyr::mutate(r = as.numeric(round(r, digits = 3))) %>%
    dplyr::filter(stringr::str_detect(tumor,"BRCA|LUAD|SKCM|COAD|SARC|PRAD|ESCA|MESO|UCEC|OV|CHOL|THCA|COAD|BLCA|STAD|KIRP|GBM|HNSC")) %>%
    dplyr::mutate(primary_site = "Breast") %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "LUAD" | tumor == "LUSC","Lung", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "STAD" | tumor == "ESCA","Esophagus/Stomach", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "BLCA","Bladder", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "SARC","Soft Tissue", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "HNSC","Head and Neck", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "UCEC","Endometrium", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "LAML","Myeloid", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "SKCM","Skin", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "GBM","Brain", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "PAAD","Pancreas", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "COAD" | tumor == "READ","Colon/Rectum", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "PRAD","Prostate", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "THCA","Thyroid", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "KIRP" | tumor == "KICH","Kidney", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "OV","Ovary", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "MESO","Pleura", as.character(primary_site))) %>%
    dplyr::mutate(primary_site = dplyr::if_else(tumor == "CHOL","Biliary Tract", as.character(primary_site)))



  coexp_target_tcga <- coexp_target_tcga %>%
    dplyr::filter(stringr::str_detect(correlation,"Very") | !is.na(oncogene) | !is.na(tumor_suppressor) | !is.na(cancer_driver)) %>%
    dplyr::arrange(desc(r))

  coexp_target_tcga_positive <-
    dplyr::filter(coexp_target_tcga, corrtype == "Positive")
  coexp_target_tcga_negative <-
    dplyr::filter(coexp_target_tcga, corrtype == "Negative") %>%
    dplyr::arrange(r)

  coexp_target_tcga <-
    dplyr::bind_rows(coexp_target_tcga_negative,
                     coexp_target_tcga_positive)


  return(coexp_target_tcga)

}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

