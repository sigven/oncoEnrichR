tcga_oncoplot_genes <-
  function(qgenes,
           qsource = "symbol",
           cstrata = "site",
           genedb = NULL,
           tcgadb = NULL,
           site = "Breast") {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    lgr::lgr$info( paste0("TCGA: generating oncoplot, tissue =  ", site))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(tcgadb))
    stopifnot(
      identical(names(tcgadb),
                c("coexpression",
                  "aberration",
                  "recurrent_variants",
                  "median_ttype_expression",
                  "pfam",
                  "maf_codes",
                  "maf",
                  "site_code",
                  "diagnosis_code",
                  "clinical_strata_code")
      ))
    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
    validate_db_df(tcgadb$site_code, dbtype = "tcga_site_code")
    validate_db_df(tcgadb$diagnosis_code, dbtype = "tcga_diagnosis_code")
    validate_db_df(tcgadb$clinical_strata_code, dbtype = "tcga_clinical_strata_code")
    stopifnot(site %in% unique(tcgadb$site_code$primary_site))
    stopifnot(qsource == "symbol" | qsource == "entrezgene")
    stopifnot(cstrata == "site" | cstrata == "site_diagnosis")
    query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
    if (qsource == 'entrezgene') {
      stopifnot(is.integer(qgenes))
      query_genes_df <- data.frame('entrezgene' = qgenes, stringsAsFactors = F)
      query_genes_df <- dplyr::inner_join(
        genedb, query_genes_df, by = "entrezgene", multiple = "all") |>
        dplyr::distinct()
    } else {
      stopifnot(is.character(qgenes))
      query_genes_df <- dplyr::inner_join(
        genedb, query_genes_df, by = "symbol", multiple = "all") |>
        dplyr::distinct()
    }

    top_mutated_genes <- tcgadb[['aberration']] |>
      dplyr::inner_join(
        dplyr::select(query_genes_df, "symbol"),
        by = c("symbol"), multiple = "all") |>
      dplyr::left_join(tcgadb[['site_code']],
                       by = "site_code", multiple = "all") |>
      dplyr::left_join(tcgadb[['diagnosis_code']],
                       by = "diagnosis_code", multiple = "all") |>
      dplyr::left_join(tcgadb[['clinical_strata_code']],
                       by = "clinical_strata_code", multiple = "all") |>
      dplyr::select(-c("site_code",
                       "diagnosis_code",
                       "clinical_strata_code")) |>
      dplyr::filter(.data$variant_type == "snv_indel" &
                      .data$primary_site == site) |>
      dplyr::filter(.data$clinical_strata == cstrata) |>

      dplyr::select(c("symbol",
                    "variant_type",
                    "primary_site",
                    "percent_mutated",
                    "samples_mutated",
                    "tot_samples",
                    "percentile")) |>
      dplyr::rename(cohort_size = "tot_samples") |>
      dplyr::arrange(dplyr::desc(.data$percent_mutated)) |>
      dplyr::distinct() |>
      utils::head(50)

    #n_omitted <- nrow(query_genes_df) - nrow(tcga_gene_stats)
    lgr::lgr$info(
      paste0("Choosing genes for oncoplot - highest SNV/InDel frequency in TCGA cohort - ", site))

    return(top_mutated_genes)

  }


tcga_aberration_matrix <- function(qgenes,
                                 qsource = "symbol",
                                 cstrata = "site",
                                 vtype = "cna_ampl",
                                 genedb = NULL,
                                 tcgadb = NULL,
                                 percentile = FALSE) {
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  lgr::lgr$info( paste0("TCGA: generating gene aberration matrix, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(tcgadb))
  stopifnot(
    identical(names(tcgadb),
              c("coexpression",
                "aberration",
                "recurrent_variants",
                "median_ttype_expression",
                "pfam",
                "maf_codes",
                "maf",
                "site_code",
                "diagnosis_code",
                "clinical_strata_code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
  validate_db_df(tcgadb$site_code, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$diagnosis_code, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$clinical_strata_code, dbtype = "tcga_clinical_strata_code")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(cstrata == "site")
  stopifnot(vtype == "cna_ampl" | vtype == "cna_homdel")
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if (qsource == 'entrezgene') {
    stopifnot(is.integer(qgenes))
    query_genes_df <-
      data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "entrezgene", multiple = "all") |>
      dplyr::distinct()
  } else {
    stopifnot(is.character(qgenes))
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "symbol", multiple = "all") |>
      dplyr::distinct()
  }

  title <- 'SNVs/InDels - TCGA'
  color <- 'steelblue'
  plotly_colors <- "Blues"
  if (vtype == 'cna_ampl') {
    plotly_colors <- "YlOrRd"
    title <- 'Copy number amplifications (sCNA) - TCGA'
    color <- 'darkgreen'
  }
  if (vtype == 'cna_homdel') {
    plotly_colors <- "YlGn"
    title <- 'Homozygous deletions (sCNA) - TCGA'
    color <- 'firebrick'
  }

  tcga_gene_stats <- tcgadb[['aberration']] |>
    dplyr::inner_join(
      dplyr::select(query_genes_df, c("symbol")),
      by = c("symbol"), multiple = "all"
    )

  ## return NULL if no query genes are found with aberration data from TCGA
  if (nrow(tcga_gene_stats) == 0) {
    lgr::lgr$info( paste0("NOTE: NO genes in query set with TCGA aberration data"))
    return(NULL)
  }

  tcga_gene_stats <- tcga_gene_stats |>
    dplyr::filter(.data$variant_type == vtype)

  ## return NULL if no query genes are found with copy number data from TCGA
  if (nrow(tcga_gene_stats) == 0) {
    lgr::lgr$info( paste0("NOTE: NO genes in query set with TCGA aberration data - type '", vtype, "'"))
    return(NULL)
  }

  tcga_gene_stats <- tcga_gene_stats |>
    dplyr::left_join(tcgadb[['site_code']],
                     by = "site_code", multiple = "all") |>
    dplyr::left_join(tcgadb[['diagnosis_code']],
                     by = "diagnosis_code", multiple = "all") |>
    dplyr::left_join(tcgadb[['clinical_strata_code']],
                     by = "clinical_strata_code", multiple = "all") |>
    dplyr::select(-c("site_code",
                     "diagnosis_code",
                     "clinical_strata_code")) |>
    dplyr::filter(.data$clinical_strata == cstrata) |>
    dplyr::filter(.data$primary_site != "Other/Unknown")

  ## return NULL if less than 3 query genes are found with copy number data from TCGA
  num_genes <- length(unique(tcga_gene_stats$symbol))
  if (num_genes < 3) {
    lgr::lgr$info( paste0("NOTE: Limited number of genes (< 3) in query set with TCGA aberration data - type '", vtype, "'"))
    return(NULL)
  }

  gene_candidates_init <- data.frame()
  tcga_ttypes <- sort(unique(tcga_gene_stats$primary_site))
  for (i in 1:length(sort(unique(tcga_gene_stats$symbol)))) {
    init <- data.frame('primary_site' <- tcga_ttypes,
                       'primary_diagnosis' = NA,
                       'symbol' = sort(unique(tcga_gene_stats$symbol))[i],
                       'clinical_strata' = cstrata,
                       'percent_mutated' = 0,
                       'percentile' = 0,
                       'variant_type' = vtype,
                       'decile' = 0,
                       stringsAsFactors = F)
    colnames(init) <- c('primary_site',
                        'primary_diagnosis',
                        'symbol',
                        'clinical_strata',
                        'percent_mutated',
                        'percentile',
                        'variant_type',
                        'decile')
    gene_candidates_init <- rbind(gene_candidates_init, init)
    i <- i + 1
  }
  gene_candidates_init <- gene_candidates_init |>
    dplyr::filter(.data$primary_site != 'Other/Unknown' &
                    .data$primary_site != "Pancancer")

  gene_aberrations <- tcga_gene_stats |>
    dplyr::filter(.data$primary_site != "Pancancer" &
                    .data$primary_site != "Other/Unknown")


  site_stats_zero <- tcga_gene_stats |>
    dplyr::select(c("primary_site",
                  "tot_samples")) |>
    dplyr::distinct() |>
    dplyr::mutate(samples_mutated = 0)

  pancan_order <- tcga_gene_stats |>
    dplyr::filter(.data$primary_site == "Pancancer") |>
    dplyr::mutate(pancancer_percent_mutated = .data$percent_mutated) |>
    dplyr::mutate(pancancer_percentile = .data$percentile) |>
    dplyr::select(c("symbol",
                  "pancancer_percent_mutated",
                  "pancancer_percentile"))

  zero_frequency_genes <-
    dplyr::anti_join(gene_candidates_init, gene_aberrations,
                     by = c("symbol", "primary_site", "variant_type")) |>
    dplyr::left_join(site_stats_zero,
                     by = c("primary_site"), multiple = "all")

  gene_aberrations <-
    dplyr::left_join(
      dplyr::bind_rows(gene_aberrations, zero_frequency_genes),
                     pancan_order, by = c("symbol"), multiple = "all") |>
    dplyr::mutate(
      pancancer_percent_mutated =
        dplyr::if_else(is.na(.data$pancancer_percent_mutated),
                       as.numeric(0),
                       as.numeric(.data$pancancer_percent_mutated)))


  gene_aberrations <- gene_aberrations |>
    dplyr::arrange(.data$pancancer_percent_mutated,
                   .data$primary_site)

  top_mutated <- gene_aberrations |>
    dplyr::arrange(dplyr::desc(.data$pancancer_percent_mutated)) |>
    dplyr::select(c("symbol")) |>
    dplyr::distinct() |>
    utils::head(75)

  gene_aberrations_top <- gene_aberrations |>
    dplyr::inner_join(top_mutated, by = "symbol", multiple = "all") |>
    dplyr::mutate(symbol = factor(.data$symbol, unique(.data$symbol)))


  gene_aberration_top_mat <- as.data.frame(
    gene_aberrations_top |>
    dplyr::select(c("symbol", "primary_site",
                  "percent_mutated")) |>
    tidyr::pivot_wider(names_from = "primary_site",
                       values_from = "percent_mutated")
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
                                  tcgadb = NULL,
                                  vtype = "snv_indel") {

  lgr::lgr$info( paste0("TCGA: collecting gene aberration data table, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(tcgadb))
  stopifnot(
    identical(names(tcgadb),
              c("coexpression",
                "aberration",
                "recurrent_variants",
                "median_ttype_expression",
                "pfam",
                "maf_codes",
                "maf",
                "site_code",
                "diagnosis_code",
                "clinical_strata_code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
  validate_db_df(tcgadb$site_code, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$diagnosis_code, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$clinical_strata_code, dbtype = "tcga_clinical_strata_code")

  stopifnot(vtype %in% c("snv_indel", "cna_homdel", "cna_ampl"))
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if (qsource == "entrezgene") {
    stopifnot(is.integer(qgenes))
    query_genes_df <- data.frame(entrezgene = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(
      genedb, query_genes_df,
      by = "entrezgene") |>
      dplyr::distinct()
  } else {
    stopifnot(is.character(qgenes))
    query_genes_df <- dplyr::inner_join(
      genedb, query_genes_df,
      by = "symbol") |>
      dplyr::distinct()
  }

  aberration_data <- tcgadb[['aberration']] |>
    dplyr::inner_join(
      dplyr::select(query_genes_df, c("symbol",
                                    "entrezgene")),
                      by=c("symbol"), multiple = "all") |>
    dplyr::left_join(tcgadb[['site_code']],
                     by = "site_code", multiple = "all") |>
    dplyr::left_join(tcgadb[['diagnosis_code']],
                     by = "diagnosis_code", multiple = "all") |>
    dplyr::left_join(tcgadb[['clinical_strata_code']],
                     by = "clinical_strata_code", multiple = "all") |>
    dplyr::select(-c("site_code",
                     "diagnosis_code",
                     "clinical_strata_code")) |>
    dplyr::filter(.data$clinical_strata == "site_diagnosis" &
                    .data$variant_type == vtype &
                    .data$primary_site != "Pancancer") |>
    dplyr::select(c("symbol",
                  "entrezgene",
                  "primary_site",
                  "primary_diagnosis",
                  "variant_type",
                  "samples_mutated",
                  "tot_samples",
                  "percent_mutated",
                  "percentile")) |>
    dplyr::rename(cohort_size = "tot_samples") |>
    dplyr::filter(!stringr::str_detect(.data$primary_diagnosis,"^Other")) |>
    dplyr::distinct() |>
    dplyr::mutate(gene = paste0(
      "<a href ='http://www.ncbi.nlm.nih.gov/gene/",
      .data$entrezgene,"' target='_blank'>",.data$symbol,"</a>")) |>
    dplyr::select(-c("entrezgene", "symbol")) |>
    dplyr::select(c("gene",
                  "variant_type",
                  "primary_site",
                  "primary_diagnosis",
                  "percent_mutated",
                  "samples_mutated",
                  "cohort_size",
                  "percentile"),
                  dplyr::everything())

  return(aberration_data)
}

tcga_coexpression <- function(qgenes,
                               qsource = "symbol",
                               genedb = NULL,
                               tcgadb = NULL) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( "TCGA: collecting co-expression data (strong negative and positive correlations)")
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(tcgadb))
  stopifnot(
    identical(names(tcgadb),
              c("coexpression",
                "aberration",
                "recurrent_variants",
                "median_ttype_expression",
                "pfam",
                "maf_codes",
                "maf",
                "site_code",
                "diagnosis_code",
                "clinical_strata_code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$coexpression, dbtype = "tcga_coexpression")
  validate_db_df(tcgadb$site_code, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$diagnosis_code, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$clinical_strata_code, dbtype = "tcga_clinical_strata_code")

  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if (qsource == "entrezgene") {
    stopifnot(is.integer(qgenes))
    query_genes_df <- data.frame(entrezgene = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(
      dplyr::select(genedb, c("entrezgene","symbol")),
      query_genes_df, by = "entrezgene", multiple = "all") |>
      dplyr::distinct()
  } else {
    stopifnot(is.character(qgenes))
    query_genes_df <- dplyr::inner_join(
      dplyr::select(genedb, c("entrezgene", "symbol")),
      query_genes_df, by = "symbol", multiple = "all") |>
      dplyr::distinct()
  }

  coexp_target_1 <- tcgadb[['coexpression']] |>
    dplyr::mutate(corrtype = dplyr::if_else(
      .data$r < 0,
      "Negative",
      as.character("Positive"))) |>
    dplyr::mutate(correlation = dplyr::case_when(
      .data$r < -0.6 & .data$r >= -0.8 ~"Strong negative",
      .data$r < -0.8 ~ "Very strong negative",
      .data$r >= 0.6 & .data$r < 0.8 ~"Strong positive",
      .data$r >= 0.8 ~ "Very strong positive",
      TRUE ~ as.character(NA))) |>
    dplyr::select(c("symbol",
                  "symbol_partner",
                  "corrtype",
                  "correlation",
                  "r",
                  "p_value",
                  "tumor")) |>
    dplyr::left_join(
      query_genes_df, by = c("symbol" = "symbol"), multiple = "all") |>
    dplyr::filter(!is.na(.data$entrezgene))

  coexp_target_tcga <- tcgadb[['coexpression']] |>
    dplyr::mutate(corrtype = dplyr::if_else(
      .data$r < 0,
      "Negative",
      as.character("Positive"))) |>
    dplyr::mutate(correlation = dplyr::case_when(
      .data$r < -0.6 & .data$r >= -0.8 ~"Strong negative",
      .data$r < -0.8 ~ "Very strong negative",
      .data$r >= 0.6 & .data$r < 0.8 ~"Strong positive",
      .data$r >= 0.8 ~ "Very strong positive",
      TRUE ~ as.character(NA))) |>
    dplyr::select(c("symbol",
                  "symbol_partner",
                  "corrtype",
                  "correlation",
                  "r",
                  "p_value",
                  "tumor")) |>
    dplyr::left_join(query_genes_df,
                     by = c("symbol_partner" = "symbol"), multiple = "all") |>
    dplyr::filter(!is.na(.data$entrezgene))

  if (NROW(coexp_target_tcga) == 0) {
    return(coexp_target_tcga)
  }
  coexp_target_tcga <- coexp_target_tcga |>
    dplyr::mutate(tmp = .data$symbol) |>
    dplyr::mutate(symbol = .data$symbol_partner) |>
    dplyr::mutate(symbol_partner = .data$tmp) |>
    dplyr::select(-c("tmp")) |>
    dplyr::bind_rows(coexp_target_1) |>
    dplyr::distinct() |>
    dplyr::left_join(
      dplyr::select(genedb,
                    c("name","oncogene","cancer_driver",
                      "tumor_suppressor",
                      "symbol","SM_tractability_category")),
      by = c("symbol_partner" = "symbol"), multiple = "all") |>
    dplyr::rename(target_gene = "symbol",
                  partner_gene = "symbol_partner",
                  partner_genename = "name",
                  target_tractability = "SM_tractability_category") |>
    dplyr::mutate(r = as.numeric(round(.data$r, digits = 3))) |>
    dplyr::filter(stringr::str_detect(.data$tumor,"BRCA|LUAD|SKCM|COAD|SARC|PRAD|ESCA|MESO|UCEC|OV|CHOL|THCA|COAD|BLCA|STAD|KIRP|GBM|HNSC")) |>
    dplyr::mutate(primary_site = "Breast") |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LUAD" | .data$tumor == "LUSC","Lung", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "STAD" | .data$tumor == "ESCA","Esophagus/Stomach", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "BLCA","Bladder", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "SARC","Soft Tissue", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "HNSC","Head and Neck", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "UCEC","Endometrium", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LAML","Myeloid", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "SKCM","Skin", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "GBM","Brain", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "PAAD","Pancreas", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "COAD" | .data$tumor == "READ","Colon/Rectum", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "PRAD","Prostate", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "THCA","Thyroid", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "KIRP" | .data$tumor == "KICH","Kidney", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "OV","Ovary", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "MESO","Pleura", as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "CHOL","Biliary Tract", as.character(.data$primary_site)))

  coexp_target_tcga <- coexp_target_tcga |>
    dplyr::filter(stringr::str_detect(.data$correlation,"Very") |
                    !is.na(.data$oncogene) |
                    !is.na(.data$tumor_suppressor) |
                    !is.na(.data$cancer_driver)) |>
    dplyr::arrange(dplyr::desc(.data$r)) |>
    dplyr::select(c("target_gene",
                  "partner_gene",
                  "correlation",
                  "r",
                  "p_value",
                  "primary_site"),
                  dplyr::everything())

  coexp_target_tcga_positive <-
    dplyr::filter(coexp_target_tcga, .data$corrtype == "Positive")
  coexp_target_tcga_negative <-
    dplyr::filter(coexp_target_tcga, .data$corrtype == "Negative") |>
    dplyr::arrange(.data$r)

  coexp_target_tcga <- as.data.frame(
    dplyr::bind_rows(coexp_target_tcga_negative,
                     coexp_target_tcga_positive)
  )

  ### remove duplicates
  duplicated_recs <- coexp_target_tcga |>
    dplyr::select(
      c("target_gene",
      "partner_gene",
      "tumor")) |>
    dplyr::inner_join(
      coexp_target_tcga, by =
        c("target_gene" = "partner_gene",
          "partner_gene" = "target_gene",
          "tumor" = "tumor"), multiple = "all"
    )

  if (NROW(duplicated_recs) > 0) {
    nonduplicated_recs <- coexp_target_tcga |>
      dplyr::anti_join(
        duplicated_recs,
        by = c("target_gene","partner_gene","tumor"))

    duplicated_recs <- duplicated_recs |>
      dplyr::filter(!(.data$tumor_suppressor == F &
                        .data$oncogene == F &
                        .data$cancer_driver == F))

    remaining <- duplicated_recs |>
      dplyr::select(
        c("target_gene",
        "partner_gene",
        "tumor")) |>
      dplyr::inner_join(
        duplicated_recs, by =
          c("target_gene" = "partner_gene",
            "partner_gene" = "target_gene",
            "tumor" = "tumor"), multiple = "all") |>
      dplyr::mutate(rn = dplyr::row_number()) |>
      dplyr::filter(.data$rn %% 2 != 0)

    deduplicated_recs <-
      duplicated_recs |> dplyr::anti_join(
      remaining, by =
        c("target_gene","partner_gene","tumor")
    )

    coexp_target_tcga <-
      nonduplicated_recs |>
      dplyr::bind_rows(deduplicated_recs)

  }

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

