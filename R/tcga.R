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
                  "coexpression_summary",
                  "aberration",
                  "median_ttype_expression",
                  "maf",
                  "code")
      ))
    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
    validate_db_df(tcgadb$code$site, dbtype = "tcga_site_code")
    validate_db_df(tcgadb$code$diagnosis, dbtype = "tcga_diagnosis_code")
    validate_db_df(tcgadb$code$clinical_strata, dbtype = "tcga_clinical_strata_code")
    stopifnot(site %in% unique(tcgadb$code$site$primary_site))
    stopifnot(qsource == "symbol" | qsource == "entrezgene")
    stopifnot(cstrata == "site" | cstrata == "site_diagnosis")
    query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
    if (qsource == 'entrezgene') {
      stopifnot(is.integer(qgenes))
      query_genes_df <- data.frame('entrezgene' = qgenes, stringsAsFactors = F)
      query_genes_df <- dplyr::inner_join(
        genedb, query_genes_df, by = "entrezgene", relationship = "many-to-many") |>
        dplyr::distinct()
    } else {
      stopifnot(is.character(qgenes))
      query_genes_df <- dplyr::inner_join(
        genedb, query_genes_df, by = "symbol", relationship = "many-to-many") |>
        dplyr::distinct()
    }

    top_mutated_genes <- tcgadb[['aberration']] |>
      dplyr::inner_join(
        dplyr::select(query_genes_df, "symbol"),
        by = c("symbol"), relationship = "many-to-many") |>
      dplyr::left_join(tcgadb[['code']][['site']],
                       by = "site_code",
                       relationship = "many-to-many") |>
      dplyr::left_join(tcgadb[['code']][['diagnosis']],
                       by = "diagnosis_code", relationship = "many-to-many") |>
      dplyr::left_join(tcgadb[['code']][['clinical_strata']],
                       by = "clinical_strata_code", relationship = "many-to-many") |>
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
      paste0("Choosing genes for oncoplot - highest SNV/InDel",
             "frequency in TCGA cohort - ", site))

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
                "coexpression_summary",
                "aberration",
                "median_ttype_expression",
                "maf",
                "code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
  validate_db_df(tcgadb$code$site, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$code$diagnosis, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$code$clinical_strata,
                 dbtype = "tcga_clinical_strata_code")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(cstrata == "site")
  stopifnot(vtype == "cna_ampl" | vtype == "cna_homdel")
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if (qsource == 'entrezgene') {
    stopifnot(is.integer(qgenes))
    query_genes_df <-
      data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "entrezgene", relationship = "many-to-many") |>
      dplyr::distinct()
  } else {
    stopifnot(is.character(qgenes))
    query_genes_df <-
      dplyr::inner_join(genedb, query_genes_df, by = "symbol", relationship = "many-to-many") |>
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
      by = c("symbol"), relationship = "many-to-many"
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
    dplyr::left_join(tcgadb[['code']][['site']],
                     by = "site_code",
                     relationship = "many-to-many") |>
    dplyr::left_join(tcgadb[['code']][['diagnosis']],
                     by = "diagnosis_code",
                     relationship = "many-to-many") |>
    dplyr::left_join(tcgadb[['code']][['clinical_strata']],
                     by = "clinical_strata_code",
                     relationship = "many-to-many") |>
    dplyr::select(-c("site_code",
                     "diagnosis_code",
                     "clinical_strata_code")) |>
    dplyr::filter(.data$clinical_strata == cstrata) |>
    dplyr::filter(.data$primary_site != "Other/Unknown")

  ## return NULL if less than 3 query genes are found with copy number data from TCGA
  num_genes <- length(unique(tcga_gene_stats$symbol))
  if (num_genes < 3) {
    lgr::lgr$info( paste0(
      "NOTE: Limited number of genes (< 3) in query set with ",
      "recurrent aberration data in TCGA - '", vtype, "'"))
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
                     by = c("primary_site"), relationship = "many-to-many")

  gene_aberrations <-
    dplyr::left_join(
      dplyr::bind_rows(gene_aberrations, zero_frequency_genes),
                     pancan_order, by = c("symbol"), relationship = "many-to-many") |>
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
    utils::head(70)

  gene_aberrations_top <- gene_aberrations |>
    dplyr::inner_join(top_mutated, by = "symbol", relationship = "many-to-many") |>
    dplyr::mutate(symbol = factor(.data$symbol, unique(.data$symbol)))


  gene_aberration_top_mat <- as.data.frame(
    gene_aberrations_top |>
    dplyr::select(c("symbol", "primary_site",
                  "percent_mutated")) |>
    dplyr::distinct() |>
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
                "coexpression_summary",
                "aberration",
                "median_ttype_expression",
                "maf",
                "code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$aberration, dbtype = "tcga_aberration")
  validate_db_df(tcgadb$code$site, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$code$diagnosis, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$code$clinical_strata, dbtype = "tcga_clinical_strata_code")

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
      dplyr::select(
        query_genes_df,
        c("symbol", "entrezgene")),
      by=c("symbol"),
      relationship = "many-to-many") |>
    dplyr::left_join(
      tcgadb[['code']][['site']],
      by = "site_code",
      relationship = "many-to-many") |>
    dplyr::left_join(
      tcgadb[['code']][['diagnosis']],
      by = "diagnosis_code",
      relationship = "many-to-many") |>
    dplyr::left_join(
      tcgadb[['code']][['clinical_strata']],
      by = "clinical_strata_code",
      relationship = "many-to-many") |>
    dplyr::select(
      -c("site_code",
         "diagnosis_code",
         "clinical_strata_code")) |>
    dplyr::filter(
      .data$clinical_strata == "site_diagnosis" &
        .data$variant_type == vtype &
        .data$primary_site != "Pancancer") |>
    dplyr::select(
      c("symbol",
        "entrezgene",
        "primary_site",
        "primary_diagnosis",
        "variant_type",
        "samples_mutated",
        "tot_samples",
        "percent_mutated",
        "percentile")) |>
    dplyr::rename(cohort_size = "tot_samples") |>
    dplyr::filter(!stringr::str_detect(
      .data$primary_diagnosis,"^Other")) |>
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
                "coexpression_summary",
                "aberration",
                "median_ttype_expression",
                "maf",
                "code")
    ))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tcgadb$coexpression, dbtype = "tcga_coexpression")
  validate_db_df(tcgadb$code$site, dbtype = "tcga_site_code")
  validate_db_df(tcgadb$code$diagnosis, dbtype = "tcga_diagnosis_code")
  validate_db_df(tcgadb$code$clinical_strata, dbtype = "tcga_clinical_strata_code")

  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  query_genes_df <-
    data.frame('symbol' = qgenes, stringsAsFactors = F)
  if (qsource == "entrezgene") {
    stopifnot(is.integer(qgenes))
    query_genes_df <- data.frame(entrezgene = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(
      dplyr::select(genedb, c("entrezgene","symbol")),
      query_genes_df, by = "entrezgene", relationship = "many-to-many") |>
      dplyr::distinct()
  } else {
    stopifnot(is.character(qgenes))
    query_genes_df <- dplyr::inner_join(
      dplyr::select(genedb, c("entrezgene", "symbol")),
      query_genes_df, by = "symbol", relationship = "many-to-many") |>
      dplyr::distinct()
  }

  coexp_target_1 <- data.frame()
  coexp_target_2 <- data.frame()

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
                  "adj_pvalue",
                  "tumor")) |>
    dplyr::inner_join(
      query_genes_df, by = c("symbol" = "symbol"),
      relationship = "many-to-many")

  if(NROW(coexp_target_1) > 0){
    coexp_target_1 <- coexp_target_1 |>
      dplyr::select(-c("entrezgene")) |>
      dplyr::left_join(
        dplyr::select(
          genedb,
          c("name","oncogene",
            "cancer_driver",
            "tumor_suppressor",
            "symbol","SM_tractability_category")),
        by = c("symbol_partner" = "symbol"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(
        target_gene = "symbol",
        partner_gene = "symbol_partner",
        partner_genename = "name",
        partner_tractability = "SM_tractability_category",
        partner_oncogene = "oncogene",
        partner_driver = "cancer_driver",
        partner_tumor_suppressor = "tumor_suppressor"
      ) |>
      dplyr::filter(
        .data$partner_oncogene == TRUE |
          .data$partner_driver == TRUE |
          .data$partner_tumor_suppressor == TRUE
      )

    if(NROW(coexp_target_1) > 0){
      coexp_target_1 <- coexp_target_1 |>
        dplyr::left_join(
          dplyr::select(
            tcgadb$coexpression_summary,
            c("symbol","expr","tumor")),
          by = c("partner_gene" = "symbol",
                 "tumor" = "tumor")) |>
        dplyr::rename(
          partner_expression = "expr"
        ) |>
        dplyr::left_join(
          dplyr::select(
            tcgadb$coexpression_summary,
            c("symbol","expr","tumor")),
          by = c("target_gene" = "symbol",
                 "tumor" = "tumor")) |>
        dplyr::rename(
          target_expression = "expr"
        )
    }
  }else{
    coexp_target_1 <- data.frame()
  }

  coexp_target_2 <- tcgadb[['coexpression']] |>
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
                    "adj_pvalue",
                    "tumor")) |>
    dplyr::inner_join(
      query_genes_df, by = c("symbol_partner" = "symbol"),
      relationship = "many-to-many")

  if(NROW(coexp_target_2) > 0){
    coexp_target_2 <- coexp_target_2 |>
      dplyr::select(-c("entrezgene")) |>
      dplyr::left_join(
        dplyr::select(
          genedb,
          c("name","oncogene",
            "cancer_driver",
            "tumor_suppressor",
            "symbol","SM_tractability_category")),
        by = c("symbol" = "symbol"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(
        target_gene = "symbol_partner",
        partner_gene = "symbol",
        partner_genename = "name",
        partner_tractability = "SM_tractability_category",
        partner_oncogene = "oncogene",
        partner_driver = "cancer_driver",
        partner_tumor_suppressor = "tumor_suppressor"
      ) |>
      dplyr::filter(
        .data$partner_oncogene == TRUE |
          .data$partner_driver == TRUE |
          .data$partner_tumor_suppressor == TRUE
      )

    if(NROW(coexp_target_2) > 0){
      coexp_target_2 <- coexp_target_2 |>
        dplyr::left_join(
          dplyr::select(
            tcgadb$coexpression_summary,
            c("symbol","expr","tumor")),
          by = c("partner_gene" = "symbol",
                 "tumor" = "tumor")) |>
        dplyr::rename(
          partner_expression = "expr"
        ) |>
        dplyr::left_join(
          dplyr::select(
            tcgadb$coexpression_summary,
            c("symbol","expr","tumor")),
          by = c("target_gene" = "symbol",
                 "tumor" = "tumor")) |>
        dplyr::rename(
          target_expression = "expr"
        )
    }
  }else{
    coexp_target_2 <- data.frame()
  }

  coexp_target_tcga <- dplyr::bind_rows(
    coexp_target_1,
    coexp_target_2) |>
    dplyr::distinct() |>
    dplyr::mutate(r = as.numeric(round(.data$r, digits = 3))) |>
    dplyr::filter(stringr::str_detect(
      .data$tumor,
      paste0(
        "BRCA|LUAD|LUSC|SKCM|COAD|SARC|",
        "PRAD|PAAD|ESCA|MESO|KIRC|UCEC|",
        "OV|CHOL|THCA|COAD|BLCA|STAD|",
        "KIRP|GBM|READ|DLBC|CESC|HNSC|",
        "KICH|TGCT|UVM|LGG|ACC|LAML"))) |>
    dplyr::mutate(primary_site = "Breast") |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LUAD" |
        .data$tumor == "LUSC","Lung",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "ACC","Adrenal Gland",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "TGCT","Testis",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "STAD" |
        .data$tumor == "ESCA","Esophagus/Stomach",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "BLCA","Bladder",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "SARC","Soft Tissue",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "HNSC","Head and Neck",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "UCEC","Endometrium",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LAML","Myeloid",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "SKCM","Skin",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "GBM" | .data$tumor == "LGG","Brain",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "CESC","Cervix",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "COAD" |
        .data$tumor == "READ","Colon/Rectum",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "PRAD","Prostate",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "THCA","Thyroid",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "PAAD","Pancreas",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "KIRP" |
        .data$tumor == "KIRC" |
        .data$tumor == "KICH","Kidney",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LIHC","Liver",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "OV","Ovary",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "MESO","Pleura",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "CHOL","Biliary Tract",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "UVM","Eye",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "DLBC","Lymphoid",
      as.character(.data$primary_site))) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      .data$tumor == "LAML","Myeloid",
      as.character(.data$primary_site))) |>
    dplyr::arrange(dplyr::desc(.data$r)) |>
    dplyr::select(
      c("target_gene",
        "partner_gene",
        "correlation",
        "r",
        "adj_pvalue",
        "primary_site",
        "tumor",
        "target_expression",
        "partner_expression",
        "partner_genename"),
      dplyr::everything())

  coexp_target_tcga_positive <-
    dplyr::filter(coexp_target_tcga, .data$corrtype == "Positive") |>
    dplyr::arrange(dplyr::desc(.data$r))
  coexp_target_tcga_negative <-
    dplyr::filter(coexp_target_tcga, .data$corrtype == "Negative") |>
    dplyr::arrange(.data$r)

  coexp_target_tcga <- as.data.frame(
    dplyr::bind_rows(coexp_target_tcga_negative,
                     coexp_target_tcga_positive) |>
      utils::head(200000)
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
          "tumor" = "tumor"), relationship = "many-to-many"
    )

  if (NROW(duplicated_recs) > 0) {
    nonduplicated_recs <- coexp_target_tcga |>
      dplyr::anti_join(
        duplicated_recs,
        by = c("target_gene","partner_gene","tumor"))

    duplicated_recs <- duplicated_recs |>
      dplyr::filter(!(.data$partner_tumor_suppressor == F &
                        .data$partner_oncogene == F &
                        .data$partner_driver == F))

    remaining <- duplicated_recs |>
      dplyr::select(
        c("target_gene",
        "partner_gene",
        "tumor")) |>
      dplyr::inner_join(
        duplicated_recs, by =
          c("target_gene" = "partner_gene",
            "partner_gene" = "target_gene",
            "tumor" = "tumor"), relationship = "many-to-many") |>
      dplyr::mutate(rn = dplyr::row_number()) |>
      dplyr::filter(.data$rn %% 2 != 0)

    deduplicated_recs <-
      duplicated_recs |> dplyr::anti_join(
      remaining, by =
        c("target_gene","partner_gene","tumor")
    )

    coexp_target_tcga <-
      nonduplicated_recs |>
      dplyr::bind_rows(deduplicated_recs) |>
      dplyr::arrange(.data$r) |>
      dplyr::distinct()

  }

  coexp_target_final <- coexp_target_tcga

  if(NROW(coexp_target_tcga) > 0){
    coexp_target_final <- data.frame()

    for(n in c('Positive','Negative')){
      coexp_data_raw <- coexp_target_tcga |>
        dplyr::filter(.data$corrtype == n)

      if(NROW(coexp_data_raw) > 0 &
         "primary_site" %in% colnames(coexp_data_raw)){
        for(ps in sort(unique(coexp_data_raw$primary_site))){
          coexp_data_ps <- coexp_data_raw |>
            dplyr::filter(.data$primary_site == ps)

          if(NROW(coexp_data_ps) == 0){
            next
          }
          coexp_data_ps <- coexp_data_ps |>
            dplyr::arrange(dplyr::desc(abs(.data$r))) |>
            utils::head(1000)

          coexp_target_final <-
            dplyr::bind_rows(
              coexp_target_final,
              coexp_data_ps)
        }
      }
    }
  }


  return(coexp_target_final)

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

