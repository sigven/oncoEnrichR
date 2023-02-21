gene_tissue_cell_spec_cat <-
  function(qgenes = NULL,
           q_id_type = "symbol",
           resolution = "tissue",
           genedb = NULL,
           hpa_enrichment_db_df = NULL,
           hpa_expr_db_df = NULL) {


    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    stopifnot(!is.null(genedb))
    stopifnot(!is.null(qgenes))
    stopifnot(!is.null(hpa_enrichment_db_df))
    stopifnot(!is.null(hpa_expr_db_df))
    stopifnot(is.data.frame(hpa_expr_db_df) &
                length(colnames(hpa_expr_db_df)) > 35 &
                typeof(hpa_expr_db_df) == "list")
    validate_db_df(genedb, dbtype = "genedb")
    stopifnot(resolution == "tissue" | resolution == "single_cell")
    stopifnot(q_id_type == "symbol" | q_id_type == "entrezgene")
    query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
    if (q_id_type == 'entrezgene') {
      stopifnot(is.integer(qgenes))
      query_genes_df <-
        data.frame('entrezgene' = qgenes, stringsAsFactors = F) |>
        dplyr::inner_join(genedb, by = "entrezgene",
                          multiple = "all") |>
        dplyr::distinct()

    } else {
      stopifnot(is.character(qgenes))
      query_genes_df <- query_genes_df |>
        dplyr::inner_join(genedb, by = "symbol",
                          multiple = "all") |>
        dplyr::distinct()
    }
    etype <- "tissue"
    edb <- "GTEx"
    specificity_categories <-
      c('Group enriched',
        'Low tissue specificity',
        'Not detected',
        'Mixed',
        'Tissue enhanced',
        'Tissue enriched')
    source = 'tissue'

    if (resolution != "tissue") {
      etype <- "cell type"
      edb <- "HPA"
      validate_db_df(hpa_enrichment_db_df,
                     dbtype = "enrichment_db_hpa_singlecell")
      specificity_categories <-
        c('Group enriched',
          'Low cell type specificity',
          'Not detected',
          'Mixed',
          'Cell type enhanced',
          'Cell type enriched')
      source = 'single_cell'
    } else {
      validate_db_df(hpa_enrichment_db_df,
                     dbtype = "enrichment_db_hpa_tissue")
    }
    lgr::lgr$info(
      paste0(edb, ": retrieving ", etype,
             " specificity category of target genes")
    )

    specificity_groups_target <- as.data.frame(
      hpa_enrichment_db_df |>
        dplyr::inner_join(
          dplyr::select(query_genes_df,
                        c("symbol",
                        "ensembl_gene_id",
                        "name")),
          by = "ensembl_gene_id", multiple = "all")
    )


    if (resolution == "tissue") {
      specificity_groups_target <-
        specificity_groups_target |>
        dplyr::mutate(
          genename = paste0(
            "<a href='https://gtexportal.org/home/gene/",
            .data$ensembl_gene_id,"' target='_blank'>",
            .data$name,"</a>")
        )
    } else {
      specificity_groups_target <-
        specificity_groups_target |>
        dplyr::mutate(
          genename = paste0(
            "<a href='https://www.proteinatlas.org/",
            .data$ensembl_gene_id,
            "-",
            .data$symbol,
            "/celltype' target='_blank'>",
            .data$name,"</a>")
        )
    }
    specificity_groups_target <- as.data.frame(
      specificity_groups_target |>
        dplyr::group_by(.data$category) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        dplyr::mutate(tot = sum(.data$n)) |>
        dplyr::mutate(pct = round((.data$n / .data$tot) * 100, digits = 2)) |>
        dplyr::mutate(
          group = paste0("Target set (n = ", .data$tot,")"
          )
        )
    )

    tot <- unique(specificity_groups_target$tot)

    for (n in specificity_categories) {
      df <- data.frame(
        'category' = n,
        'pct' = 0,
        'n' = 0,
        'tot' = 0,
        group = paste0("Target set (n = ",
                       formatC(tot, format="f",
                               big.mark = ",", digits=0),")"))
      if (nrow(dplyr::inner_join(
        df,
        specificity_groups_target,
        by = "category", multiple = "all")) == 0) {
        specificity_groups_target <-
          specificity_groups_target |>
          dplyr::bind_rows(df)

      }
    }

    specificity_groups_all <- as.data.frame(
      hpa_enrichment_db_df |>
        dplyr::group_by(.data$category) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        dplyr::mutate(tot = sum(.data$n)) |>
        dplyr::mutate(pct = round((.data$n / .data$tot) * 100, digits = 2)) |>
        dplyr::mutate(group =
                        paste0("All HPA proteins (n = ",
                               formatC(.data$tot, format="f",
                                       big.mark = ",",
                                       digits=0),")"
                        )
        )
    )

    category_df <- dplyr::bind_rows(
      specificity_groups_all,
      specificity_groups_target
    )

    exp_dist_df <- data.frame()
    exp_dist_df <- hpa_expr_db_df
    exp_dist_df$ensembl_gene_id <-
      rownames(exp_dist_df)

    if (resolution == "tissue") {
      exp_dist_df <- as.data.frame(
        exp_dist_df |>
          tidyr::pivot_longer(cols = !tidyr::starts_with("ensembl"),
                              names_to = "tissue",
                              values_to = "nTPM") |>
          dplyr::inner_join(
            dplyr::select(
              query_genes_df,
              c("symbol",
                "ensembl_gene_id")),
            by = "ensembl_gene_id", multiple = "all") |>
          dplyr::mutate(exp = round(log2(.data$nTPM + 1), digits = 3)) |>
          dplyr::mutate(exp_measure = "log2(nTPM + 1)")
      )
    } else {
      exp_dist_df <- as.data.frame(
        exp_dist_df |>
          tidyr::pivot_longer(cols = !tidyr::starts_with("ensembl"),
                              names_to = "cell_type",
                              values_to = "nTPM") |>
          dplyr::inner_join(
            dplyr::select(
              query_genes_df,
              c("symbol",
                "ensembl_gene_id")),
            by = "ensembl_gene_id", multiple = "all") |>
          dplyr::mutate(exp = round(log2(.data$nTPM + 1), digits = 3)) |>
          dplyr::mutate(exp_measure = "log2(nTPM  + 1)")
      )
    }

    return(list('category_df' = category_df,
                'exp_dist_df' = exp_dist_df))

  }


gene_tissue_cell_enrichment <-
  function(qgenes_entrez = NULL,
           background_entrez = NULL,
           genedb = NULL,
           hpa_enrichment_db_df = NULL,
           hpa_enrichment_db_SE = NULL,
           resolution = "tissue") {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    stopifnot(!is.null(genedb))
    stopifnot(!is.null(hpa_enrichment_db_df))
    stopifnot(!is.null(hpa_enrichment_db_SE))
    stopifnot(!is.null(qgenes_entrez) &
                is.integer(qgenes_entrez))
    stopifnot(resolution == "tissue" | resolution == "single_cell")
    validate_db_df(genedb, dbtype = "genedb")

    #if (!("GSEABase" %in% (.packages(all.available = T)))) {
    #  suppressPackageStartupMessages(library(GSEABase))
    #}

    etype <- "tissues"
    edb <- "GTEx"
    if (resolution != "tissue") {
      etype <- "cell types"
      edb <- "HPA"
      validate_db_df(hpa_enrichment_db_df,
                     dbtype = "enrichment_db_hpa_singlecell")
    } else {
      validate_db_df(hpa_enrichment_db_df,
                     dbtype = "enrichment_db_hpa_tissue")
    }
    lgr::lgr$info(
      paste0(edb, ": estimating enrichment of ", etype,
             " in target set with TissueEnrich"))

    df <- data.frame('entrezgene' = as.integer(qgenes_entrez),
                     stringsAsFactors = F) |>
      dplyr::inner_join(
        dplyr::select(genedb,
                      c("entrezgene",
                      "ensembl_gene_id",
                      "symbol",
                      "name",
                      "cancer_max_rank")),
        by = "entrezgene", multiple = "all")

    if (resolution == "tissue") {
      df <- df |>
        dplyr::mutate(
          genename = paste0(
            "<a href='https://gtexportal.org/home/gene/",
            .data$ensembl_gene_id,"' target='_blank'>",
            .data$name,"</a>")
        )
    } else {
      df <- df |>
        dplyr::mutate(
          genename = paste0(
            "<a href='https://www.proteinatlas.org/",
            .data$ensembl_gene_id,
            "-",
            .data$symbol,
            "/celltype' target='_blank'>",
            .data$name,"</a>")
        )
    }

    bg <- hpa_enrichment_db_df

    q <- bg |>
      dplyr::select("ensembl_gene_id") |>
      dplyr::inner_join(
        df, by = "ensembl_gene_id", multiple = "all") |>
      dplyr::distinct()
    query_ensembl <- q$ensembl_gene_id

    specificities_per_gene <- bg |>
      dplyr::filter(!is.na(.data$ensembl_gene_id)) |>
      dplyr::inner_join(
        df, by = "ensembl_gene_id", multiple = "all")

    if (nrow(specificities_per_gene) > 0) {

      if (resolution == "tissue") {
        specificities_per_gene <- specificities_per_gene |>
          dplyr::select(c("symbol",
                         "genename",
                        "category",
                        "tissue",
                        "cancer_max_rank")) |>
          dplyr::arrange(.data$category,
                         dplyr::desc(.data$cancer_max_rank))
      } else {
        specificities_per_gene <- specificities_per_gene |>
          dplyr::select(c("symbol",
                        "genename",
                        "category",
                        "cell_type",
                        "cancer_max_rank")) |>
          dplyr::arrange(.data$category,
                         dplyr::desc(.data$cancer_max_rank))
      }
    }

    background_ensembl <- bg$ensembl_gene_id
    if (!is.null(background_entrez)) {

      df <-
        data.frame('entrezgene' = as.integer(background_entrez),
                   stringsAsFactors = F) |>
        dplyr::inner_join(
          dplyr::select(
            genedb, c("ensembl_gene_id", "entrezgene")),
          by = "entrezgene", multiple = "all"
        )
      bg <- bg |>
        dplyr::inner_join(
          df, by = "ensembl_gene_id", multiple = "all") |>
        dplyr::distinct()
      background_ensembl <- bg$ensembl_gene_id

      background_ensembl <- unique(
        c(background_ensembl, query_ensembl)
      )
    }

    gene_id_TE <- GSEABase::ENSEMBLIdentifier()
    gs_query <- GSEABase::GeneSet(
      geneIds = query_ensembl,
      organism = "Homo Sapiens",
      geneIdType = gene_id_TE)

    gs_background <- GSEABase::GeneSet(
      geneIds = background_ensembl,
      organism = "Homo Sapiens",
      geneIdType = gene_id_TE)

    te_output <- NULL

    ## get pre-defined tissue-specificities (SummarisedExperiement) - created with
    ## tissueEnrich::teGeneRetrieval on expression data sets
    se <- hpa_enrichment_db_SE
    #se <- oeDB$tissuecelldb[[resolution]][['te_SE']]

    ## perform tissue enrichment analysis for query dataset (gs_query),
    ## using gs_background as the background dataset, and the annotated
    ## tissue/cell-type-specific genes in se as basis for detection
    ## of enriched tissues/cell types
    suppressMessages(
      te_output <- TissueEnrich::teEnrichmentCustom(
        inputGenes = gs_query,
        tissueSpecificGenes = se,
        backgroundGenes = gs_background)
    )

    enrichment_df <- data.frame()

    if (!is.null(te_output)) {
      if (!is.null(te_output[[1]])) {
        se_enrich_output <- te_output[[1]]
        enrichment_df <- stats::setNames(
          data.frame(SummarizedExperiment::assay(se_enrich_output),
                     row.names = SummarizedExperiment::rowData(se_enrich_output)[,1]),
          SummarizedExperiment::colData(se_enrich_output)[,1])

        enrichment_df$Tissue <- rownames(enrichment_df)
        rownames(enrichment_df) <- NULL
        enrichment_df <- enrichment_df |>
          dplyr::filter(.data$Tissue != "All") |>
          dplyr::rename(fold_change = "fold.change",
                        tissue = "Tissue",
                        tissue_specific_genes = "Tissue.Specific.Genes",
                        log10_pvalue = "Log10PValue") |>
          dplyr::select(c("tissue",
                        "tissue_specific_genes",
                        "fold_change",
                        "log10_pvalue")) |>
          dplyr::arrange(dplyr::desc(.data$fold_change))

        if (resolution == "single_cell") {
          enrichment_df <- enrichment_df |>
            dplyr::rename(cell_type = "tissue",
                          celltype_specific_genes =
                            "tissue_specific_genes")
        }

      }
    }

    return(list('per_gene' = specificities_per_gene,
                'per_type' = enrichment_df))

  }

