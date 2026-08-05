gene_tissue_celltype_specificity <-
  function(qgenes = NULL,
           genedb = NULL,
           tissuecelldb = NULL,
           resolution = "tissue") {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(tissuecelldb))
    stopifnot(resolution == "tissue" | resolution == "celltype")
    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(tissuecelldb[['specificity']],
                   dbtype = "tissuecelldb_specificity")

    ctype <- "tissue"
    url_suffix <- "tissue"
    if (resolution != "tissue") {
      ctype <- "cell type"
      url_suffix <- "celltype"
    }

    lgr::lgr$info(
      glue::glue(
        "HPA: Retrieving {ctype} expression specificity ",
        "(tau index) for target genes"))

    query_genes_df <- data.frame(
      'symbol' = qgenes, stringsAsFactors = F) |>
      dplyr::inner_join(
        dplyr::select(genedb,
                      c("symbol",
                        "entrezgene",
                        "ensembl_gene_id",
                        "name",
                        "cancer_max_rank")),
        by = "symbol", relationship = "many-to-many") |>
      dplyr::distinct()

    per_gene <- tissuecelldb[['specificity']] |>
      dplyr::inner_join(
        query_genes_df, by = "ensembl_gene_id",
        relationship = "many-to-many") |>
      dplyr::mutate(
        genename = glue::glue(
          "<a href='https://www.proteinatlas.org/",
          "{.data$ensembl_gene_id}-{.data$symbol}/{url_suffix}' ",
          "target='_blank'>{.data$name}</a>")) |>
      dplyr::rename(top_context = "max_context") |>
      dplyr::select(
        c("symbol",
          "genename",
          "cancer_max_rank",
          "tau",
          "top_context",
          "top_contexts",
          "top_expressions")) |>
      dplyr::arrange(dplyr::desc(.data$tau))

    tau_distribution <- dplyr::bind_rows(
      tissuecelldb[['specificity']] |>
        dplyr::select("tau") |>
        dplyr::mutate(gene_set = "All protein-coding genes"),
      per_gene |>
        dplyr::select("tau") |>
        dplyr::mutate(gene_set = "Query set")
    )

    return(
      list(
        'per_gene' = per_gene,
        'tau_distribution' = tau_distribution
      )
    )

  }


gene_tissue_celltype_enrichment <-
  function(qgenes = NULL,
           genedb = NULL,
           tissuecelldb = NULL,
           resolution = "tissue",
           background_entrez = NULL,
           p_value_adjustment_method = "BH") {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(tissuecelldb))
    stopifnot(resolution == "tissue" | resolution == "celltype")
    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(tissuecelldb[['enrichment_membership']],
                   dbtype = "tissuecelldb_enrichment_membership")

    ctype <- "tissues"
    url_suffix <- "tissue"
    if (resolution != "tissue") {
      ctype <- "cell types"
      url_suffix <- "celltype"
    }

    lgr::lgr$info(
      glue::glue(
        "HPA: Estimating enrichment of {ctype} in target set ",
        "(hypergeometric test of highly expressed genes)"))

    query_genes_df <- data.frame(
      'symbol' = qgenes, stringsAsFactors = F) |>
      dplyr::inner_join(
        dplyr::select(genedb,
                      c("symbol",
                        "entrezgene",
                        "ensembl_gene_id",
                        "name",
                        "cancer_max_rank")),
        by = "symbol", relationship = "many-to-many") |>
      dplyr::distinct()

    membership <- tissuecelldb[['enrichment_membership']]
    background_summary <- tissuecelldb[['enrichment_background_summary']]

    ## Restrict background (and its per-context gene counts) to a
    ## custom set of genes, if provided - mirrors the 'bgset' concept
    ## used elsewhere in oncoEnrichR (e.g. GO/pathway enrichment)
    if (!is.null(background_entrez)) {
      background_ensembl <- data.frame(
        'entrezgene' = as.integer(background_entrez),
        stringsAsFactors = F) |>
        dplyr::inner_join(
          dplyr::select(genedb, c("entrezgene", "ensembl_gene_id")),
          by = "entrezgene", relationship = "many-to-many") |>
        dplyr::distinct() |>
        dplyr::pull(.data$ensembl_gene_id)
      background_ensembl <- unique(
        c(background_ensembl, query_genes_df$ensembl_gene_id))

      membership <- membership |>
        dplyr::filter(.data$ensembl_gene_id %in% background_ensembl)
      background_summary <- membership |>
        dplyr::count(.data$category, name = "n_bg_genes") |>
        dplyr::mutate(n_total_genes = length(background_ensembl))
    }

    ## Per-gene membership: which contexts each query gene is
    ## "highly expressed" in (top quantile, precomputed)
    per_gene <- membership |>
      dplyr::inner_join(
        query_genes_df, by = "ensembl_gene_id",
        relationship = "many-to-many") |>
      dplyr::mutate(
        genename = glue::glue(
          "<a href='https://www.proteinatlas.org/",
          "{.data$ensembl_gene_id}-{.data$symbol}/{url_suffix}' ",
          "target='_blank'>{.data$name}</a>")) |>
      dplyr::select(
        c("symbol", "genename", "cancer_max_rank", "category")) |>
      dplyr::arrange(.data$category, dplyr::desc(.data$cancer_max_rank))

    ## Per-context hypergeometric enrichment test - fully vectorized
    ## across all contexts at once (no per-context loop): count query
    ## genes per context, join precomputed background sizes, and run
    ## a single phyper() call over all contexts
    n_input_genes <- nrow(query_genes_df)

    genes_contributing <- per_gene |>
      dplyr::group_by(.data$category) |>
      dplyr::summarise(
        genes_contributing =
          paste(sort(unique(.data$symbol)), collapse = ", "),
        .groups = "drop")

    per_type <- membership |>
      dplyr::filter(
        .data$ensembl_gene_id %in% query_genes_df$ensembl_gene_id) |>
      dplyr::count(.data$category, name = "n_input_in_context") |>
      dplyr::right_join(background_summary, by = "category") |>
      dplyr::mutate(
        n_input_in_context = dplyr::if_else(
          is.na(.data$n_input_in_context), 0L, .data$n_input_in_context),
        n_input_genes = n_input_genes) |>
      dplyr::mutate(
        p_value = stats::phyper(
          q = .data$n_input_in_context - 1,
          m = .data$n_bg_genes,
          n = .data$n_total_genes - .data$n_bg_genes,
          k = .data$n_input_genes,
          lower.tail = FALSE),
        fold_enrichment = dplyr::if_else(
          .data$n_input_in_context == 0, 0,
          .data$n_input_in_context /
            (.data$n_input_genes * (.data$n_bg_genes / .data$n_total_genes)))
      ) |>
      dplyr::mutate(
        p_adjusted = stats::p.adjust(
          .data$p_value, method = p_value_adjustment_method)
      ) |>
      dplyr::left_join(genes_contributing, by = "category") |>
      dplyr::mutate(
        genes_contributing = dplyr::if_else(
          is.na(.data$genes_contributing), "", .data$genes_contributing)) |>
      dplyr::rename(
        n_query_genes = "n_input_in_context",
        n_background_genes = "n_bg_genes") |>
      dplyr::select(
        c("category", "n_query_genes", "n_background_genes",
          "fold_enrichment", "p_adjusted", "p_value",
          "genes_contributing")) |>
      dplyr::arrange(.data$p_value)

    if (resolution == "tissue") {
      per_gene <- dplyr::rename(per_gene, tissue = "category")
      per_type <- dplyr::rename(per_type, tissue = "category")
    } else {
      per_gene <- dplyr::rename(per_gene, cell_type = "category")
      per_type <- dplyr::rename(per_type, cell_type = "category")
    }

    return(
      list(
        'per_gene' = per_gene,
        'per_type' = per_type
      )
    )

  }
