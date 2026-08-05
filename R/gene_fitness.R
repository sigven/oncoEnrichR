get_fitness_lof_scores <-
  function(qgenes,
           qsource = "symbol",
           cellmodeldb = NULL,
           genedb = NULL,
           max_fitness_score = -2) {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    lgr::lgr$info( paste0(
      "Cell Model Passports (DepMap): retrieval of genes ",
      "associated with loss-of-fitness in cancer cell lines"))
    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(cellmodeldb))
    stopifnot(!is.null(cellmodeldb[['fitness_data']]))
    validate_db_df(cellmodeldb[['fitness_data']][['scores']],
                   dbtype = "fitness_scores")
    validate_db_df(cellmodeldb[['fitness_data']][['models']],
                   dbtype = "fitness_models")
    validate_db_df(cellmodeldb[['fitness_data']][['genes']],
                   dbtype = "fitness_genes")
    stopifnot(!is.null(genedb))
    validate_db_df(genedb, dbtype = "genedb")

    target_genes <- data.frame(
      "symbol" = qgenes, stringsAsFactors = F)

    fitness_lof_results <- list()
    fitness_lof_results[["targets"]] <- data.frame()
    fitness_lof_results[["n_targets"]] <- 0
    fitness_lof_results[["n_essential"]] <- 0
    fitness_lof_results[["common_essential"]] <- data.frame()

    fitness_genes <- as.data.frame(
      target_genes |>
        dplyr::inner_join(
          cellmodeldb[['fitness_data']][['genes']],
          by = c("symbol"),
          relationship = "many-to-many")
    )
    if (NROW (fitness_genes) > 0){
      fitness_lof_hits <- as.data.frame(
        fitness_genes |>
          dplyr::inner_join(
            cellmodeldb[['fitness_data']][['scores']],
            by = c("gene_id_cmpassports"),
            relationship = "many-to-many") |>
          dplyr::inner_join(
            cellmodeldb[['fitness_data']][['models']],
            by = c("model_id"), relationship = "many-to-many") |>
          dplyr::left_join(
            dplyr::select(
              genedb, c("entrezgene", "name")),
            by = c("entrezgene"),
          ) |>
          dplyr::rename(genename = "name") |>
          dplyr::arrange(.data$scaled_BF)
      )
    } else {
      fitness_lof_hits <- data.frame()
    }


    if (nrow(fitness_lof_hits) > 0) {

      fitness_lof_results[["targets"]] <- as.data.frame(
        fitness_lof_hits |>
          dplyr::mutate(
            model_link_ps = paste0(
              "<a href='https://cellmodelpassports.sanger.ac.uk/passports/",
              .data$model_id,"' target='_blank'>",
              stringr::str_replace_all(
                .data$model_name,"\\.","-"),"</a>")) |>
          dplyr::mutate(
            symbol_link_ps = paste0(
              "<a href='https://www.ncbi.nlm.nih.gov/gene/",
              .data$entrezgene,"' target='_blank'>",
              .data$symbol,"</a>")) |>
          dplyr::select(c("symbol",
                          "symbol_link_ps",
                          "genename",
                          "essential_gene",
                          "model_name",
                          "tissue",
                          "model_link_ps",
                          "cancer_type",
                          "sample_site",
                          "tissue_status",
                          "scaled_BF")) |>
          dplyr::filter(.data$scaled_BF <= max_fitness_score)
      )

      gene_pr_tissue_stats <- as.data.frame(
        fitness_lof_results[["targets"]] |>
          dplyr::group_by(
            .data$symbol, .data$tissue) |>
          dplyr::summarise(
            n_models_dependent_tissue = dplyr::n(),
            .groups = "drop")
      )

      gene_stats <- as.data.frame(
        gene_pr_tissue_stats |>
          dplyr::mutate(tissue_number = paste0(
            .data$tissue,":",.data$n_models_dependent_tissue)) |>
          dplyr::arrange(
            .data$symbol,
            dplyr::desc(.data$n_models_dependent_tissue)) |>
          dplyr::group_by(.data$symbol) |>
          dplyr::summarise(
            model_dependencies = paste0(
              .data$tissue_number,
              collapse = ", "),
            n_models_dependent_total = sum(
              .data$n_models_dependent_tissue),
            .groups = "drop")
      )

      fitness_lof_results[['targets']] <-
        fitness_lof_results[['targets']] |>
        dplyr::left_join(
          gene_pr_tissue_stats,
          by = c("symbol","tissue"),
          relationship = "many-to-many") |>
        dplyr::left_join(
          gene_stats,
          by = "symbol",
          relationship = "many-to-many") |>
        dplyr::select(
          c("symbol",
            "symbol_link_ps",
            "genename",
            "essential_gene",
            "model_name",
            "scaled_BF",
            "tissue",
            "model_link_ps",
            "cancer_type",
            "sample_site",
            "tissue_status",
            "model_dependencies",
            "n_models_dependent_tissue",
            "n_models_dependent_total")) |>
        dplyr::distinct()

      fitness_lof_results[["n_targets"]] <-
        NROW(gene_stats)
      fitness_lof_results[["n_essential"]] <-
        fitness_lof_results[['targets']] |>
        dplyr::filter(.data$essential_gene == TRUE) |>
        dplyr::select(c("symbol")) |>
        dplyr::distinct() |>
        nrow()

      if(fitness_lof_results[['n_essential']] > 0){
        lgr::lgr$info( paste0("Detected n = ",
                              fitness_lof_results[['n_essential']],
                              " common essential genes associated with loss-of-fitness"))
        fitness_lof_results[['common_essential']] <-
          fitness_lof_results[['targets']] |>
          dplyr::filter(.data$essential_gene == TRUE) |>
          dplyr::select(c("symbol_link_ps",
                          "genename",
                          "model_dependencies",
                          "n_models_dependent_total")) |>
          dplyr::arrange(dplyr::desc(.data$n_models_dependent_total)) |>
          dplyr::distinct()
      } else {
        lgr::lgr$info( "No essential genes associated with loss-of-fitness detected")
      }
      fitness_lof_results[['targets']]$genename <-
        NULL
      fitness_lof_results[['targets']]$model_dependencies <-
        NULL
    }

    return(fitness_lof_results)
  }


get_target_priority_scores <-
  function(qgenes,
           qsource = "symbol",
           cellmodeldb = NULL) {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    lgr::lgr$info( paste0("Cell Model Passports (DepMap): retrieval of ",
                      "prioritized targets from loss-of-fitness screens ",
                      "in cancer cell lines"))
    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(cellmodeldb))
    stopifnot(!is.null(cellmodeldb[['target_priority_scores']]))
    validate_db_df(cellmodeldb[['target_priority_scores']],
                                 dbtype = "target_priority_scores")

    target_genes <- data.frame(
      "symbol" = qgenes, stringsAsFactors = F)

    prioritized_targets <- list()
    prioritized_targets[["targets"]] <- data.frame()
    prioritized_targets[["n_pri_targets"]] <- 0

    targets_crispr_priority <- as.data.frame(
      cellmodeldb[['target_priority_scores']] |>
        dplyr::inner_join(
          target_genes,
          by = c("symbol"),
          relationship = "many-to-many")
      )

    prioritized_targets[["n_pri_targets"]] <-
      length(unique(targets_crispr_priority$symbol))

    if (nrow(targets_crispr_priority) > 0) {

      prioritized_targets[["targets"]] <- as.data.frame(
        targets_crispr_priority |>
        dplyr::mutate(symbol = factor(.data$symbol, levels =
                 levels(cellmodeldb[['target_priority_scores']]$symbol))) |>
        dplyr::arrange(dplyr::desc(.data$priority_score)) |>
        dplyr::select(c("symbol",
                       "tumor_type",
                       "priority_score"))

      )
    }

    return(prioritized_targets)

  }

