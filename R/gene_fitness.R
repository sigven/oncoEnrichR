get_fitness_lof_scores <- function(qgenes,
                                  qsource = "symbol",
                                  depmapdb = NULL,
                                  max_fitness_score = -2) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( paste0("DepMap/Project Score: retrieval of genes ",
                    "associated with loss-of-fitness in cancer cell lines"))
  stopifnot(!is.null(qgenes))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(depmapdb))
  stopifnot(!is.null(depmapdb[['fitness_scores']]))
  validate_db_df(depmapdb[['fitness_scores']],
                               dbtype = "fitness_scores")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  fitness_lof_results <- list()
  fitness_lof_results[["targets"]] <- data.frame()
  fitness_lof_results[["n_targets"]] <- 0

  fitness_lof_hits <- as.data.frame(
    target_genes |>
    dplyr::inner_join(
      depmapdb[['fitness_scores']],
      by = c("symbol"), relationship = "many-to-many") |>
    dplyr::arrange(.data$scaled_BF)
  )
  if (nrow(fitness_lof_hits) > 0) {

    fitness_lof_results[["targets"]] <- as.data.frame(
      fitness_lof_hits |>
        dplyr::mutate(
          model_link_ps = paste0(
            "<a href='https://score.depmap.sanger.ac.uk/model/",
            .data$model_id,"?&scoreMax=0' target='_blank'>",
            stringr::str_replace_all(
              .data$model_name,"\\.","-"),"</a>")) |>
        dplyr::mutate(
          symbol_link_ps = paste0(
            "<a href='https://score.depmap.sanger.ac.uk/gene/",
            .data$gene_id_project_score,"' target='_blank'>",
            .data$symbol,"</a>")) |>
        dplyr::select(c("symbol",
                      "symbol_link_ps",
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
        dplyr::group_by(.data$symbol, .data$tissue) |>
        dplyr::summarise(n_gene_tissue = dplyr::n(),
                         .groups = "drop")
    )

    gene_stats <- as.data.frame(
      gene_pr_tissue_stats |>
        dplyr::group_by(.data$symbol) |>
        dplyr::summarise(n_gene = sum(.data$n_gene_tissue),
                         .groups = "drop")
    )

    fitness_lof_results[['targets']] <- fitness_lof_results[['targets']] |>
      dplyr::left_join(gene_pr_tissue_stats,
                       by = c("symbol","tissue"), relationship = "many-to-many") |>
      dplyr::left_join(gene_stats, by = "symbol", relationship = "many-to-many") |>
      dplyr::select(c("symbol",
                    "symbol_link_ps",
                    "model_name",
                    "scaled_BF",
                    "tissue",
                    "model_link_ps",
                    "cancer_type",
                    "sample_site",
                    "tissue_status",
                    "n_gene_tissue",
                    "n_gene"))

    fitness_lof_results[["n_targets"]] <- nrow(gene_stats)
  }

  return(fitness_lof_results)
}

get_target_priority_scores <-
  function(qgenes,
           qsource = "symbol",
           depmapdb = NULL) {

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

    lgr::lgr$info( paste0("DepMap/Project Score: retrieval of ",
                      "prioritized targets from loss-of-fitness screens ",
                      "in cancer cell lines"))
    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(depmapdb))
    stopifnot(!is.null(depmapdb[['target_priority_scores']]))
    validate_db_df(depmapdb[['target_priority_scores']],
                                 dbtype = "target_priority_scores")

    target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

    prioritized_targets <- list()
    prioritized_targets[["targets"]] <- data.frame()
    prioritized_targets[["n_pri_targets"]] <- 0

    targets_crispr_priority <- as.data.frame(
      depmapdb[['target_priority_scores']] |>
        dplyr::inner_join(target_genes, by = c("symbol"), relationship = "many-to-many")
      )

    prioritized_targets[["n_pri_targets"]] <-
      length(unique(targets_crispr_priority$symbol))

    if (nrow(targets_crispr_priority) > 0) {

      prioritized_targets[["targets"]] <- as.data.frame(
        targets_crispr_priority |>
        dplyr::mutate(symbol = factor(.data$symbol, levels =
                 levels(depmapdb[['target_priority_scores']]$symbol))) |>
        dplyr::arrange(dplyr::desc(.data$priority_score)) |>
        dplyr::select(c("symbol",
                       "tumor_type",
                       "priority_score"))

      )
    }

    return(prioritized_targets)

  }

