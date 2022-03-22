get_fitness_lof_scores <- function(qgenes,
                                  qsource = "symbol",
                                  projectscoredb = NULL,
                                  logger = NULL) {

  stopifnot(!is.null(logger))
  log4r_info(logger, paste0("Project Score: retrieval of genes ",
                    "associated with loss-of-fitness in cancer cell lines"))
  stopifnot(!is.null(qgenes))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(projectscoredb))
  stopifnot(!is.null(projectscoredb[['fitness_scores']]))
  validate_db_df(projectscoredb[['fitness_scores']],
                               dbtype = "fitness_scores")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  fitness_lof_results <- list()
  fitness_lof_results[["targets"]] <- data.frame()
  fitness_lof_results[["n_targets"]] <- 0

  fitness_lof_hits <- as.data.frame(
    target_genes %>%
    dplyr::inner_join(projectscoredb[['fitness_scores']], by = c("symbol"))
  )
  if (nrow(fitness_lof_hits) > 0) {

    fitness_lof_hits <- as.data.frame(
      fitness_lof_hits %>%
        dplyr::mutate(
          model_link_cmp = paste0(
            "<a href='https://cellmodelpassports.sanger.ac.uk/passports/",
            .data$model_id,"' target='_blank'>",
            stringr::str_replace_all(.data$model_name,"\\.","-"),"</a>")) %>%
        dplyr::mutate(
          symbol_link_ps = paste0(
            "<a href='https://score.depmap.sanger.ac.uk/gene/",
            .data$gene_id_project_score,"' target='_blank'>", .data$symbol,"</a>")) %>%
      dplyr::select(.data$symbol,
                    .data$symbol_link_ps,
                    .data$model_name,
                    .data$tissue,
                    .data$model_link_cmp) %>%
      dplyr::group_by(.data$symbol,
                      .data$symbol_link_ps,
                      .data$tissue) %>%
      dplyr::summarise(n_gene_tissue = dplyr::n(),
                       cell_lines = paste(.data$model_name, collapse = ", "),
                       cmp_link = paste(.data$model_link_cmp, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cell_lines = stringr::str_replace_all(.data$cell_lines, "\\.", "-"))
    )

    total <- as.data.frame(
      fitness_lof_hits %>%
        dplyr::group_by(.data$symbol) %>%
        dplyr::summarise(n_gene = sum(.data$n_gene_tissue))
      )

    fitness_lof_results[["n_targets"]] <- nrow(total)
    fitness_lof_results[["targets"]] <- fitness_lof_hits %>%
      dplyr::left_join(total, by = c("symbol"))
  }

  return(fitness_lof_results)
}

get_target_priority_scores <-
  function(qgenes,
           qsource = "symbol",
           projectscoredb = NULL,
           logger = NULL){

    stopifnot(!is.null(logger))
    log4r_info(logger, paste0("Project Score: retrieval of ",
                      "prioritized targets from loss-of-fitness screens ",
                      "in cancer cell lines"))
    stopifnot(!is.null(qgenes))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(projectscoredb))
    stopifnot(!is.null(projectscoredb[['target_priority_scores']]))
    validate_db_df(projectscoredb[['target_priority_scores']],
                                 dbtype = "target_priority_scores")

    target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

    prioritized_targets <- list()
    prioritized_targets[["targets"]] <- data.frame()
    prioritized_targets[["n_pri_targets"]] <- 0

    targets_crispr_priority <- as.data.frame(
      projectscoredb[['target_priority_scores']] %>%
        dplyr::inner_join(target_genes, by = c("symbol"))
      )

    prioritized_targets[["n_pri_targets"]] <-
      length(unique(targets_crispr_priority$symbol))

    if(nrow(targets_crispr_priority) > 0){

      prioritized_targets[["targets"]] <- as.data.frame(
        targets_crispr_priority %>%
        dplyr::mutate(symbol = factor(.data$symbol, levels =
                 levels(projectscoredb[['target_priority_scores']]$symbol))) %>%
        dplyr::arrange(dplyr::desc(.data$priority_score)) %>%
        dplyr::select(.data$symbol,
                      .data$tumor_type,
                      .data$priority_score) # %>%

      )
    }

    return(prioritized_targets)

  }

