get_crispr_lof_scores <- function(qgenes,
                                  qsource = "symbol",
                                  projectscoredb = NULL) {

  oncoEnrichR:::log4r_info(paste0("Project Score (CRISPR/Cas9 screen): retrieval of genes ",
                    "associated with loss-of-fitness in cancer cell lines"))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(projectscoredb))
  stopifnot(!is.null(projectscoredb[['fitness_scores']]))
  oncoEnrichR:::validate_db_df(projectscoredb[['fitness_scores']],
                               dbtype = "fitness_scores")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  crispr_lof_results <- list()
  crispr_lof_results[["targets"]] <- data.frame()
  crispr_lof_results[["n_targets"]] <- 0

  crispr_lof_hits <- as.data.frame(
    target_genes %>%
    dplyr::inner_join(projectscoredb[['fitness_scores']], by = c("symbol"))
  )
  if (nrow(crispr_lof_hits) > 0) {

    crispr_lof_hits <- as.data.frame(
      crispr_lof_hits %>%
        dplyr::mutate(
          model_link_cmp = paste0(
            "<a href='https://cellmodelpassports.sanger.ac.uk/passports/",
            model_id,"' target='_blank'>",
            stringr::str_replace_all(model_name,"\\.","-"),"</a>")) %>%
        dplyr::mutate(
          symbol_link_ps = paste0(
            "<a href='https://score.depmap.sanger.ac.uk/gene/",
            gene_id_project_score,"' target='_blank'>",symbol,"</a>")) %>%
      dplyr::select(symbol, symbol_link_ps, model_name,
                    tissue, model_link_cmp) %>%
      dplyr::group_by(symbol, symbol_link_ps, tissue) %>%
      dplyr::summarise(n_gene_tissue = dplyr::n(),
                       cell_lines = paste(model_name, collapse = ", "),
                       cmp_link = paste(model_link_cmp, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cell_lines = stringr::str_replace_all(cell_lines, "\\.", "-"))
    )

    total <- as.data.frame(
      crispr_lof_hits %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(n_gene = sum(n_gene_tissue))
      )

    crispr_lof_results[["n_targets"]] <- nrow(total)
    crispr_lof_results[["targets"]] <- crispr_lof_hits %>%
      dplyr::left_join(total, by = c("symbol"))
  }

  return(crispr_lof_results)
}

get_target_priority_scores <-
  function(qgenes,
           qsource = "symbol",
           projectscoredb = NULL){

    oncoEnrichR:::log4r_info(paste0("Project Score (CRISPR/Cas9 screen): retrieval of ",
                      "prioritized targets from loss-of-fitness screens ",
                      "in cancer cell lines"))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(projectscoredb))
    stopifnot(!is.null(projectscoredb[['target_priority_scores']]))
    oncoEnrichR:::validate_db_df(projectscoredb[['target_priority_scores']],
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
        dplyr::mutate(symbol = factor(symbol, levels =
                 levels(projectscoredb[['target_priority_scores']]$symbol))) %>%
        dplyr::arrange(symbol) %>%
        dplyr::select(symbol, tumor_type, priority_score) # %>%

      )
    }

    return(prioritized_targets)

  }

