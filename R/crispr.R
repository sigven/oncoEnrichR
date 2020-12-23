get_crispr_lof_scores <- function(qgenes,
                                  qsource = "symbol",
                                  projectscoredb = NULL) {

  rlogging::message("Project Score (CRISPR/Cas9 screen): retrieval of genes associated with loss-of-fitness in cancer cell lines")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(projectscoredb))
  oncoEnrichR:::validate_db_df(projectscoredb, dbtype = "projectscoredb")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  crispr_lof_results <- list()
  crispr_lof_results[["hits_df"]] <- data.frame()
  crispr_lof_results[["n_genes_with_hits"]] <- 0

  crispr_lof_hits <- as.data.frame(
    target_genes %>%
    dplyr::inner_join(projectscoredb, by = c("symbol"))
  )
  if (nrow(crispr_lof_hits) > 0) {

    crispr_lof_hits <- as.data.frame(
      crispr_lof_hits %>%
      dplyr::select(symbol, symbol_link_ps, model_name,
                    tissue, model_name, model_link_cmp) %>%
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

    crispr_lof_results[["n_genes_with_hits"]] <- nrow(total)
    crispr_lof_results[["hits_df"]] <- crispr_lof_hits %>%
      dplyr::left_join(total, by = c("symbol"))
  }

  return(crispr_lof_results)
}
