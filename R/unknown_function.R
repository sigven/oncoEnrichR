
get_genes_unknown_function <- function(qgenes,
                                       genedb = NULL,
                                       logger = NULL){

  stopifnot(!is.null(logger))
  log4r_info(logger, "Retrieval of genes with unknown/poorly defined function in target set")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)
  poorly_defined_genes <- genedb %>%
    dplyr::filter(!is.na(.data$unknown_function_rank)) %>%
    dplyr::select(.data$symbol,
                  .data$genename,
                  .data$num_go_terms,
                  .data$gene_summary,
                  .data$unknown_function_rank,
                  .data$has_gene_summary)


  results <- data.frame()
  pct <- 0
  if(nrow(target_genes) > 0){
    results <- poorly_defined_genes %>%
      dplyr::inner_join(target_genes, by = "symbol") %>%
      dplyr::arrange(.data$unknown_function_rank) %>%
      dplyr::filter(.data$unknown_function_rank <= 5)

    pct <- round(as.numeric(NROW(results) / NROW(target_genes)) * 100, digits = 2)
  }

  log4r_info(logger, paste0("Detected n = ", nrow(results),
                    " (", pct,"%) target genes with unknown/poorly defined function"))
  return(results)

}
