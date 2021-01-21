
get_genes_unknown_function <- function(qgenes, genedb = NULL,
                                       poorly_defined_genes = NULL){

  rlogging::message("Retrieval of genes with unknown/poorly defined function in target set")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(poorly_defined_genes))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR:::validate_db_df(poorly_defined_genes, dbtype = "pdf")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  results <- data.frame()
  pct <- 0
  if(nrow(target_genes) > 0){
    results <- poorly_defined_genes %>%
      dplyr::inner_join(target_genes, by = "symbol")

    pct <- round(as.numeric(NROW(results) / NROW(target_genes)) * 100, digits = 2)
  }

  rlogging::message("Detected n = ", nrow(results),
                    " (", pct,"%) target genes with unknown/poorly defined function")
  return(results)

}
