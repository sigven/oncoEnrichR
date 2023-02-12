
get_genes_unknown_function <- function(qgenes = NULL,
                                       genedb = NULL) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( "Retrieval of genes with unknown/poorly defined function in target set")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)
  poorly_defined_genes <- genedb |>
    dplyr::filter(!is.na(.data$unknown_function_rank)) |>
    dplyr::select(c("symbol",
                  "genename",
                  "num_go_terms",
                  "go_term_link",
                  "gene_summary",
                  "unknown_function_rank",
                  "has_gene_summary"))


  results <- data.frame()
  pct <- 0
  if (nrow(target_genes) > 0) {
    results <- poorly_defined_genes |>
      dplyr::inner_join(
        target_genes, by = "symbol", multiple = "all") |>
      dplyr::arrange(.data$unknown_function_rank) |>
      dplyr::filter(.data$unknown_function_rank <= 6)

    pct <- round(as.numeric(NROW(results) / NROW(target_genes)) * 100, digits = 2)
  }

  lgr::lgr$info( paste0("Detected n = ", nrow(results),
                    " (", pct,"%) target genes with unknown/poorly defined function"))
  return(results)

}
