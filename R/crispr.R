get_crispr_lof_scores <- function(qgenes, qsource = "symbol", projectscoredb = NULL){

  rlogging::message("Project Score (CRISPR/Cas9 screen): retrieval of genes associated with loss-of-fitness in cancer cell lines")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(projectscoredb))
  oncoEnrichR::validate_db_df(projectscoredb, dbtype = "projectscoredb")

  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F)

  crispr_lof_results <- list()
  crispr_lof_results[['df']] <- data.frame()
  crispr_lof_results[['plot']] <- NULL
  crispr_lof_results[['n_genes_with_hits']] <- 0

  crispr_lof_hits <- as.data.frame(
    target_genes %>%
    dplyr::inner_join(projectscoredb, by = c("symbol"))
  )
  if(nrow(crispr_lof_hits) > 0){

    crispr_lof_hits <- as.data.frame(
      crispr_lof_hits %>%
      dplyr::select(symbol, model_name, tissue, model_name, model_id_cmp) %>%
      dplyr::group_by(symbol,tissue) %>%
      dplyr::summarise(n_gene_tissue = dplyr::n(),
                       cell_lines = paste(model_name, collapse=", "),
                       cmp_link = paste(model_id_cmp, collapse=", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cell_lines = stringr::str_replace_all(cell_lines,"\\.","-"))
    )

    total <- as.data.frame(
      crispr_lof_hits %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(n_gene = sum(n_gene_tissue))
      )

    crispr_lof_results[['n_genes_with_hits']] <- nrow(total)
    #crispr_lof_results[['plot_height']] <- crispr_lof_results[['plot_height']] + as.integer((nrow(total) - 20)/ 8.5)
    crispr_lof_hits <- dplyr::left_join(crispr_lof_hits, total, by=c("symbol"))

    p <- crispr_lof_hits %>%
      dplyr::mutate(symbol = forcats::fct_reorder(symbol, n_gene)) %>%
      ggplot2::ggplot(ggplot2::aes(x = symbol, y = n_gene_tissue, fill = tissue)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::ylab("Number of cell lines with loss-of-fitness in CRISPR/Cas9 drop-out screen") +
      ggplot2::xlab("") +
      ggplot2::scale_fill_manual(values = pals::stepped(15)) +
      #ggplot2::geom_text(data = total, aes(x = symbol, y = n_total + 1, label = n_total)) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(family = "Helvetica", size = 11),
        legend.text = ggplot2::element_text(family = "Helvetica", size = 11),
        axis.title.x = ggplot2::element_text(family = "Helvetica", size = 12),
        legend.title = ggplot2::element_blank(),
        #set thickness of axis ticks
        axis.ticks = ggplot2::element_line(size=0.2),
        #remove plot background
        plot.background = ggplot2::element_blank()
      )

    crispr_lof_results[['df']] <- crispr_lof_hits
    crispr_lof_results[['plot']] <- p
  }

  return(crispr_lof_results)
}
