
annotate_subcellular_compartments <-
  function(qgenes,
           minimum_confidence = 1,
           show_cytosol = F,
           genedb = NULL,
           oeDB = NULL,
           logger = NULL){

    stopifnot(is.character(qgenes))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(logger))
    stopifnot(!is.null(oeDB$subcelldb))
    stopifnot(!is.null(oeDB$subcelldb$comppidb))
    stopifnot(is.numeric(minimum_confidence))
    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(oeDB$subcelldb$comppidb, dbtype = "comppidb")

    log4r_info(logger, "ComPPI: retrieval of subcellular compartments for target set")


    target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F)

    target_compartments <- list()
    target_compartments[["all"]] <- data.frame()
    target_compartments[["grouped"]] <- data.frame()

    target_compartments_all <- as.data.frame(
      dplyr::inner_join(oeDB$subcelldb$comppidb, target_genes, by = c("symbol")) %>%
        dplyr::filter(.data$confidence >= minimum_confidence)
    )

    if(nrow(target_compartments_all) > 0){
      target_compartments_all <- as.data.frame(
        target_compartments_all %>%
          #dplyr::select(-c(go_ontology,uniprot_acc)) %>%
          dplyr::left_join(
            oeDB$subcelldb[['go_gganatogram_map']],
            by = "go_id") %>%
          dplyr::left_join(
            dplyr::select(genedb,
                          .data$symbol,
                          .data$entrezgene,
                          .data$genename),
            by = c("symbol")) %>%
          dplyr::mutate(
            genelink =
              paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                     .data$entrezgene,"' target='_blank'>", .data$symbol, "</a>")) %>%
          dplyr::mutate(
            compartment =
              paste0("<a href='http://amigo.geneontology.org/amigo/term/",
                     .data$go_id,"' target='_blank'>", .data$go_term, "</a>"))
      )

      target_compartments[["grouped"]] <- as.data.frame(
        target_compartments_all %>%
          dplyr::group_by(.data$go_id, .data$go_term, .data$compartment) %>%
          dplyr::summarise(targets = paste(unique(.data$symbol),
                                           collapse = ", "),
                           targetlinks = paste(unique(.data$genelink),
                                               collapse = ", "),
                           n = dplyr::n()) %>%
          dplyr::arrange(dplyr::desc(.data$n)) %>%
          dplyr::ungroup() %>%
          dplyr::select(-c(.data$go_id, .data$go_term))
      )

      target_compartments[["all"]] <- target_compartments_all %>%
        dplyr::select(-c(.data$entrezgene, .data$go_id,
                         .data$go_term, .data$genelink)) %>%
        dplyr::select(.data$symbol,
                      .data$genename,
                      .data$compartment,
                      dplyr::everything())

      n_genes <- length(unique(target_compartments_all$symbol))
      target_compartments[["anatogram"]] <-
        gganatogram::cell_key$cell %>%
        dplyr::select(-.data$value)

      anatogram_values <- as.data.frame(
        target_compartments_all %>%
          dplyr::select(.data$symbol, .data$ggcompartment) %>%
          dplyr::filter(!is.na(.data$ggcompartment)) %>%
          dplyr::distinct() %>%
          dplyr::group_by(.data$ggcompartment) %>%
          dplyr::summarise(n_comp = dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(value = round((.data$n_comp / n_genes) * 100, digits = 4)) %>%
          dplyr::rename(organ = .data$ggcompartment) %>%
          dplyr::select(.data$organ, .data$value)
      )

      target_compartments[["anatogram"]] <-
        target_compartments[["anatogram"]] %>%
        dplyr::left_join(anatogram_values, by = "organ") %>%
        dplyr::mutate(value = dplyr::if_else(
          is.na(.data$value), 0 , as.numeric(.data$value)))


    }
  return(target_compartments)


}
