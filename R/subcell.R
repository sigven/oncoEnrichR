
annotate_subcellular_compartments <-
  function(query_entrez,
           compartments_min_confidence = 3,
           compartments_min_channels = 1,
           compartments_channel_types =
             c("Experimental","Text mining","Knowledge"),
           show_cytosol = F,
           genedb = NULL,
           compartments = NULL,
           go_gganatogram_map = NULL) {

    stopifnot(!is.null(query_entrez))
    stopifnot(is.integer(query_entrez))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(compartments))
    stopifnot(!is.null(go_gganatogram_map))
    stopifnot(is.numeric(compartments_min_confidence))
    stopifnot(is.numeric(compartments_min_channels))

    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(compartments, dbtype = "compartments")

    validate_db_df(go_gganatogram_map, dbtype = "go_gganatogram")

    lgr::lgr$info( "COMPARTMENTS: retrieval of subcellular compartments for target set")

    target_genes <- data.frame(
      'entrezgene' = query_entrez, stringsAsFactors = F) |>
      dplyr::left_join(
        dplyr::select(genedb,
                      c("symbol",
                      "entrezgene",
                      "genename")),
        by = c("entrezgene"), multiple = "all")

    target_compartments <- list()
    target_compartments[["all"]] <- data.frame()
    target_compartments[["grouped"]] <- data.frame()

    target_compartments_all <- as.data.frame(
      dplyr::inner_join(
        compartments, target_genes, by = c("entrezgene"), multiple = "all") |>
        dplyr::filter(!is.na(.data$symbol)) |>
        dplyr::filter(.data$confidence >= compartments_min_confidence) |>
        dplyr::filter(
          .data$annotation_channel %in% compartments_channel_types
        ) |>
        dplyr::mutate(
          channel_confidence = paste(
            .data$annotation_channel,
            .data$confidence, sep =":"
          )
        ) |>
        dplyr::group_by(
          .data$entrezgene,
          .data$symbol,
          .data$genename,
          .data$go_id,
          .data$go_term
        ) |>
        dplyr::summarise(
          minimum_confidence = min(.data$confidence),
          supporting_channels = paste(
            sort(unique(.data$annotation_channel)),
            collapse ="|"),
          supporting_channels_confidence = paste(
            sort(unique(.data$channel_confidence)),
            collapse =", "),
          supporting_sources = paste(
            sort(unique(.data$annotation_source)),
            collapse = ", "
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          n_supporting_channels = stringr::str_count(
            .data$supporting_channels, "\\|") + 1
        ) |>
        dplyr::filter(
          .data$n_supporting_channels >= compartments_min_channels
        ) |>
        dplyr::arrange(
          .data$symbol,
          dplyr::desc(.data$minimum_confidence),
          dplyr::desc(.data$n_supporting_channels)
        )
    )

    if (nrow(target_compartments_all) > 0) {
      target_compartments_all <- as.data.frame(
        target_compartments_all |>
          dplyr::arrange(
            dplyr::desc(.data$minimum_confidence),
            .data$symbol) |>
          dplyr::left_join(
            go_gganatogram_map,
            by = "go_id", multiple = "all") |>
          dplyr::mutate(
            genelink =
              paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                     .data$entrezgene,"' target='_blank'>", .data$symbol, "</a>")) |>
          dplyr::mutate(
            compartment =
              paste0("<a href='http://amigo.geneontology.org/amigo/term/",
                     .data$go_id,"' target='_blank'>", .data$go_term, "</a>")) |>
          utils::head(2500)
      )

      target_compartments[["grouped"]] <- as.data.frame(
        target_compartments_all |>
          dplyr::group_by(.data$go_id, .data$go_term, .data$compartment) |>
          dplyr::summarise(
            targets = paste(unique(.data$symbol),
                            collapse = ", "),
            targetlinks = paste(unique(.data$genelink),
                                collapse = ", "),
            n = dplyr::n()) |>
          dplyr::arrange(dplyr::desc(.data$n)) |>
          dplyr::ungroup() |>
          dplyr::select(-c("go_id", "go_term"))
      )

      target_compartments[["all"]] <- target_compartments_all |>
        dplyr::select(-c("entrezgene", "go_id",
                         "go_term", "genelink")) |>
        dplyr::select(c("symbol",
                       "genename",
                       "compartment"),
                      dplyr::everything())

      n_genes <- length(unique(target_compartments_all$symbol))
      target_compartments[["anatogram"]] <-
        gganatogram::cell_key$cell |>
        dplyr::select(-c("value"))

      anatogram_values <- as.data.frame(
        target_compartments_all |>
          dplyr::select(c("symbol", "ggcompartment")) |>
          dplyr::filter(!is.na(.data$ggcompartment)) |>
          dplyr::distinct() |>
          dplyr::group_by(.data$ggcompartment) |>
          dplyr::summarise(n_comp = dplyr::n()) |>
          dplyr::ungroup() |>
          dplyr::mutate(value = round((.data$n_comp / n_genes) * 100, digits = 4)) |>
          dplyr::rename(organ = "ggcompartment") |>
          dplyr::select(c("organ","value"))
      )

      target_compartments[["anatogram"]] <-
        target_compartments[["anatogram"]] |>
        dplyr::left_join(anatogram_values, by = "organ", multiple = "all") |>
        dplyr::mutate(value = dplyr::if_else(
          is.na(.data$value), 0 , as.numeric(.data$value)))


    }
  return(target_compartments)


}
