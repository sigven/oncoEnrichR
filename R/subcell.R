
annotate_subcellular_compartments <-
  function(query_entrez,
           compartments_min_confidence = 3,
           compartments_min_channels = 1,
           show_cytosol = F,
           genedb = NULL,
           compartments = NULL){

    stopifnot(!is.null(query_entrez))
    stopifnot(is.integer(query_entrez))
    stopifnot(!is.null(genedb))
    stopifnot(!is.null(compartments))
    stopifnot(is.numeric(compartments_min_confidence))
    stopifnot(is.numeric(compartments_min_channels))

    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(compartments, dbtype = "compartments")
    lgr::lgr$info( "COMPARTMENTS: retrieval of subcellular compartments for target set")

    target_genes <- data.frame(
      'entrezgene' = query_entrez, stringsAsFactors = F) |>
      dplyr::left_join(
        dplyr::select(genedb,
                      c("symbol",
                      "entrezgene",
                      "cancer_max_rank",
                      "genename")),
        by = c("entrezgene"), relationship = "many-to-many")

    target_compartments <- list()
    target_compartments[["all"]] <- data.frame()
    target_compartments[["grouped"]] <- data.frame()
    target_compartments[["comp_density"]] <- data.frame()

    target_compartments_all <- as.data.frame(
      dplyr::inner_join(
        compartments, target_genes,
        by = c("entrezgene"), relationship = "many-to-many") |>
        dplyr::filter(!is.na(.data$symbol)) |>
        dplyr::filter(.data$confidence >= compartments_min_confidence) |>
        dplyr::mutate(
          channel_confidence = paste(
            .data$annotation_channel,
            .data$confidence, sep =":"
          )
        ) |>
        dplyr::group_by(
          .data$entrezgene,
          .data$symbol,
          .data$cancer_max_rank,
          .data$genename,
          .data$go_id,
          .data$go_term
        ) |>
        dplyr::summarise(
          maximum_confidence = max(.data$confidence),
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
          dplyr::desc(.data$cancer_max_rank),
          .data$symbol,
          dplyr::desc(.data$maximum_confidence),
          dplyr::desc(.data$n_supporting_channels)
        )
    )

    if (nrow(target_compartments_all) > 0) {
      target_compartments_all <- as.data.frame(
        target_compartments_all |>
          dplyr::arrange(
            dplyr::desc(.data$maximum_confidence),
            dplyr::desc(.data$cancer_max_rank),
            .data$symbol) |>
          dplyr::left_join(
            dplyr::select(
              oncoEnrichR::subcell_map$map,
              c("id","name","go_id","subcellular_location_id")),
                by = "go_id",
            relationship = "many-to-many"
          ) |>
          ## ignore extracellular space and secreted
          dplyr::filter(
            .data$id != "SL0112" &
              .data$id != "SL0243"
          ) |>
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
          dplyr::arrange(
            .data$go_id,
            .data$go_term,
            .data$compartment,
            dplyr::desc(.data$cancer_max_rank),
          ) |>
          dplyr::group_by(
            .data$go_id, .data$go_term, .data$compartment) |>
          dplyr::summarise(
            targets = paste(head(unique(.data$symbol),75),
                            collapse = ", "),
            targetlinks = paste(head(unique(.data$genelink),75),
                                collapse = ", "),
            n = dplyr::n(), .groups = "drop") |>
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

      target_compartments[['comp_density']] <- as.data.frame(
        target_compartments_all |>
          dplyr::select(c("symbol", "id","name")) |>
          dplyr::filter(!is.na(.data$id)) |>
          dplyr::distinct() |>
          dplyr::group_by(.data$id, .data$name) |>
          dplyr::summarise(
            n_comp = dplyr::n(),
            genes = paste(.data$symbol, collapse=", "),
            .groups = "drop") |>
          dplyr::ungroup() |>
          dplyr::filter(
            stringr::str_detect(.data$genes,",")
          ) |>
          dplyr::mutate(proportion = round(
            (.data$n_comp / n_genes), digits = 4)) |>
          dplyr::arrange(dplyr::desc(.data$proportion))
      )

      target_compartments[['comp_density']]$bin <-
        cut(
          target_compartments[['comp_density']]$proportion,
          breaks = seq(0, 1, by = 0.1),
          include.lowest = TRUE, labels = FALSE
        )

      target_compartments[['comp_density']] <-
        target_compartments[['comp_density']] |>
        dplyr::left_join(
          oncoEnrichR::color_palette$subcell_compartments,
          by = "bin"
        )

    }
  return(target_compartments)


}
