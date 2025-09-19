
annotate_tf_targets <- function(qgenes,
                                genedb = NULL,
                                tf_target_interactions = NULL,
                                regulatory_min_resources = 2) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(tf_target_interactions))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tf_target_interactions, dbtype = "collectri")
  stopifnot(is.numeric(regulatory_min_resources))
  stopifnot(regulatory_min_resources >= 1 &
              regulatory_min_resources <= 3)

  lgr::lgr$info(
    paste0("Collectri: retrieval of regulatory interactions involving members of target set"))

  target_genes <- data.frame("target" = qgenes, stringsAsFactors = F) |>
    dplyr::left_join(tf_target_interactions,
                     by = c("target"), relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$n_resources) &
                  .data$n_resources >= regulatory_min_resources)

  if (nrow(target_genes) > 0) {
    target_genes <- target_genes |>
      dplyr::mutate(
        queryset_overlap = paste("TARGET",
                                 .data$n_resources, sep = "_"))
  }

  tf_genes <- data.frame("regulator" = qgenes, stringsAsFactors = F) |>
    dplyr::left_join(tf_target_interactions,
                     by = c("regulator"),
                     relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$n_resources) &
                    .data$n_resources >= regulatory_min_resources)

  if (nrow(tf_genes) > 0) {
    tf_genes <- tf_genes |>
      dplyr::mutate(
        queryset_overlap = paste("TF", .data$n_resources, sep = "_"))
  }

  query_tf_target_interactions <- data.frame()

  if (nrow(target_genes) > 0 | nrow(tf_genes) > 0) {
    query_tf_target_interactions <- as.data.frame(
      tf_genes |>
        dplyr::bind_rows(target_genes) |>
        dplyr::group_by(.data$regulator, .data$target, .data$mode_of_regulation,
                        .data$interaction_sources, .data$n_resources,
                        .data$tf_target_literature_support) |>
        dplyr::summarise(queryset_overlap = paste(.data$queryset_overlap, collapse = ";"),
                         .groups = "drop") |>
        dplyr::left_join(
          dplyr::select(genedb,
                        c("symbol", "genename",
                        "cancer_max_rank")),
          by = c("regulator" = "symbol"), relationship = "many-to-many") |>
        dplyr::rename(regulator_name = "genename",
                      regulator_cancer_max_rank = "cancer_max_rank") |>
        dplyr::left_join(
          dplyr::select(genedb,
                        c("symbol", "genename",
                          "cancer_max_rank")),
          by = c("target" = "symbol"), relationship = "many-to-many") |>
        dplyr::rename(target_name = "genename",
                      target_cancer_max_rank = "cancer_max_rank") |>
        dplyr::mutate(queryset_overlap = stringr::str_replace(
          .data$queryset_overlap, "(0|1|2|3|4|5);","")) |>
        dplyr::mutate(queryset_overlap = factor(
          .data$queryset_overlap,
          levels = c("TF_TARGET_1","TF_TARGET_2","TF_TARGET_3",
                     "TF_TARGET_4","TF_TARGET_5","TF_TARGET_6",
                      "TF_1","TF_2","TF_3","TF_4","TF_5","TF_6",
                     "TARGET_1","TARGET_2","TARGET_3",
                      "TARGET_4","TARGET_5","TARGET_6"))) |>
        dplyr::arrange(.data$queryset_overlap, .data$n_resources,
                       dplyr::desc(.data$regulator_cancer_max_rank),
                       dplyr::desc(.data$target_cancer_max_rank)) |>
        dplyr::rename(literature_support = "tf_target_literature_support") |>
        dplyr::select(c("regulator",
                        "regulator_name",
                        "target", "target_name",
                        "n_resources",
                        "mode_of_regulation",
                        "literature_support",
                        "interaction_sources",
                        "queryset_overlap",
                        "regulator_cancer_max_rank",
                        "target_cancer_max_rank"))
    )
  }


  return(query_tf_target_interactions)

}

retrieve_tf_target_network <- function(tf_target_interactions = NULL) {


  stopifnot(!is.null(tf_target_interactions))
  validate_db_df(tf_target_interactions, dbtype = "tf_target_interactions")

  tf_target_network <- list()
  tf_target_network[['nodes']] <- data.frame()
  tf_target_network[['edges']] <- data.frame()

  if (NROW(tf_target_interactions) > 0) {

    complete_interactions <- tf_target_interactions |>
      dplyr::filter(stringr::str_detect(.data$queryset_overlap,"TF_TARGET_")) |>
      utils::head(200)

    if (NROW(complete_interactions) > 0) {

      tf_target_network[['edges']] <- complete_interactions |>
        dplyr::select(c("regulator",
                        "target",
                        "n_resources",
                        "mode_of_regulation")) |>
        dplyr::rename(from = "regulator",
                      to = "target") |>
        dplyr::mutate(title = paste0(.data$mode_of_regulation, ": num resources ",
                                     .data$n_resources)) |>
        dplyr::mutate(arrows = "to") |>
        dplyr::mutate(length = dplyr::case_when(
          .data$n_resources == 6 | .data$n_resources == 5 ~ 150,
          .data$n_resources == 4 ~ 200,
          .data$n_resources == 3 ~ 250,
          .data$n_resources == 2 ~ 300,
          .data$n_resources == 1 ~ 200,
          TRUE ~ 150
        )) |>
        dplyr::mutate(dashes = T) |>
        dplyr::mutate(color = dplyr::case_when(
          .data$mode_of_regulation == "Stimulation" ~ "darkgreen",
          .data$mode_of_regulation == "Repression" ~ "darkred",
         )) |>
        dplyr::arrange(.data$n_resources) |>
        utils::head(150)

      all_nodes <-
      #tf_target_network[['nodes']] <-
        data.frame('id' = complete_interactions$regulator,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$regulator,
                   'title' = stringr::str_trim(
                     textclean::replace_html(
                       complete_interactions$regulator_name)
                   ),
                   stringsAsFactors = F) |>
      dplyr::bind_rows(
        data.frame('id' = complete_interactions$target,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$target,
                   'title' = stringr::str_trim(
                     textclean::replace_html(
                       complete_interactions$target_name)
                   ),
                   stringsAsFactors = F)) |>
        dplyr::distinct()

      nodeset1 <- all_nodes |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], c("from")),
          by = c("id" = "from"), relationship = "many-to-many") |>
        dplyr::distinct()

      nodeset2 <- all_nodes |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], c("to")),
          by = c("id" = "to"), relationship = "many-to-many") |>
        dplyr::distinct()

      tf_target_network[['nodes']] <-
        nodeset1 |>
        dplyr::bind_rows(nodeset2) |>
        dplyr::distinct()


      regulators <- tf_target_network[['edges']] |>
        dplyr::select(c("from")) |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['nodes']], c("id")),
          by = c("from" = "id"), relationship = "many-to-many") |>
        dplyr::rename(id = "from") |>
        dplyr::distinct()

      targets <- tf_target_network[['edges']] |>
        dplyr::select(c("to")) |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['nodes']], c("id")),
          by = c("to" = "id"), relationship = "many-to-many") |>
        dplyr::rename(id = "to") |>
        dplyr::distinct()

      target_and_regulator <- targets |>
        dplyr::inner_join(regulators, by = "id", relationship = "many-to-many")

      if (nrow(target_and_regulator) > 0) {
        regulators <- regulators |>
          dplyr::anti_join(target_and_regulator, by = "id") |>
          dplyr::mutate(group = "Regulator")

        targets <- targets |>
          dplyr::anti_join(target_and_regulator, by = "id") |>
          dplyr::mutate(group = "Target")

        target_and_regulator <- target_and_regulator |>
          dplyr::mutate(group = "Regulator/target")
      } else {
        regulators <- regulators |>
          dplyr::mutate(group = "Regulator")

        targets <- targets |>
          dplyr::mutate(group = "Target")
      }

      all_shapes <- regulators |>
        dplyr::bind_rows(targets) |>
        dplyr::bind_rows(target_and_regulator)

      tf_target_network[['nodes']] <- tf_target_network[['nodes']] |>
        dplyr::left_join(all_shapes, by = "id", relationship = "many-to-many")

    }

  }

  return(tf_target_network)

}
