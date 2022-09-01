
annotate_tf_targets <- function(qgenes,
                                genedb = NULL,
                                tf_target_interactions = NULL,
                                collection = "global",
                                min_confidence_reg_interaction = "D"){

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(is.character(collection))
  stopifnot(collection == "global" | collection == "pancancer")

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(typeof(tf_target_interactions) == "list")
  stopifnot(!is.null(tf_target_interactions[[collection]]))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(tf_target_interactions[[collection]], dbtype = "dorothea")
  stopifnot(min_confidence_reg_interaction %in% c("A","B","C","D"))

  lgr::lgr$info(
    paste0("DoRothEA: retrieval of regulatory interactions involving members of target set - ", collection))

  exclude_level_regex <- "E"
  if(min_confidence_reg_interaction == "C"){
    exclude_level_regex <- "(D|E)"
  }
  if(min_confidence_reg_interaction == "B"){
    exclude_level_regex <- "(C|D|E)"
  }
  if(min_confidence_reg_interaction == "A"){
    exclude_level_regex <- "(B|C|D|E)"
  }

  target_genes <- data.frame("target" = qgenes, stringsAsFactors = F) |>
    dplyr::left_join(tf_target_interactions[[collection]],
                     by = c("target")) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$confidence_level)) |>
    dplyr::filter(!stringr::str_detect(
      .data$confidence_level, exclude_level_regex))

  if(nrow(target_genes) > 0){
    target_genes <- target_genes |>
      dplyr::mutate(
        queryset_overlap = paste("TARGET",
                                 .data$confidence_level, sep = "_"))
  }

  tf_genes <- data.frame("regulator" = qgenes, stringsAsFactors = F) |>
    dplyr::left_join(tf_target_interactions[[collection]],
                     by = c("regulator")) |>
    dplyr::distinct() |>
    dplyr::filter(!stringr::str_detect(
      .data$confidence_level, exclude_level_regex))

  if(nrow(tf_genes) > 0){
    tf_genes <- tf_genes |>
      dplyr::mutate(queryset_overlap = paste("TF", .data$confidence_level, sep = "_"))
  }

  query_tf_target_interactions <- data.frame()

  if(nrow(target_genes) > 0 | nrow(tf_genes) > 0){
    query_tf_target_interactions <- as.data.frame(
      tf_genes |>
        dplyr::bind_rows(target_genes) |>
        dplyr::group_by(.data$regulator, .data$target, .data$mode_of_regulation,
                        .data$interaction_sources, .data$confidence_level,
                        .data$tf_target_literature_support) |>
        dplyr::summarise(queryset_overlap = paste(.data$queryset_overlap, collapse=";"),
                         .groups = "drop") |>
        dplyr::left_join(dplyr::select(genedb,
                                       .data$symbol, .data$genename,
                                       .data$cancer_max_rank),
                         by = c("regulator" = "symbol")) |>
        dplyr::rename(regulator_name = .data$genename,
                      regulator_cancer_max_rank = .data$cancer_max_rank) |>
        dplyr::left_join(dplyr::select(genedb,
                                       .data$symbol, .data$genename,
                                       .data$cancer_max_rank),
                         by = c("target" = "symbol")) |>
        dplyr::rename(target_name = .data$genename,
                      target_cancer_max_rank = .data$cancer_max_rank) |>
        dplyr::mutate(queryset_overlap = stringr::str_replace(
          .data$queryset_overlap, "(A|B|C|D);","")) |>
        dplyr::mutate(queryset_overlap = factor(
          .data$queryset_overlap, levels = c("TF_TARGET_A","TF_TARGET_B","TF_TARGET_C",
                                       "TF_TARGET_D","TF_A","TF_B","TF_C","TF_D",
                                       "TARGET_A","TARGET_B","TARGET_C","TARGET_D"))) |>
        dplyr::arrange(.data$queryset_overlap, .data$confidence_level,
                       dplyr::desc(.data$regulator_cancer_max_rank),
                       dplyr::desc(.data$target_cancer_max_rank)) |>
        dplyr::rename(literature_support = .data$tf_target_literature_support) |>
        dplyr::select(.data$regulator, .data$regulator_name,
                      .data$target, .data$target_name,
                      .data$confidence_level, .data$mode_of_regulation,
                      .data$literature_support, .data$interaction_sources,
                      .data$queryset_overlap, .data$regulator_cancer_max_rank,
                      .data$target_cancer_max_rank)
    )
  }


  return(query_tf_target_interactions)

}

retrieve_tf_target_network <- function(tf_target_interactions = NULL){


  stopifnot(!is.null(tf_target_interactions))
  validate_db_df(tf_target_interactions, dbtype = "tf_target_interactions")

  tf_target_network <- list()
  tf_target_network[['nodes']] <- data.frame()
  tf_target_network[['edges']] <- data.frame()

  if(NROW(tf_target_interactions) > 0){

    complete_interactions <- tf_target_interactions |>
      dplyr::filter(stringr::str_detect(.data$queryset_overlap,"TF_TARGET_")) |>
      utils::head(200)

    if(NROW(complete_interactions) > 0){

      tf_target_network[['edges']] <- complete_interactions |>
        dplyr::select(.data$regulator, .data$target,
                      .data$confidence_level, .data$mode_of_regulation) |>
        dplyr::rename(from = .data$regulator,
                      to = .data$target) |>
        dplyr::mutate(title = paste0(.data$mode_of_regulation, ": confidence level ",
                                     .data$confidence_level)) |>
        dplyr::mutate(arrows = "to") |>
        dplyr::mutate(length = dplyr::case_when(
          .data$confidence_level == "A" ~ 150,
          .data$confidence_level == "B" ~ 200,
          .data$confidence_level == "C" ~ 250,
          .data$confidence_level == "D" ~ 300,
        )) |>
        dplyr::mutate(dashes = T) |>
        dplyr::mutate(color = dplyr::case_when(
          .data$mode_of_regulation == "Stimulation" ~ "darkgreen",
          .data$mode_of_regulation == "Repression" ~ "darkred",
         )) |>
        dplyr::arrange(.data$confidence_level) |>
        utils::head(150)

      all_nodes <-
      #tf_target_network[['nodes']] <-
        data.frame('id' = complete_interactions$regulator,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$regulator,
                   'title' = stringr::str_trim(
                     textclean::replace_html(complete_interactions$regulator_name)
                   ),
                   stringsAsFactors = F) |>
      dplyr::bind_rows(
        data.frame('id' = complete_interactions$target,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$target,
                   'title' = stringr::str_trim(
                     textclean::replace_html(complete_interactions$target_name)
                   ),
                   stringsAsFactors = F)) |>
        dplyr::distinct()

      nodeset1 <- all_nodes |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], .data$from),
          by = c("id" = "from")) |>
        dplyr::distinct()

      nodeset2 <- all_nodes |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], .data$to),
          by = c("id" = "to")) |>
        dplyr::distinct()

      tf_target_network[['nodes']] <-
        nodeset1 |>
        dplyr::bind_rows(nodeset2) |>
        dplyr::distinct()


      regulators <- tf_target_network[['edges']] |>
        dplyr::select(.data$from) |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['nodes']], .data$id),
          by = c("from" = "id")
          ) |>
        dplyr::rename(id = .data$from) |>
        dplyr::distinct()

      targets <- tf_target_network[['edges']] |>
        dplyr::select(.data$to) |>
        dplyr::inner_join(
          dplyr::select(tf_target_network[['nodes']], .data$id),
          by = c("to" = "id")
        ) |>
        dplyr::rename(id = .data$to) |>
        dplyr::distinct()

      target_and_regulator <- targets |>
        dplyr::inner_join(regulators, by = "id")

      if(nrow(target_and_regulator) > 0){
        regulators <- regulators |>
          dplyr::anti_join(target_and_regulator, by = "id") |>
          dplyr::mutate(group = "Regulator")

        targets <- targets |>
          dplyr::anti_join(target_and_regulator, by = "id") |>
          dplyr::mutate(group = "Target")

        target_and_regulator <- target_and_regulator |>
          dplyr::mutate(group = "Regulator/target")
      }else{
        regulators <- regulators |>
          dplyr::mutate(group = "Regulator")

        targets <- targets |>
          dplyr::mutate(group = "Target")
      }

      all_shapes <- regulators |>
        dplyr::bind_rows(targets) |>
        dplyr::bind_rows(target_and_regulator)

      tf_target_network[['nodes']] <- tf_target_network[['nodes']] |>
        dplyr::left_join(all_shapes, by = "id")

    }

  }

  return(tf_target_network)

}
