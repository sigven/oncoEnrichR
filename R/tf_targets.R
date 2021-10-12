
annotate_tf_targets <- function(qgenes,
                                genedb = NULL,
                                tf_target_interactions = NULL,
                                collection = "global",
                                min_confidence_reg_interaction = "D"){

  stopifnot(is.character(collection))
  stopifnot(collection == "global" | collection == "pancancer")
  oncoEnrichR:::log4r_info(
    paste0("DoRothEA: retrieval of regulatory interactions involving members of target set - ", collection))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(tf_target_interactions[[collection]]))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR:::validate_db_df(tf_target_interactions[[collection]], dbtype = "dorothea")

  exclude_level_regex <- "E"
  if(min_confidence_reg_interaction == "C"){
    exclude_level_regex <- "D|E"
  }
  if(min_confidence_reg_interaction == "B"){
    exclude_level_regex <- "C|D|E"
  }
  if(min_confidence_reg_interaction == "A"){
    exclude_level_regex <- "B|C|D|E"
  }

  target_genes <- data.frame("target" = qgenes, stringsAsFactors = F) %>%
    dplyr::left_join(tf_target_interactions[[collection]],
                     by = c("target")) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(confidence_level)) %>%
    dplyr::filter(!stringr::str_detect(confidence_level, exclude_level_regex))

  if(nrow(target_genes) > 0){
    target_genes <- target_genes %>%
      dplyr::mutate(queryset_overlap = paste("TARGET", confidence_level, sep = "_"))
  }

  tf_genes <- data.frame("regulator" = qgenes, stringsAsFactors = F) %>%
    dplyr::left_join(tf_target_interactions[[collection]],
                     by = c("regulator")) %>%
    dplyr::distinct() %>%
    dplyr::filter(!stringr::str_detect(confidence_level, exclude_level_regex))

  if(nrow(tf_genes) > 0){
    tf_genes <- tf_genes %>%
      dplyr::mutate(queryset_overlap = paste("TF", confidence_level, sep = "_"))
  }

  query_tf_target_interactions <- data.frame()

  if(nrow(target_genes) > 0 | nrow(tf_genes) > 0){
    query_tf_target_interactions <- as.data.frame(
      tf_genes %>%
        dplyr::bind_rows(target_genes) %>%
        dplyr::group_by(regulator, target, mode_of_regulation,
                        interaction_sources, confidence_level,
                        tf_target_literature_support) %>%
        dplyr::summarise(queryset_overlap = paste(queryset_overlap, collapse=";"),
                         .groups = "drop") %>%
        dplyr::left_join(dplyr::select(oncoEnrichR::genedb$all,
                                       symbol, genename),
                         by = c("regulator" = "symbol")) %>%
        dplyr::rename(regulator_name = genename) %>%
        dplyr::left_join(dplyr::select(oncoEnrichR::genedb$all,
                                       symbol, genename),
                         by = c("target" = "symbol")) %>%
        dplyr::rename(target_name = genename) %>%
        dplyr::mutate(queryset_overlap = stringr::str_replace(
          queryset_overlap, "(A|B|C|D);","")) %>%
        dplyr::mutate(queryset_overlap = factor(
          queryset_overlap, levels = c("TF_TARGET_A","TF_TARGET_B","TF_TARGET_C",
                                       "TF_TARGET_D","TF_A","TF_B","TF_C","TF_D",
                                       "TARGET_A","TARGET_B","TARGET_C","TARGET_D"))) %>%
        dplyr::arrange(queryset_overlap, confidence_level) %>%
        dplyr::rename(literature_support = tf_target_literature_support) %>%
        dplyr::select(regulator, regulator_name, target, target_name,
                      confidence_level, mode_of_regulation,
                      literature_support, interaction_sources,
                      queryset_overlap)
    )
  }


  return(query_tf_target_interactions)

}

retrieve_tf_target_network <- function(tf_target_interactions = NULL){


  tf_target_network <- list()
  tf_target_network[['nodes']] <- data.frame()
  tf_target_network[['edges']] <- data.frame()

  if(NROW(tf_target_interactions) > 0){

    complete_interactions <- tf_target_interactions %>%
      dplyr::filter(stringr::str_detect(queryset_overlap,"TF_TARGET_"))

    if(NROW(complete_interactions) > 0){

      tf_target_network[['edges']] <- complete_interactions %>%
        dplyr::select(regulator, target, confidence_level, mode_of_regulation) %>%
        dplyr::rename(from = regulator, to = target) %>%
        dplyr::mutate(title = paste0(mode_of_regulation, ": confidence level ",
                                     confidence_level)) %>%
        dplyr::mutate(arrows = "to") %>%
        dplyr::mutate(value = dplyr::case_when(
          confidence_level == "A" ~ 2,
          confidence_level == "B" ~ 1.75,
          confidence_level == "C" ~ 1.5,
          confidence_level == "D" ~ 1.25
        )) %>%
        dplyr::mutate(group = dplyr::case_when(
          mode_of_regulation == "Stimulation" ~ "A",
          mode_of_regulation == "Repression" ~ "B"
        )) %>%
        dplyr::arrange(confidence_level) %>%
        head(100)

      all_nodes <-
      #tf_target_network[['nodes']] <-
        data.frame('id' = complete_interactions$regulator,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$regulator,
                   'title' = stringr::str_trim(
                     textclean::replace_html(complete_interactions$regulator_name)
                   ),
                   stringsAsFactors = F) %>%
      dplyr::bind_rows(
        data.frame('id' = complete_interactions$target,
                   'shadow' = F,
                   'size' = 25,
                   'label' = complete_interactions$target,
                   'title' = stringr::str_trim(
                     textclean::replace_html(complete_interactions$target_name)
                   ),
                   stringsAsFactors = F)) %>%
        dplyr::distinct()

      nodeset1 <- all_nodes %>%
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], from),
          by = c("id" = "from")) %>%
        dplyr::distinct()

      nodeset2 <- all_nodes %>%
        dplyr::inner_join(
          dplyr::select(tf_target_network[['edges']], to),
          by = c("id" = "to")) %>%
        dplyr::distinct()

      tf_target_network[['nodes']] <-
        nodeset1 %>%
        dplyr::bind_rows(nodeset2) %>%
        dplyr::distinct()

    }

  }

  return(tf_target_network)

}
