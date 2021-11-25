
annotate_ligand_receptor_interactions <- function(qgenes,
                                                  genedb = NULL,
                                                  ligandreceptordb = NULL){

  log4r_info("CellChatDB: retrieval of ligand-receptor interactions")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(ligandreceptordb))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(ligandreceptordb$db, dbtype = "ligand_receptor")

  query_df <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  all_query_ligand_receptors <- as.data.frame(
    ligandreceptordb$xref %>%
      dplyr::left_join(
        dplyr::filter(dplyr::select(genedb,
                                    .data$entrezgene,
                                    .data$symbol),
                      !is.na(.data$symbol)),
        by=c("symbol")) %>%
      dplyr::filter(!is.na(.data$entrezgene)) %>%
      dplyr::inner_join(query_df, by = "symbol")
  )

  query_hits <- list()
  query_hits[['ligand']] <- all_query_ligand_receptors %>%
    dplyr::filter(class == "ligand") %>%
    dplyr::select(.data$interaction_id)
  query_hits[['receptor']] <- all_query_ligand_receptors %>%
    dplyr::filter(class == "receptor") %>%
    dplyr::select(.data$interaction_id)

  valid_query_hits <-
    query_hits[['ligand']] %>%
    dplyr::inner_join(query_hits[['receptor']],
                      by = "interaction_id")

  ligand_receptor_results <- list()
  ligand_receptor_results[['cell_cell_contact']] <-
    data.frame()
  ligand_receptor_results[['secreted_signaling']] <-
    data.frame()
  ligand_receptor_results[['ecm_receptor']] <-
    data.frame()

  if(NROW(valid_query_hits) > 0){
    ligand_receptor_results[['cell_cell_contact']] <-
      as.data.frame(
        ligandreceptordb$db %>%
          dplyr::filter(.data$annotation == "Cell-Cell Contact") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(.data$interaction_id,
                       .data$annotation,
                       .data$interaction_members)) %>%
      dplyr::distinct()

    ligand_receptor_results[['secreted_signaling']] <-
      as.data.frame(
        ligandreceptordb$db %>%
          dplyr::filter(.data$annotation == "Secreted Signaling") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(.data$interaction_id,
                       .data$annotation,
                       .data$interaction_members)) %>%
      dplyr::distinct()


    ligand_receptor_results[['ecm_receptor']] <-
      as.data.frame(
        ligandreceptordb$db %>%
          dplyr::filter(.data$annotation == "ECM-Receptor") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(.data$interaction_id,
                       .data$annotation,
                       .data$interaction_members)) %>%
      dplyr::distinct()

  }

  return(ligand_receptor_results)

}
