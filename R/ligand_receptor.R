
annotate_ligand_receptor_interactions <- function(qgenes,
                                                  genedb = NULL,
                                                  ligandreceptordb = NULL){

  oncoEnrichR:::log4r_info("CellChatDB: retrieval of ligand-receptor interactions")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(ligandreceptordb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR:::validate_db_df(ligandreceptordb$db, dbtype = "ligand_receptor")

  query_df <- data.frame("symbol" = qgenes, stringsAsFactors = F)

  all_query_ligand_receptors <- as.data.frame(
    ligandreceptordb$xref %>%
      dplyr::left_join(
        dplyr::filter(dplyr::select(genedb, entrezgene, symbol),
                      !is.na(symbol)),
        by=c("symbol")) %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::inner_join(query_df, by = "symbol")
  )

  query_hits <- list()
  query_hits[['ligand']] <- all_query_ligand_receptors %>%
    dplyr::filter(class == "ligand") %>%
    dplyr::select(interaction_id)
  query_hits[['receptor']] <- all_query_ligand_receptors %>%
    dplyr::filter(class == "receptor") %>%
    dplyr::select(interaction_id)

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
          dplyr::filter(annotation == "Cell-Cell Contact") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(interaction_id, annotation, interaction_members))

    ligand_receptor_results[['secreted_signaling']] <-
      as.data.frame(
        ligandreceptordb$db %>%
          dplyr::filter(annotation == "Secreted Signaling") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(interaction_id, annotation, interaction_members))


    ligand_receptor_results[['ecm_receptor']] <-
      as.data.frame(
        ligandreceptordb$db %>%
          dplyr::filter(annotation == "ECM-Receptor") %>%
          dplyr::inner_join(valid_query_hits, by = "interaction_id")
      ) %>%
      dplyr::select(-c(interaction_id, annotation, interaction_members))

  }

  return(ligand_receptor_results)

}
