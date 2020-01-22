
annotate_subcellular_compartments <- function(qgenes, genedb = NULL, comppidb = NULL){

  rlogging::message("ComPPI: retrieval of subcellular compartments for query set")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(comppidb))
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR::validate_db_df(comppidb, dbtype = "comppidb")

  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F)

  target_compartments <- list()
  target_compartments[['all']] <- data.frame()
  target_compartments[['grouped']] <- data.frame()

  target_compartments_all <- as.data.frame(
    dplyr::inner_join(comppidb, target_genes, by=c("symbol"))
  )

  if(nrow(target_compartments_all) > 0){
    target_compartments_all <- as.data.frame(
      target_compartments_all %>%
      dplyr::select(-c(go_ontology,uniprot_acc)) %>%
      dplyr::left_join(dplyr::select(genedb, symbol, genename),by=c("symbol")) %>%
      dplyr::mutate(genelink = paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",entrezgene,"' target='_blank'>",symbol,"</a>")) %>%
      dplyr::mutate(compartment = paste0('<a href=\'http://amigo.geneontology.org/amigo/term/',go_id,'\' target=\'_blank\'>',go_term,'</a>'))
    )

    target_compartments[['grouped']] <- as.data.frame(
      target_compartments_all %>%
        dplyr::group_by(go_id, go_term, compartment) %>%
        dplyr::summarise(targets = paste(unique(symbol),collapse=", "),
                         targetlinks = paste(unique(genelink),collapse=", "),
                         n = dplyr::n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(go_id,go_term))
    )

    target_compartments[['all']] <- target_compartments_all %>%
      dplyr::select(-c(entrezgene,go_id,go_term,genelink)) %>%
      dplyr::select(symbol, genename, compartment, dplyr::everything())

  }
  return(target_compartments)


}
