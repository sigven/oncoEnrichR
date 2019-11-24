
annotate_protein_complex <- function(qgenes, genedb = NULL, corum_db = NULL, uniprot_acc = NULL){

  rlogging::message("CORUM: retrieval of protein complex information for query set")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(corum_db))
  stopifnot(!is.null(uniprot_acc))
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR::validate_db_df(corum_db, dbtype = "corum")
  oncoEnrichR::validate_db_df(uniprot_acc, dbtype = "uniprot_acc")


  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F) %>%
    dplyr::left_join(uniprot_acc, by = "symbol") %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(uniprot_acc))

  results <- data.frame()
  if(nrow(target_genes) > 0){
    targets_uniprot_acc <- paste(target_genes$uniprot_acc,collapse=",")
    query_api <- paste0('http://omnipathdb.org/complexes?proteins=',targets_uniprot_acc,'&databases=CORUM')
    results <- read.table(query_api,header=T,sep="\t",quote="",stringsAsFactors = F,na.strings=c("","NA"))
    #results <- as.data.frame(readr::read_tsv(query_api))

    qtarget <- target_genes %>%
      dplyr::select(entrezgene,symbol) %>%
      dplyr::rename(target_gene = symbol)

    if(nrow(results) > 0){
      results <- as.data.frame(results %>%
        tidyr::separate_rows(identifiers,sep=";") %>%
        dplyr::filter(stringr::str_detect(identifiers,"CORUM")) %>%
        dplyr::mutate(complex_id = stringr::str_replace(identifiers,"CORUM:","")) %>%
        dplyr::select(components, complex_id) %>%
        dplyr::rename(uniprot_acc = components) %>%
        dplyr::mutate(complex_id = as.numeric(complex_id)) %>%
        dplyr::left_join(corum_db,by=c("complex_id")) %>%
        dplyr::filter(!is.na(complex_name)) %>%
        dplyr::distinct() %>%
        tidyr::separate_rows(uniprot_acc,sep="_") %>%
        dplyr::left_join(dplyr::select(uniprot_acc, entrezgene, uniprot_acc), by=c("uniprot_acc" = "uniprot_acc")) %>%
        dplyr::left_join(qtarget,by=c("entrezgene")) %>%
        dplyr::left_join(dplyr::filter(dplyr::select(genedb, entrezgene, symbol),!is.na(entrezgene)),by=c("entrezgene")) %>%
        dplyr::mutate(genelink = paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",entrezgene,"' target='_blank'>",symbol,"</a>")) %>%
        dplyr::group_by(complex_name, disease_comment, complex_comment, protein_complex_purification_method, citation_link) %>%
        dplyr::summarise(target_genes = paste(unique(target_gene),collapse=","), complex_genes = paste(unique(genelink),collapse=",")) %>%
        dplyr::mutate(target_genes = stringr::str_replace_all(target_genes,",NA$|^NA,","")) %>%
        dplyr::mutate(target_genes = stringr::str_replace_all(target_genes,",NA,",",")) %>%
        dplyr::rename(purification_method = protein_complex_purification_method, citation = citation_link) %>%
        dplyr::select(complex_name, target_genes, citation, complex_genes, disease_comment, complex_comment, purification_method) %>%
        dplyr::mutate(num_target_members = stringr::str_count(target_genes,",")) %>%
        dplyr::arrange(desc(num_target_members)) %>%
        dplyr::select(-num_target_members)
      )

    }
  }
  return(results)

}
