verify_query_genes <- function(qgenes, qsource = "symbol", ignore_unknown = F, genedb = NULL, uniprot_acc = NULL){

  stopifnot(qsource == "symbol" | qsource == "entrezgene" |
            qsource == "uniprot_acc" | qsource == "ensembl_gene_id")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(uniprot_acc))
  oncoEnrichR::validate_db_df(genedb, dbtype = "gene")
  oncoEnrichR::validate_db_df(uniprot_acc, dbtype = "uniprot_acc")

  target_genes <- data.frame('qid' = qgenes, stringsAsFactors = F)
  gdb <- genedb %>% dplyr::select(symbol, entrezgene, ensembl_gene_id) %>%
    dplyr::distinct()

  not_found <- data.frame()
  found <- data.frame()

  result <- list()
  result[['found']] <- found
  result[['not_found']] <- not_found
  result[['success']] <- 1


  if(qsource == 'entrezgene'){
    target_genes <- target_genes %>%
      dplyr::left_join(gdb, by = c("qid" = "entrezgene")) %>%
      dplyr::mutate(entrezgene = qid) %>%
      dplyr::distinct()
  }
  if(qsource == 'symbol'){
    target_genes <- dplyr::left_join(target_genes, gdb, by = c("qid" = "symbol")) %>%
      dplyr::mutate(symbol = qid) %>%
      dplyr::distinct()

  }
  if(qsource == 'uniprot_acc'){
    target_genes <- target_genes %>%
      dplyr::left_join(uniprot_acc, by = c("qid" = "uniprot_acc")) %>%
      dplyr::mutate(uniprot_acc = qid) %>%
      dplyr::distinct()
  }
  if(qsource == 'ensembl_gene_id'){
    target_genes <- target_genes %>%
      dplyr::left_join(gdb, by = c("qid" = "ensembl_gene_id")) %>%
      dplyr::mutate(ensembl_gene_id = qid) %>%
      dplyr::distinct()
  }

  found <- target_genes %>%
    dplyr::filter(!is.na(symbol) & !is.na(entrezgene) & !is.na(ensembl_gene_id))

  not_found <- target_genes %>%
    dplyr::filter(is.na(symbol) | is.na(entrezgene) | is.na(ensembl_gene_id))

  result[['found']] <- found
  result[['not_found']] <- not_found

  if(nrow(not_found) > 0){
    if(ignore_unknown == T){
      rlogging::message(paste0("WARNING: query gene identifiers NOT found: ",paste0(not_found$qid,collapse=", ")))
    }else{
      rlogging::warning(paste0("ERROR: query gene identifiers NOT found: ",paste0(not_found$qid,collapse=", ")))
      result[['success']] <- -1
    }
  }else{
    if(nrow(found) == length(qgenes)){
      rlogging::message(paste0('SUCCESS: Identified all genes (n = ',nrow(found),') in query set'))
    }
    else{
      rlogging::warning(paste0("ERROR: query gene identifiers NOT found: ",paste0(target_genes$qid,collapse=", ")," - wrong query_source (",qsource,")?"),"\n")
        result[['success']] <- -1
    }
  }

  result[['found']]$qid <- NULL
  result[['not_found']]$qid <- NULL

  return(result)

}

validate_db_df <- function(df, dbtype = "genedb"){
  stopifnot(is.data.frame(df))
  if(dbtype == "genedb"){
    for(var in c('symbol','entrezgene','p_oncogene','tsgene','cdriver','tcga_driver','ensembl_gene_id','name',
               'gencode_gene_biotype','ot_tractability_compound','signaling_pw','genename','targeted_drugs')){
      stopifnot(var %in% colnames(df))
    }
  }
  if(dbtype == "corum"){
    for(var in c('complex_id','complex_name','protein_complex_purification_method','complex_comment',
                 'disease_comment','citation','citation_link')){
      stopifnot(var %in% colnames(df))
    }
  }
  if(dbtype == "uniprot_acc"){
    for(var in c('symbol','entrezgene','uniprot_acc')){
      stopifnot(var %in% colnames(df))
    }
  }

  if(dbtype == "comppidb"){
    for(var in c('symbol','entrezgene','uniprot_acc','go_id','go_term','go_ontology','annotation_source','annotation_type')){
      stopifnot(var %in% colnames(df))
    }
  }

  if(dbtype == "ppi_nodes"){
    for(var in c('symbol','entrezgene','genename','name','gencode_gene_biotype',
                 'ot_tractability_compound','signaling_pw','query_node','cdriver',
                 'id','tsgene','p_oncogene')){
      stopifnot(var %in% colnames(df))
    }
  }

  if(dbtype == "projectscoredb"){
    for(var in c('symbol','model_name','loss_of_fitness','model_id','model_id_cmp','synonyms',
                 'model_type','tissue','cancer_type','cancer_type_detail',
                 'sample_site','gender','ethnicity')){
      stopifnot(var %in% colnames(df))
    }
  }

  if(dbtype == "ppi_edges"){
    for(var in c('preferredName_A','preferredName_B','entrezgene_a',
                 'entrezgene_b','p_oncogene_A','p_oncogene_B','tsgene_A','tsgene_B','tcga_driver_A',
                 'tcga_driver_B','cdriver_A','cdriver_B','query_node_A','query_node_B','weight',
                 'from','to','fscore','tscore','score','ascore','pscore','nscore','dscore','escore',
                 'interaction_symbol')){
      stopifnot(var %in% colnames(df))
    }
  }


}
