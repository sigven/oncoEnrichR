

target_disease_associations <-
  function(qgenes,
           genedb = NULL,
           min_association_score = 0.3){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
      dplyr::distinct()

  rlogging::message(paste0("Open Targets Platform: annotation of protein targets to disease phenotypes (minimum association score =  ",min_association_score,")"))

  target_assocs <- target_genes %>%
    dplyr::left_join(oncoEnrichR::otdb, by="ensembl_gene_id")

  result <- list()
  result[['target']] <- target_genes
  result[['target_assoc']] <- target_assocs
  result[['assoc_pr_gene']] <- list()
  result[['assoc_pr_gene']][['any']] <- data.frame()
  result[['assoc_pr_gene']][['cancer']] <- data.frame()


  tmp <- target_assocs %>%
    dplyr::filter(cancer_phenotype == T &
                    ot_association_score >= min_association_score)

  if(nrow(tmp) > 0){
    result[['assoc_pr_gene']][['cancer']] <-
      as.data.frame(
        target_assocs %>%
          dplyr::filter((cancer_phenotype == T |
                           stringr::str_detect(tolower(efo_name),"carcinoma|cancer|tumor")) &
                          ot_association_score >= min_association_score) %>%
          dplyr::arrange(desc(ot_association_score)) %>%
          dplyr::group_by(symbol) %>%
          dplyr::summarise(n_cancer_phenotypes = dplyr::n(),
                           ot_cancer_rank = as.numeric(max(ot_association_score)),
                           ot_cancer_links = paste(unique(ot_link), collapse=", "),
                           ot_cancer_diseases = paste(unique(efo_name),collapse=", "))
      )
  }

  tmp2 <- target_assocs %>%
    dplyr::filter(ot_association_score >= min_association_score &
                    is.na(cancer_phenotype))

  if(nrow(tmp2) > 0){
    result[['assoc_pr_gene']][['other']] <-
      as.data.frame(
        target_assocs %>%
          dplyr::filter(ot_association_score >= min_association_score &
                          is.na(cancer_phenotype )) %>%
          dplyr::arrange(desc(ot_association_score)) %>%
          dplyr::group_by(symbol) %>%
          dplyr::summarise(ot_links = paste(unique(ot_link), collapse=", "),
                           ot_diseases = paste(unique(efo_name),collapse=", "))
      )
  }

  rm(tmp)
  rm(tmp2)
  result[['target']] <- result[['target']] %>%
    dplyr::left_join(result[['assoc_pr_gene']][['cancer']],by="symbol") %>%
    dplyr::left_join(result[['assoc_pr_gene']][['other']],by="symbol") %>%
    dplyr::mutate(ot_cancer_rank = dplyr::if_else(is.na(ot_cancer_rank),
                                                  as.numeric(0),
                                                  as.numeric(ot_cancer_rank))) %>%
    dplyr::arrange(desc(ot_cancer_rank), desc(n_cancer_phenotypes)) %>%
    dplyr::select(symbol, genename, ensembl_gene_id, oncogene, tumor_suppressor,
                  cancergene_support, ot_cancer_diseases,
                  ot_cancer_links, ot_cancer_rank, ot_diseases, ot_links,
                  gene_function_description) %>%
    dplyr::rename(cancergene_evidence = cancergene_support,
                  disease_associations = ot_diseases,
                  disease_association_links = ot_links,
                  cancer_associations = ot_cancer_diseases,
                  cancer_association_links = ot_cancer_links) %>%
    dplyr::distinct()

  return(result)
}

target_drug_associations <- function(qgenes,
                                    genedb = NULL){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes,
                             stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
    dplyr::distinct()

  rlogging::message(paste0("Open Targets Platform: annotation of protein targets to targeted drugs"))

  result <- list()
  result[['target']] <- target_genes %>%
    dplyr::select(symbol, genename, targeted_cancer_drugs_lp,
                  targeted_cancer_drugs_ep) %>%
    dplyr::filter(!is.na(targeted_cancer_drugs_lp) |
                    !is.na(targeted_cancer_drugs_ep))

  return(result)

}

