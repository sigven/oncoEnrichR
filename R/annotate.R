

target_disease_associations <-
  function(qgenes,
           genedb = NULL,
           show_top_diseases_only = FALSE,
           min_association_score = 0.2){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
      dplyr::distinct()

  oncoEnrichR:::log4r_info(paste0("Open Targets Platform: annotation of protein targets to disease phenotypes (minimum association score =  ",min_association_score,")"))

  target_assocs <- target_genes %>%
    dplyr::left_join(oncoEnrichR::otdb$all, by=c("ensembl_gene_id","symbol"))

  result <- list()
  result[['target']] <- target_genes
  result[['target_assoc']] <- target_assocs
  result[['assoc_pr_gene']] <- list()
  result[['assoc_pr_gene']][['any']] <- data.frame()
  result[['assoc_pr_gene']][['cancer']] <- data.frame()
  result[['ttype_matrix']] <- matrix()

  num_top_disease_terms <- 20

  tmp <- target_assocs %>%
    dplyr::filter(
      !is.na(primary_site) &
        ot_association_score >= min_association_score
    ) %>%
    dplyr::select(symbol, genename,
                  ensembl_gene_id,
                  oncogene,
                  tumor_suppressor,
                  cancergene_support,
                  ot_association_score,
                  disease_efo_id,
                  efo_name,
                  primary_site,
                  ot_link,
                  gene_summary)

  if(nrow(tmp) > 0){

    gene_targetset_cancer_rank <- as.data.frame(
      tmp %>%
        dplyr::group_by(symbol, primary_site) %>%
        dplyr::summarise(
          mean_tissue_score = mean(ot_association_score),
          .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(
          target_cancer_score = sum(mean_tissue_score),
          .groups = "drop") %>%
        dplyr::arrange(desc(target_cancer_score)) %>%
        dplyr::mutate(
          targetset_cancer_prank = round(
            dplyr::percent_rank(target_cancer_score),
            digits = 2)
        ) %>%
        dplyr::select(symbol, targetset_cancer_prank)
    )

    result[['assoc_pr_gene']][['cancer']] <- tmp %>%
      dplyr::left_join(gene_targetset_cancer_rank,
                        by = "symbol") %>%
      dplyr::arrange(desc(targetset_cancer_prank))

    if(show_top_diseases_only){
      result[['assoc_pr_gene']][['cancer']] <- as.data.frame(
        result[['assoc_pr_gene']][['cancer']] %>%
          dplyr::group_by(symbol) %>%
          dplyr::summarise(
            #n_cancer_phenotypes = dplyr::n(),
            # ot_cancer_rank =
            #   as.numeric(max(ot_association_score)),
            targetset_cancer_prank =
              as.numeric(mean(targetset_cancer_prank)),
            ot_cancer_links =
              paste(head(ot_link, num_top_disease_terms), collapse=", "),
            ot_cancer_diseases =
              paste(head(stringr::str_to_title(efo_name),
                         num_top_disease_terms),collapse=", "),
            .groups = "drop"
          )
      )
    }else{

      result[['assoc_pr_gene']][['cancer']] <- as.data.frame(
        result[['assoc_pr_gene']][['cancer']] %>%
          dplyr::group_by(symbol) %>%
          dplyr::summarise(
            #n_cancer_phenotypes = dplyr::n(),
            # ot_cancer_rank =
            #   as.numeric(max(ot_association_score)),
            # ot_cancer_rank =
            #   as.numeric(mean(ot_association_score)),
            targetset_cancer_prank =
              as.numeric(mean(targetset_cancer_prank)),
            ot_cancer_links =
              paste(unique(ot_link), collapse=", "),
            ot_cancer_diseases =
              paste(unique(stringr::str_to_title(efo_name)),
                           collapse=", "),
            .groups = "drop"
          )
      )
    }
  }

  tmp2 <- target_assocs %>%
    dplyr::filter(ot_association_score >= min_association_score &
                    is.na(cancer_phenotype))

  if(nrow(tmp2) > 0){

    gene_targetset_disease_rank <- as.data.frame(
      tmp2 %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(
          total_disease_score = sum(ot_association_score),
          .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          targetset_disease_prank = round(
            dplyr::percent_rank(total_disease_score),
            digits = 2)
        ) %>%
        dplyr::select(symbol, targetset_disease_prank)
    )


    result[['assoc_pr_gene']][['other']] <- tmp2 %>%
      dplyr::left_join(gene_targetset_disease_rank,
                        by = "symbol") %>%
      dplyr::arrange(desc(targetset_disease_prank))

    if(show_top_diseases_only){
      result[['assoc_pr_gene']][['other']] <- as.data.frame(
        result[['assoc_pr_gene']][['other']] %>%
          dplyr::group_by(symbol) %>%
          dplyr::summarise(
            targetset_disease_prank =
              as.numeric(mean(targetset_disease_prank)),
            ot_links =
              paste(head(ot_link, num_top_disease_terms), collapse=", "),
            ot_diseases =
              paste(head(stringr::str_to_title(efo_name),
                         num_top_disease_terms), collapse=", "),
            .groups = "drop"
          )
      )
    }else{
      result[['assoc_pr_gene']][['other']] <- as.data.frame(
          result[['assoc_pr_gene']][['other']] %>%
            dplyr::group_by(symbol) %>%
            dplyr::summarise(
              targetset_disease_prank =
                as.numeric(mean(targetset_disease_prank)),
              ot_links =
                paste(unique(ot_link), collapse=", "),
              ot_diseases =
                paste(unique(stringr::str_to_title(efo_name)),
                      collapse=", "),
              .groups = "drop"
            )
        )
    }
  }

  rm(tmp)
  rm(tmp2)
  result[['target']] <- result[['target']] %>%
    dplyr::left_join(result[['assoc_pr_gene']][['cancer']],by="symbol") %>%
    dplyr::left_join(result[['assoc_pr_gene']][['other']],by="symbol") %>%
    dplyr::mutate(
      targetset_cancer_prank =
        dplyr::if_else(is.na(targetset_cancer_prank),
                       as.numeric(0),
                       as.numeric(targetset_cancer_prank))) %>%
    dplyr::arrange(desc(targetset_cancer_prank)) %>%
    dplyr::select(symbol, genename,
                  ensembl_gene_id,
                  oncogene, tumor_suppressor,
                  cancergene_support,
                  ot_cancer_diseases,
                  ot_cancer_links,
                  targetset_cancer_prank,
                  targetset_disease_prank,
                  ot_diseases,
                  ot_links,
                  gene_summary) %>%
    dplyr::rename(cancergene_evidence = cancergene_support,
                  disease_associations = ot_diseases,
                  disease_association_links = ot_links,
                  cancer_associations = ot_cancer_diseases,
                  cancer_association_links = ot_cancer_links) %>%
    dplyr::distinct()


  ttype_rank_df <- dplyr::select(target_genes, symbol) %>%
    dplyr::left_join(
      dplyr::select(oncoEnrichR::otdb$site_rank, primary_site,
                    symbol, tissue_assoc_rank),
      by = "symbol"
    )

  if(nrow(ttype_rank_df) > 0){

    result[['ttype_matrix']] <- as.data.frame(
      ttype_rank_df %>%
        dplyr::filter(!is.na(primary_site) & !is.na(tissue_assoc_rank))
    )

    if(nrow(result[['ttype_matrix']]) > 0){

      result[['ttype_matrix']] <- as.data.frame(
        result[['ttype_matrix']] %>%
        dplyr::filter(!(primary_site == "Adrenal Gland" |
                          primary_site == "Eye" |
                          primary_site == "Vulva/Vagina" |
                          primary_site == "Ampulla of Vater" |
                          primary_site == "Peritoneum")) %>%
        dplyr::mutate(tissue_assoc_rank =
                        tissue_assoc_rank * 100) %>%
        dplyr::mutate(
          symbol = factor(
            symbol, levels = result$target$symbol
          )
        ) %>%
        dplyr::arrange(symbol) %>%
        tidyr::pivot_wider(names_from = primary_site,
                           values_from = tissue_assoc_rank)
      )
    }

    rownames(result[['ttype_matrix']]) <-
      result[['ttype_matrix']]$symbol
    result[['ttype_matrix']]$symbol <- NULL
    result[['ttype_matrix']] <- as.matrix(result[['ttype_matrix']])
  }



  return(result)
}

target_drug_associations <- function(qgenes,
                                     cancerdrugdb = NULL,
                                    genedb = NULL){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(cancerdrugdb))
  stopifnot(!is.null(cancerdrugdb[['tractability']]))
  validate_db_df(genedb, dbtype = "genedb")
  #validate_db_df(genedb, dbtype = "drug_tractability")

  target_genes <- data.frame('symbol' = qgenes,
                             stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
    dplyr::distinct()

  oncoEnrichR:::log4r_info(paste0("Open Targets Platform: annotation of protein targets to targeted drugs"))

  result <- list()
  result[['target_drugs']] <- data.frame()
  result[['tractability_ab']] <- data.frame()
  result[['tractability_sm']] <- data.frame()

  result[['target_drugs']] <- target_genes %>%
    dplyr::select(symbol, genename,
                  targeted_cancer_drugs_lp,
                  targeted_cancer_drugs_ep) %>%
    dplyr::filter(!is.na(targeted_cancer_drugs_lp) |
                    !is.na(targeted_cancer_drugs_ep))

  oncoEnrichR:::log4r_info(paste0("Open Targets Platform: annotation of target tractabilities (druggability)"))

  result[['tractability_ab']] <- target_genes %>%
    dplyr::select(ensembl_gene_id) %>%
    dplyr::inner_join(cancerdrugdb[['tractability']][['ab']],
                      by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::arrange(AB_tractability_category)

  result[['tractability_sm']] <- target_genes %>%
    dplyr::select(ensembl_gene_id) %>%
    dplyr::inner_join(cancerdrugdb[['tractability']][['sm']],
                      by = "ensembl_gene_id") %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::arrange(SM_tractability_category)

  return(result)

}

