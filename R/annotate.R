

target_disease_associations <-
  function(qgenes,
           genedb = NULL,
           oeDB = NULL,
           show_top_diseases_only = FALSE,
           min_association_score = 0.1,
           logger = NULL){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(oeDB))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(logger))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
      dplyr::distinct()

  log4r_info(logger, paste0(
    "Open Targets Platform: annotation of protein targets to disease phenotypes (minimum association score =  ",
    min_association_score,")"))

  target_assocs <- target_genes %>%
    dplyr::left_join(oeDB$otdb$all, by=c("ensembl_gene_id"))

  result <- list()
  result[['target']] <- target_genes
  result[['target_assoc']] <- target_assocs
  result[['assoc_pr_gene']] <- list()
  result[['assoc_pr_gene']][['other']] <- data.frame()
  result[['assoc_pr_gene']][['cancer']] <- data.frame()
  result[['ttype_matrix']] <- matrix()

  num_top_disease_terms <- 20

  tmp <- target_assocs %>%
    dplyr::filter(
      !is.na(.data$primary_site) &
        .data$ot_association_score >= min_association_score
    ) %>%
    dplyr::select(.data$symbol,
                  .data$genename,
                  .data$ensembl_gene_id,
                  .data$oncogene,
                  .data$tumor_suppressor,
                  .data$cancergene_support,
                  .data$ot_association_score,
                  .data$disease_efo_id,
                  .data$efo_name,
                  .data$primary_site,
                  #.data$ot_link,
                  .data$gene_summary)

  if(nrow(tmp) > 0){

    tmp <- tmp %>%
      dplyr::mutate(
        ot_link = paste0(
          "<a href='https://platform.opentargets.org/evidence/",
          .data$ensembl_gene_id,"/",
          stringr::str_replace(.data$disease_efo_id,":","_"),
          "' target=\"_blank\">",
          stringr::str_to_title(.data$efo_name),"</a>"))

    gene_targetset_cancer_rank <- as.data.frame(
      tmp %>%
        dplyr::group_by(.data$symbol, .data$primary_site) %>%
        dplyr::summarise(
          mean_tissue_score = mean(.data$ot_association_score),
          .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$symbol) %>%
        dplyr::summarise(
          target_cancer_score = sum(.data$mean_tissue_score),
          .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$target_cancer_score)) %>%
        dplyr::mutate(
          targetset_cancer_prank = round(
            dplyr::percent_rank(.data$target_cancer_score),
            digits = 2)
        ) %>%
        dplyr::select(.data$symbol,
                      .data$targetset_cancer_prank)
    )

    result[['assoc_pr_gene']][['cancer']] <- tmp %>%
      dplyr::left_join(gene_targetset_cancer_rank,
                        by = "symbol") %>%
      dplyr::arrange(dplyr::desc(.data$targetset_cancer_prank))

    if(show_top_diseases_only){
      result[['assoc_pr_gene']][['cancer']] <- as.data.frame(
        result[['assoc_pr_gene']][['cancer']] %>%
          dplyr::group_by(.data$symbol) %>%
          dplyr::summarise(
            targetset_cancer_prank =
              as.numeric(mean(.data$targetset_cancer_prank)),
            ot_cancer_links =
              paste(utils::head(.data$ot_link,
                                num_top_disease_terms), collapse=", "),
            ot_cancer_diseases =
              paste(utils::head(stringr::str_to_title(.data$efo_name),
                                num_top_disease_terms), collapse=", "),
            .groups = "drop"
          )
      )
    }else{

      result[['assoc_pr_gene']][['cancer']] <- as.data.frame(
        result[['assoc_pr_gene']][['cancer']] %>%
          dplyr::group_by(.data$symbol) %>%
          dplyr::summarise(
            targetset_cancer_prank =
              as.numeric(mean(.data$targetset_cancer_prank)),
            ot_cancer_links =
              paste(unique(.data$ot_link), collapse=", "),
            ot_cancer_diseases =
              paste(unique(stringr::str_to_title(.data$efo_name)),
                           collapse=", "),
            .groups = "drop"
          )
      )
    }
  }

  tmp2 <- target_assocs %>%
    dplyr::filter(.data$ot_association_score >= min_association_score &
                    is.na(.data$cancer_phenotype))

  if(nrow(tmp2) > 0){

    tmp2 <- tmp2 %>%
      dplyr::mutate(
        ot_link = paste0(
          "<a href='https://platform.opentargets.org/evidence/",
          .data$ensembl_gene_id,"/",
          stringr::str_replace(.data$disease_efo_id,":","_"),
          "' target=\"_blank\">",
          stringr::str_to_title(.data$efo_name),"</a>"))

    gene_targetset_disease_rank <- as.data.frame(
      tmp2 %>%
        dplyr::group_by(.data$symbol) %>%
        dplyr::summarise(
          total_disease_score = sum(.data$ot_association_score),
          .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          targetset_disease_prank = round(
            dplyr::percent_rank(.data$total_disease_score),
            digits = 2)
        ) %>%
        dplyr::select(.data$symbol,
                      .data$targetset_disease_prank)
    )


    result[['assoc_pr_gene']][['other']] <- tmp2 %>%
      dplyr::left_join(gene_targetset_disease_rank,
                        by = "symbol") %>%
      dplyr::arrange(dplyr::desc(.data$targetset_disease_prank))

    if(show_top_diseases_only){
      result[['assoc_pr_gene']][['other']] <- as.data.frame(
        result[['assoc_pr_gene']][['other']] %>%
          dplyr::group_by(.data$symbol) %>%
          dplyr::summarise(
            targetset_disease_prank =
              as.numeric(mean(.data$targetset_disease_prank)),
            ot_links =
              paste(utils::head(.data$ot_link, num_top_disease_terms),
                    collapse=", "),
            ot_diseases =
              paste(utils::head(stringr::str_to_title(.data$efo_name),
                         num_top_disease_terms), collapse=", "),
            .groups = "drop"
          )
      )
    }else{
      result[['assoc_pr_gene']][['other']] <- as.data.frame(
          result[['assoc_pr_gene']][['other']] %>%
            dplyr::group_by(.data$symbol) %>%
            dplyr::summarise(
              targetset_disease_prank =
                as.numeric(mean(.data$targetset_disease_prank)),
              ot_links =
                paste(unique(.data$ot_link), collapse=", "),
              ot_diseases =
                paste(unique(stringr::str_to_title(.data$efo_name)),
                      collapse=", "),
              .groups = "drop"
            )
        )
    }
  }

  rm(tmp)
  rm(tmp2)

  if(nrow(result[['assoc_pr_gene']][['cancer']]) > 0){
    result[['target']] <- result[['target']] %>%
      dplyr::left_join(
        result[['assoc_pr_gene']][['cancer']],by="symbol"
        )
  }else{
    result[['target']]$targetset_cancer_prank <- NA
    result[['target']]$ot_cancer_links <- NA
    result[['target']]$ot_cancer_diseases <- NA
  }

  if(nrow(result[['assoc_pr_gene']][['other']]) > 0){
    result[['target']] <- result[['target']] %>%
      dplyr::left_join(
        result[['assoc_pr_gene']][['other']],by="symbol"
      )
  }else{
    result[['target']]$targetset_disease_prank <- NA
    result[['target']]$ot_links <- NA
    result[['target']]$ot_diseases <- NA
  }

  result[['target']] <- result[['target']] %>%
    dplyr::mutate(
      targetset_cancer_prank =
        dplyr::if_else(
          is.na(.data$targetset_cancer_prank),
          as.numeric(0),
          as.numeric(.data$targetset_cancer_prank)
          )
      ) %>%
    dplyr::arrange(dplyr::desc(.data$targetset_cancer_prank)) %>%
    dplyr::select(.data$symbol,
                  .data$genename,
                  .data$ensembl_gene_id,
                  .data$oncogene,
                  .data$tumor_suppressor,
                  .data$targetset_cancer_prank,
                  .data$cancergene_support,
                  .data$ot_cancer_diseases,
                  .data$ot_cancer_links,
                  .data$targetset_disease_prank,
                  .data$ot_diseases,
                  .data$ot_links,
                  .data$gene_summary) %>%
    dplyr::rename(cancergene_evidence = .data$cancergene_support,
                  disease_associations = .data$ot_diseases,
                  disease_association_links = .data$ot_links,
                  cancer_associations = .data$ot_cancer_diseases,
                  cancer_association_links = .data$ot_cancer_links) %>%
    dplyr::distinct()


  ttype_rank_df <- dplyr::select(target_genes,
                                 .data$symbol,
                                 .data$ensembl_gene_id) %>%
    dplyr::left_join(
      dplyr::select(oeDB$otdb$site_rank,
                    .data$primary_site,
                    .data$ensembl_gene_id,
                    .data$tissue_assoc_rank),
      by = "ensembl_gene_id"
    ) %>%
    dplyr::select(-.data$ensembl_gene_id)

  if(nrow(ttype_rank_df) > 0){

    result[['ttype_matrix']] <- as.data.frame(
      ttype_rank_df %>%
        dplyr::filter(!is.na(.data$primary_site) &
                        !is.na(.data$tissue_assoc_rank))
    )

    if(nrow(result[['ttype_matrix']]) > 0){

      result[['ttype_matrix']] <- as.data.frame(
        result[['ttype_matrix']] %>%
        dplyr::filter(!(.data$primary_site == "Adrenal Gland" |
                          .data$primary_site == "Eye" |
                          .data$primary_site == "Vulva/Vagina" |
                          .data$primary_site == "Ampulla of Vater" |
                          .data$primary_site == "Peritoneum")) %>%
        dplyr::mutate(tissue_assoc_rank =
                        .data$tissue_assoc_rank * 100) %>%
        dplyr::mutate(
          symbol = factor(
            .data$symbol, levels = result$target$symbol
          )
        ) %>%
        dplyr::arrange(.data$symbol) %>%
        tidyr::pivot_wider(names_from = .data$primary_site,
                           values_from = .data$tissue_assoc_rank)
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
                                    genedb = NULL,
                                    logger = NULL){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(logger))
  stopifnot(!is.null(cancerdrugdb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes,
                             stringsAsFactors = F) %>%
    dplyr::inner_join(genedb, by = "symbol") %>%
    dplyr::distinct()

  log4r_info(logger, paste0("Open Targets Platform: annotation of protein targets to targeted drugs (cancer indications only)"))

  result <- list()
  result[['target_drugs']] <- data.frame()
  result[['tractability_ab']] <- data.frame()
  result[['tractability_sm']] <- data.frame()

  result[['target_drugs']] <- target_genes %>%
    dplyr::select(.data$symbol, .data$genename,
                  .data$targeted_cancer_drugs_lp,
                  .data$targeted_cancer_drugs_ep) %>%
    dplyr::filter(!is.na(.data$targeted_cancer_drugs_lp) |
                    !is.na(.data$targeted_cancer_drugs_ep))

  log4r_info(logger, paste0("Open Targets Platform: annotation of target tractabilities (druggability)"))

  result[['tractability_ab']] <- target_genes %>%
    dplyr::select(.data$ensembl_gene_id,
                  .data$symbol,
                  .data$AB_tractability_category,
                  .data$AB_tractability_support) %>%
    dplyr::mutate(symbol = paste0(
      "<a href='https://platform.opentargets.org/target/",
      ensembl_gene_id,"' target='_blank'>", symbol, "</a>"
    )) %>%
    dplyr::select(-.data$ensembl_gene_id) %>%
    dplyr::arrange(.data$AB_tractability_category,
                   dplyr::desc(nchar(.data$AB_tractability_support)))

  result[['tractability_sm']] <- target_genes %>%
    dplyr::select(.data$ensembl_gene_id,
                  .data$symbol,
                  .data$SM_tractability_category,
                  .data$SM_tractability_support) %>%
    dplyr::mutate(symbol = paste0(
      "<a href='https://platform.opentargets.org/target/",
      ensembl_gene_id,"' target='_blank'>", symbol, "</a>"
    )) %>%
    dplyr::select(-.data$ensembl_gene_id) %>%
    dplyr::arrange(.data$SM_tractability_category,
                   dplyr::desc(nchar(.data$SM_tractability_support)))

  return(result)

}

