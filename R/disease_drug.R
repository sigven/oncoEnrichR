

target_disease_associations <-
  function(qgenes,
           genedb = NULL,
           otdb_all = NULL,
           otdb_gene_rank = NULL,
           show_top_diseases_only = FALSE,
           min_association_score = 0.05){

    lgr::lgr$appenders$console$set_layout(
      lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(otdb_all))
    stopifnot(!is.null(otdb_gene_rank))
    stopifnot(!is.null(genedb))
    stopifnot(is.numeric(min_association_score))
    stopifnot(min_association_score <= 1 && min_association_score > 0)

    validate_db_df(genedb, dbtype = "genedb")
    validate_db_df(otdb_all, dbtype = "opentarget_disease_assoc")
    validate_db_df(otdb_gene_rank, dbtype = "opentarget_disease_site_rank")

    target_genes <- data.frame('symbol' = qgenes, stringsAsFactors = F) |>
      dplyr::inner_join(genedb, by = "symbol") |>
      dplyr::distinct()

    lgr::lgr$info( paste0(
      "Open Targets Platform: annotation of protein targets to disease phenotypes (minimum association score =  ",
      min_association_score,")"))

    target_assocs <- target_genes |>
      dplyr::left_join(otdb_all, by=c("ensembl_gene_id"))

    result <- list()
    result[['target']] <- target_genes
    result[['assoc_pr_gene']] <- list()
    result[['assoc_pr_gene']][['other']] <- data.frame()
    result[['assoc_pr_gene']][['cancer']] <- data.frame()
    result[['ttype_matrix']] <- matrix()

    num_disease_terms <- 1000
    if(show_top_diseases_only){
      num_disease_terms <- 20
    }

    tmp <- list()
    tmp[['cancer_assocs']] <- target_assocs |>
      dplyr::filter(
        !is.na(.data$primary_site) &
          .data$ot_association_score >= min_association_score
      ) |>
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
                    .data$ot_datatype_support,
                    .data$gene_summary) |>
      dplyr::arrange(.data$symbol,
                     .data$ot_association_score)

    if(nrow(tmp[['cancer_assocs']]) > 0){

      tmp[['cancer_assocs']] <- tmp[['cancer_assocs']] |>
        dplyr::mutate(
          ot_link = paste0(
            "<a href='https://platform.opentargets.org/evidence/",
            .data$ensembl_gene_id,"/",
            stringr::str_replace(.data$disease_efo_id,":","_"),
            "' target=\"_blank\">",
            stringr::str_to_title(.data$efo_name),"</a> (", .data$ot_datatype_support,")"))


      otdb_global_cancer_rank <- otdb_gene_rank |>
        dplyr::select(.data$ensembl_gene_id,
                      .data$global_assoc_rank) |>
        dplyr::distinct()

      gene_targetset_cancer_rank <- as.data.frame(
        tmp[['cancer_assocs']] |>
          dplyr::left_join(
            otdb_global_cancer_rank, by = "ensembl_gene_id"
          ) |>
          dplyr::select(.data$symbol, .data$global_assoc_rank) |>
          dplyr::distinct() |>
          dplyr::arrange(dplyr::desc(.data$global_assoc_rank)) |>

          ## within targetset rank
          dplyr::mutate(
            targetset_cancer_rank = round(
              dplyr::percent_rank(.data$global_assoc_rank),
              digits = 5)
          )
      )

      result[['assoc_pr_gene']][['cancer']] <- as.data.frame(
        tmp[['cancer_assocs']] |>
          dplyr::left_join(gene_targetset_cancer_rank,
                           by = "symbol") |>
          dplyr::arrange(dplyr::desc(.data$targetset_cancer_rank),
                         dplyr::desc(.data$ot_association_score)) |>
          dplyr::select(-.data$primary_site) |>
          dplyr::distinct() |>
          dplyr::group_by(.data$symbol) |>
          dplyr::summarise(
            targetset_cancer_rank =
              as.numeric(mean(.data$targetset_cancer_rank)),
            global_cancer_rank =
              as.numeric(mean(.data$global_assoc_rank)),
            ot_cancer_links =
              paste(utils::head(.data$ot_link,
                                num_disease_terms), collapse=", "),
            ot_cancer_diseases =
              paste(utils::head(stringr::str_to_title(.data$efo_name),
                                num_disease_terms), collapse=", "),
            .groups = "drop"
          ) |>
          dplyr::arrange(
            dplyr::desc(.data$targetset_cancer_rank))

      )
    }

    tmp[['disease_assocs']] <- target_assocs |>
      dplyr::filter(
        .data$ot_association_score >= min_association_score &
          is.na(.data$cancer_phenotype)) |>
      ## ignore "genetic disorder"
      dplyr::filter(.data$disease_efo_id != "EFO:0000508")

    if(nrow(tmp[['disease_assocs']]) > 0){

      tmp[['disease_assocs']] <- tmp[['disease_assocs']] |>
        dplyr::mutate(
          ot_link = paste0(
            "<a href='https://platform.opentargets.org/evidence/",
            .data$ensembl_gene_id,"/",
            stringr::str_replace(.data$disease_efo_id,":","_"),
            "' target=\"_blank\">",
            stringr::str_to_title(.data$efo_name),"</a> (", .data$ot_datatype_support,")"))

      gene_targetset_disease_rank <- as.data.frame(
        tmp[['disease_assocs']] |>
          dplyr::group_by(.data$symbol) |>
          dplyr::summarise(
            total_disease_score = sum(.data$ot_association_score),
            .groups = "drop") |>
          dplyr::ungroup() |>
          dplyr::mutate(
            targetset_disease_rank = round(
              dplyr::percent_rank(.data$total_disease_score),
              digits = 2)
          ) |>
          dplyr::select(.data$symbol,
                        .data$targetset_disease_rank)
      )


      result[['assoc_pr_gene']][['other']] <- as.data.frame(
        tmp[['disease_assocs']] |>
          dplyr::left_join(gene_targetset_disease_rank,
                           by = "symbol") |>
          dplyr::arrange(dplyr::desc(.data$targetset_disease_rank),
                         dplyr::desc(.data$ot_association_score)) |>

          dplyr::group_by(.data$symbol) |>
          dplyr::summarise(
            targetset_disease_rank =
              as.numeric(mean(.data$targetset_disease_rank)),
            ot_links =
              paste(utils::head(.data$ot_link, num_disease_terms),
                    collapse=", "),
            ot_diseases =
              paste(utils::head(stringr::str_to_title(.data$efo_name),
                                num_disease_terms), collapse=", "),
            .groups = "drop"
          )
      )
    }

    rm(tmp)

    if(nrow(result[['assoc_pr_gene']][['cancer']]) > 0){
      result[['target']] <- result[['target']] |>
        dplyr::left_join(
          result[['assoc_pr_gene']][['cancer']],by="symbol"
        )
    }else{
      result[['target']]$targetset_cancer_rank <- NA
      result[['target']]$global_cancer_rank <- NA
      result[['target']]$ot_cancer_links <- NA
      result[['target']]$ot_cancer_diseases <- NA
    }

    if(nrow(result[['assoc_pr_gene']][['other']]) > 0){
      result[['target']] <- result[['target']] |>
        dplyr::left_join(
          result[['assoc_pr_gene']][['other']],by="symbol"
        )
    }else{
      result[['target']]$targetset_disease_rank <- NA
      result[['target']]$ot_links <- NA
      result[['target']]$ot_diseases <- NA
    }

    result[['target']] <- result[['target']] |>
      dplyr::mutate(
        targetset_cancer_rank = dplyr::if_else(
          is.na(.data$targetset_cancer_rank),
          as.numeric(0),
          as.numeric(.data$targetset_cancer_rank)
        ),
        targetset_disease_rank = dplyr::if_else(
          is.na(.data$targetset_disease_rank),
          as.numeric(0),
          as.numeric(.data$targetset_disease_rank)
        ),
        global_cancer_rank = dplyr::if_else(
          is.na(.data$global_cancer_rank),
          as.numeric(0),
          as.numeric(.data$global_cancer_rank)
        )) |>
      dplyr::arrange(dplyr::desc(.data$targetset_cancer_rank)) |>
      dplyr::select(.data$symbol,
                    .data$genename,
                    .data$ensembl_gene_id,
                    .data$oncogene,
                    .data$tumor_suppressor,
                    .data$cancer_driver,
                    .data$targetset_cancer_rank,
                    .data$global_cancer_rank,
                    .data$cancergene_support,
                    .data$ot_cancer_diseases,
                    .data$ot_cancer_links,
                    .data$targetset_disease_rank,
                    .data$ot_diseases,
                    .data$ot_links,
                    .data$gene_summary) |>
      dplyr::rename(cancergene_evidence = .data$cancergene_support,
                    disease_associations = .data$ot_diseases,
                    disease_association_links = .data$ot_links,
                    cancer_associations = .data$ot_cancer_diseases,
                    cancer_association_links = .data$ot_cancer_links) |>
      dplyr::distinct()


    ttype_rank_df <- dplyr::select(target_genes,
                                   .data$symbol,
                                   .data$ensembl_gene_id) |>
      dplyr::left_join(
        dplyr::select(otdb_gene_rank,
                      .data$primary_site,
                      .data$ensembl_gene_id,
                      .data$tissue_assoc_rank),
        by = "ensembl_gene_id"
      ) |>
      dplyr::select(-.data$ensembl_gene_id)

    if(nrow(ttype_rank_df) > 0){

      result[['ttype_matrix']] <- as.data.frame(
        ttype_rank_df |>
          dplyr::filter(!is.na(.data$primary_site) &
                          !is.na(.data$tissue_assoc_rank))
      )

      if(nrow(result[['ttype_matrix']]) > 0){

        result[['ttype_matrix']] <- as.data.frame(
          result[['ttype_matrix']] |>
            dplyr::filter(!(.data$primary_site == "Adrenal Gland" |
                              .data$primary_site == "Eye" |
                              .data$primary_site == "Vulva/Vagina" |
                              .data$primary_site == "Ampulla of Vater" |
                              .data$primary_site == "Peritoneum")) |>
            dplyr::mutate(tissue_assoc_rank =
                            .data$tissue_assoc_rank * 100) |>
            dplyr::mutate(
              symbol = factor(
                .data$symbol, levels = unique(result$target$symbol)
              )
            ) |>
            dplyr::arrange(.data$symbol) |>
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
                                    genedb = NULL){

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")

  target_genes <- data.frame('symbol' = qgenes,
                             stringsAsFactors = F) |>
    dplyr::inner_join(genedb, by = "symbol") |>
    dplyr::distinct()

  lgr::lgr$info( paste0("Open Targets Platform: annotation of protein targets to targeted drugs (cancer indications only)"))

  result <- list()
  result[['target_drugs']] <- data.frame()
  result[['tractability_ab']] <- data.frame()
  result[['tractability_sm']] <- data.frame()

  result[['target_drugs']] <- target_genes |>
    dplyr::select(.data$symbol, .data$genename,
                  .data$targeted_cancer_drugs_lp,
                  .data$targeted_cancer_drugs_ep) |>
    dplyr::filter(!is.na(.data$targeted_cancer_drugs_lp) |
                    !is.na(.data$targeted_cancer_drugs_ep))

  lgr::lgr$info( paste0("Open Targets Platform: annotation of target tractabilities (druggability)"))

  result[['tractability_ab']] <- target_genes |>
    dplyr::select(.data$ensembl_gene_id,
                  .data$symbol,
                  .data$AB_tractability_category,
                  .data$AB_tractability_support) |>
    dplyr::mutate(symbol = paste0(
      "<a href='https://platform.opentargets.org/target/",
      .data$ensembl_gene_id,"' target='_blank'>", .data$symbol, "</a>"
    )) |>
    dplyr::select(-.data$ensembl_gene_id) |>
    dplyr::arrange(.data$AB_tractability_category,
                   dplyr::desc(nchar(.data$AB_tractability_support)))

  result[['tractability_sm']] <- target_genes |>
    dplyr::select(.data$ensembl_gene_id,
                  .data$symbol,
                  .data$SM_tractability_category,
                  .data$SM_tractability_support) |>
    dplyr::mutate(symbol = paste0(
      "<a href='https://platform.opentargets.org/target/",
      .data$ensembl_gene_id,"' target='_blank'>", .data$symbol, "</a>"
    )) |>
    dplyr::select(-.data$ensembl_gene_id) |>
    dplyr::arrange(.data$SM_tractability_category,
                   dplyr::desc(nchar(.data$SM_tractability_support)))

  return(result)

}

