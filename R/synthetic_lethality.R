annotate_synleth_paralog_pairs <- function(
  qgenes,
  genedb = NULL,
  slparalogdb = NULL,
  logger = NULL){

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(logger))
  stopifnot(is.data.frame(slparalogdb))
  validate_db_df(slparalogdb, dbtype = "slparalog")
  validate_db_df(genedb, dbtype = "genedb")

  log4r_info(logger,
             paste0("Annotation of membership in predicted synthetic lethal interactions - De Kegel et al., Cell Systems, 2021"))


  paralog_synleth_interactions <- list()
  paralog_synleth_interactions[['both_in_pair']] <- data.frame()
  paralog_synleth_interactions[['single_pair_member']] <- data.frame()

  targetA_interactions <- as.data.frame(
    data.frame("target" = qgenes, stringsAsFactors = F) %>%
      dplyr::inner_join(slparalogdb,
                        by = c("target" = "symbol_A1")) %>%
      dplyr::rename(symbol_A = .data$target) %>%
      dplyr::left_join(
        dplyr::select(genedb, .data$entrezgene, .data$genename),
        by = c("entrezgene_A1" = "entrezgene")
      ) %>%
      dplyr::rename(genename_A = .data$genename) %>%
      dplyr::left_join(
        dplyr::select(genedb, .data$entrezgene, .data$genename),
        by = c("entrezgene_A2" = "entrezgene")
      ) %>%
      dplyr::rename(genename_B = .data$genename, symbol_B = .data$symbol_A2) %>%
      dplyr::select(-c(.data$entrezgene_A1, .data$entrezgene_A2)) %>%
      dplyr::select(.data$symbol_A, .data$genename_A, .data$symbol_B,
                    .data$genename_B, dplyr::everything()) %>%
      dplyr::arrange(dplyr::desc(.data$prediction_score))
  )

  targetB_interactions <- as.data.frame(
    data.frame("target" = qgenes, stringsAsFactors = F) %>%
      dplyr::inner_join(slparalogdb,
                        by = c("target" = "symbol_A2")) %>%
      dplyr::rename(symbol_B = .data$target) %>%
      dplyr::left_join(
        dplyr::select(genedb, .data$entrezgene, .data$genename),
        by = c("entrezgene_A2" = "entrezgene")
      ) %>%
      dplyr::rename(genename_B = .data$genename) %>%
      dplyr::left_join(
        dplyr::select(genedb, .data$entrezgene, .data$genename),
        by = c("entrezgene_A1" = "entrezgene")
      ) %>%
      dplyr::rename(genename_A = .data$genename,
                    symbol_A = .data$symbol_A1) %>%
      dplyr::select(-c(.data$entrezgene_A1, .data$entrezgene_A2)) %>%
      dplyr::select(.data$symbol_A, .data$genename_A, .data$symbol_B,
                    .data$genename_B, dplyr::everything()) %>%
      dplyr::arrange(dplyr::desc(.data$prediction_score))
  )

  if(NROW(targetB_interactions) > 0 &
     NROW(targetA_interactions) > 0){
    paralog_synleth_interactions[['both_in_pair']] <- as.data.frame(
      dplyr::select(targetA_interactions,
                    .data$symbol_A,
                    .data$symbol_B) %>%
        dplyr::inner_join(targetB_interactions,
                          by = c("symbol_A","symbol_B"))
    )

    if(NROW(paralog_synleth_interactions[['both_in_pair']]) > 0){

      targetA_interactions <- targetA_interactions %>%
        dplyr::anti_join(paralog_synleth_interactions[['both_in_pair']],
                         by = c("symbol_A", "symbol_B"))

      targetB_interactions <- targetB_interactions %>%
        dplyr::anti_join(paralog_synleth_interactions[['both_in_pair']],
                         by = c("symbol_A", "symbol_B")) %>%
        dplyr::rename(tmp_symbol = .data$symbol_A,
                      tmp_genename = .data$genename_A) %>%
        dplyr::mutate(symbol_A = .data$symbol_B,
                      genename_A = .data$genename_B,
                      symbol_B = .data$tmp_symbol,
                      genename_B = .data$tmp_genename) %>%
        dplyr::select(-c(.data$tmp_symbol, .data$tmp_genename)) %>%
        dplyr::select(.data$symbol_A,
                      .data$genename_A,
                      .data$symbol_B,
                      .data$genename_B, dplyr::everything())

      paralog_synleth_interactions[['single_pair_member']] <-
        targetA_interactions %>%
        dplyr::bind_rows(targetB_interactions) %>%
        dplyr::arrange(dplyr::desc(.data$prediction_score))

    }
  }else{

    if(NROW(targetB_interactions) > 0){
      paralog_synleth_interactions[['single_pair_member']] <-
        targetB_interactions %>%
        dplyr::rename(tmp_symbol = .data$symbol_A,
                      tmp_genename = .data$genename_A) %>%
        dplyr::mutate(symbol_A = .data$symbol_B,
                      genename_A = .data$genename_B,
                      symbol_B = .data$tmp_symbol,
                      genename_B = .data$tmp_genename) %>%
        dplyr::select(-c(.data$tmp_symbol, .data$tmp_genename)) %>%
        dplyr::select(.data$symbol_A,
                      .data$genename_A,
                      .data$symbol_B,
                      .data$genename_B, dplyr::everything()) %>%
        dplyr::arrange(dplyr::desc(.data$prediction_score))
    }
    if(NROW(targetA_interactions) > 0){
      paralog_synleth_interactions[['single_pair_member']] <-
        targetA_interactions %>%
        dplyr::arrange(dplyr::desc(.data$prediction_score))
    }

  }

  rm(targetA_interactions)
  rm(targetB_interactions)

  log4r_info(logger,
             paste0("Found n = ", NROW(paralog_synleth_interactions[['both_in_pair']]),
                    " predicted interactions of synthetic lethality for which both members are present in the query set"))
  log4r_info(logger,
             paste0("Found n = ", NROW(paralog_synleth_interactions[['single_pair_member']]),
                    " predicted interactions of synthetic lethality for whih a single member is present in the query set"))

  return(paralog_synleth_interactions)

}
