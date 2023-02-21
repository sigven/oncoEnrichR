annotate_synleth_paralog_pairs <- function(
  qgenes,
  genedb = NULL,
  slparalogdb = NULL) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(is.data.frame(slparalogdb))
  validate_db_df(slparalogdb, dbtype = "slparalog")
  validate_db_df(genedb, dbtype = "genedb")

  lgr::lgr$info(
             paste0("Annotation of membership in predicted synthetic lethal interactions - De Kegel et al., Cell Systems, 2021"))


  paralog_synleth_interactions <- list()
  paralog_synleth_interactions[['both_in_pair']] <- data.frame()
  paralog_synleth_interactions[['single_pair_member']] <- data.frame()

  targetA_interactions <- as.data.frame(
    data.frame("target" = qgenes, stringsAsFactors = F) |>
      dplyr::inner_join(slparalogdb,
                        by = c("target" = "symbol_A1"), multiple = "all") |>
      dplyr::rename(gene_A = "target") |>
      dplyr::left_join(
        dplyr::select(genedb, c("entrezgene", "genename")),
        by = c("entrezgene_A1" = "entrezgene"), multiple = "all"
      ) |>
      dplyr::rename(genename_A = "genename") |>
      dplyr::left_join(
        dplyr::select(genedb, c("entrezgene", "genename")),
        by = c("entrezgene_A2" = "entrezgene"), multiple = "all"
      ) |>
      dplyr::rename(genename_B = "genename",
                    gene_B = "symbol_A2") |>
      dplyr::select(-c("entrezgene_A1", "entrezgene_A2")) |>
      dplyr::select(c("gene_A",
                    "genename_A",
                    "gene_B",
                    "genename_B"),
                    dplyr::everything()) |>
      dplyr::arrange(dplyr::desc(.data$prediction_score))
  )

  targetB_interactions <- as.data.frame(
    data.frame("target" = qgenes, stringsAsFactors = F) |>
      dplyr::inner_join(slparalogdb,
                        by = c("target" = "symbol_A2"), multiple = "all") |>
      dplyr::rename(gene_B = "target") |>
      dplyr::left_join(
        dplyr::select(genedb, c("entrezgene", "genename")),
        by = c("entrezgene_A2" = "entrezgene"), multiple = "all"
      ) |>
      dplyr::rename(genename_B = "genename") |>
      dplyr::left_join(
        dplyr::select(genedb, c("entrezgene", "genename")),
        by = c("entrezgene_A1" = "entrezgene"), multiple = "all"
      ) |>
      dplyr::rename(genename_A = "genename",
                    gene_A = "symbol_A1") |>
      dplyr::select(-c("entrezgene_A1", "entrezgene_A2")) |>
      dplyr::select(c("gene_A",
                    "genename_A",
                    "gene_B",
                    "genename_B"),
                    dplyr::everything()) |>
      dplyr::arrange(dplyr::desc(.data$prediction_score))
  )

  if (NROW(targetB_interactions) > 0 &
     NROW(targetA_interactions) > 0) {
    paralog_synleth_interactions[['both_in_pair']] <- as.data.frame(
      dplyr::select(targetA_interactions,
                    c("gene_A",
                    "gene_B")) |>
        dplyr::inner_join(targetB_interactions,
                          by = c("gene_A","gene_B"), multiple = "all")
    )

    if (NROW(paralog_synleth_interactions[['both_in_pair']]) > 0) {

      targetA_interactions <- targetA_interactions |>
        dplyr::anti_join(paralog_synleth_interactions[['both_in_pair']],
                         by = c("gene_A", "gene_B"))

      targetB_interactions <- targetB_interactions |>
        dplyr::anti_join(paralog_synleth_interactions[['both_in_pair']],
                         by = c("gene_A", "gene_B")) |>
        dplyr::rename(tmp_symbol = "gene_A",
                      tmp_genename = "genename_A") |>
        dplyr::mutate(gene_A = .data$gene_B,
                      genename_A = .data$genename_B,
                      gene_B = .data$tmp_symbol,
                      genename_B = .data$tmp_genename) |>
        dplyr::select(-c("tmp_symbol", "tmp_genename")) |>
        dplyr::select(c("gene_A",
                      "genename_A",
                      "gene_B",
                      "genename_B"),
                      dplyr::everything())

      paralog_synleth_interactions[['single_pair_member']] <-
        targetA_interactions |>
        dplyr::bind_rows(targetB_interactions) |>
        dplyr::arrange(dplyr::desc(.data$prediction_score))

    }
  } else {

    if (NROW(targetB_interactions) > 0) {
      paralog_synleth_interactions[['single_pair_member']] <-
        targetB_interactions |>
        dplyr::rename(tmp_symbol = "gene_A",
                      tmp_genename = "genename_A") |>
        dplyr::mutate(gene_A = .data$gene_B,
                      genename_A = .data$genename_B,
                      gene_B = .data$tmp_symbol,
                      genename_B = .data$tmp_genename) |>
        dplyr::select(-c("tmp_symbol", "tmp_genename")) |>
        dplyr::select(c("gene_A",
                      "genename_A",
                      "gene_B",
                      "genename_B"),
                      dplyr::everything()) |>
        dplyr::arrange(dplyr::desc(.data$prediction_score))
    }
    if (NROW(targetA_interactions) > 0) {
      paralog_synleth_interactions[['single_pair_member']] <-
        targetA_interactions |>
        dplyr::arrange(dplyr::desc(.data$prediction_score))
    }

  }

  rm(targetA_interactions)
  rm(targetB_interactions)

  lgr::lgr$info(
             paste0("Found n = ", NROW(paralog_synleth_interactions[['both_in_pair']]),
                    " predicted interactions of synthetic lethality for which both members are present in the query set"))
  lgr::lgr$info(
             paste0("Found n = ", NROW(paralog_synleth_interactions[['single_pair_member']]),
                    " predicted interactions of synthetic lethality for whih a single member is present in the query set"))

  return(paralog_synleth_interactions)

}
