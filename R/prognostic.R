hpa_prognostic_genes <- function(qgenes,
                                 q_id_type = "symbol",
                                genedb = NULL,
                                oeDB = NULL,
                                logger = NULL){
  stopifnot(!is.null(logger))
  log4r_info(logger,
    paste0("Human Protein Atlas: retrieving prognostic ",
           "associations (gene expression) to cancer"))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(oeDB$hpa))
  validate_db_df(genedb, dbtype = "genedb")
  stopifnot(q_id_type == "symbol" | q_id_type == "entrezgene")
  stopifnot(is.character(qgenes))
  stopifnot(is.data.frame(oeDB$hpa))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(q_id_type == 'entrezgene'){
    query_genes_df <-
      data.frame('entrezgene' = qgenes, stringsAsFactors = F) %>%
      dplyr::inner_join(genedb, by = "entrezgene") %>%
      dplyr::distinct()
  }else{
    query_genes_df <- query_genes_df %>%
      dplyr::inner_join(genedb, by = "symbol") %>%
      dplyr::distinct()
  }

  prognostic_associations <- oeDB$hpa %>%
    dplyr::filter(stringr::str_detect(.data$property,"pathology_prognostics")) %>%
    tidyr::separate(.data$value, c("evidence_direction","is_prognostic","p_value"),
                    sep="\\|") %>%
    dplyr::mutate(p_value = as.numeric(.data$p_value)) %>%
    dplyr::mutate(p_value = dplyr::if_else(.data$p_value == 0, 1e-16, .data$p_value)) %>%
    dplyr::mutate(log10_p_value = round(-log10(.data$p_value), digits = 2)) %>%
    dplyr::mutate(property = stringr::str_replace(.data$property,
                                                  "pathology_prognostics_","")) %>%
    dplyr::inner_join(dplyr::select(genedb,
                                    .data$symbol,
                                    .data$ensembl_gene_id),
                      by=c("ensembl_gene_id")) %>%
    dplyr::mutate(hpa_link = paste0("<a href='https://www.proteinatlas.org/",
                                    .data$ensembl_gene_id,
                                    "-",
                                    .data$symbol,
                                    "/pathology/",
                                    stringr::str_replace_all(.data$property,"_","+"),"'",
                                    " target='_blank'>",
                                    stringr::str_to_title(
                                      stringr::str_replace_all(.data$property,"_"," ")
                                      ),
                                    "</a>")) %>%
    dplyr::mutate(property = dplyr::case_when(
      stringr::str_detect(.data$property,"breast") ~ "Breast",
      stringr::str_detect(.data$property,"colorect") ~ "Colon/Rectum",
      stringr::str_detect(.data$property,"pancreat") ~ "Pancreas",
      stringr::str_detect(.data$property,"lung") ~ "Lung",
      stringr::str_detect(.data$property,"prostate") ~ "Prostate",
      stringr::str_detect(.data$property,"testis") ~ "Testis",
      stringr::str_detect(.data$property,"liver") ~ "Liver",
      stringr::str_detect(.data$property,"stomach") ~ "Stomach",
      stringr::str_detect(.data$property,"renal") ~ "Kidney",
      stringr::str_detect(.data$property,"thyroid") ~ "Thyroid",
      stringr::str_detect(.data$property,"ovarian") ~ "Ovary",
      stringr::str_detect(.data$property,"head_and_neck") ~ "Head/Neck",
      stringr::str_detect(.data$property,"endometr") ~ "Endometrium",
      stringr::str_detect(.data$property,"cervical") ~ "Cervix",
      stringr::str_detect(.data$property,"melanom") ~ "Skin",
      stringr::str_detect(.data$property,"glioma") ~ "CNS/Brain",
      stringr::str_detect(.data$property,"urothelial") ~ "Bladder",
      TRUE ~ as.character("Other"))) %>%
    dplyr::arrange(dplyr::desc(.data$log10_p_value)) %>%
    dplyr::mutate(percentile_rank_all =
                    round(dplyr::percent_rank(.data$log10_p_value) * 100, digits = 1))

  all_prog_assocs <- data.frame()

  for(tissue in unique(prognostic_associations$property)){
    assocs <- dplyr::filter(prognostic_associations,
                            .data$property == tissue) %>%
      dplyr::arrange(dplyr::desc(.data$log10_p_value)) %>%
      dplyr::mutate(percentile_rank_site = round(
        dplyr::percent_rank(.data$log10_p_value) * 100, digits = 1))

    all_prog_assocs <- dplyr::bind_rows(all_prog_assocs, assocs)
  }



  n_query_associations <- 0
  prognostoc_query_associations <- 0

  if(nrow(all_prog_assocs) > 0){
    prognostic_query_associations <- all_prog_assocs %>%
      dplyr::inner_join(dplyr::select(query_genes_df, .data$symbol,
                                      .data$ensembl_gene_id, .data$name),
                        by = c("symbol","ensembl_gene_id")) %>%
      dplyr::rename(primary_site = .data$property,
                    tumor_types = .data$hpa_link) %>%
      dplyr::select(-.data$is_prognostic) %>%
      dplyr::select(.data$symbol,
                    .data$primary_site,
                    .data$evidence_direction,
                    .data$tumor_types,
                    .data$p_value,
                    .data$percentile_rank_site,
                    .data$percentile_rank_all,
                    .data$log10_p_value) %>%
      dplyr::arrange(dplyr::desc(.data$log10_p_value))

    n_query_associations <- NROW(prognostic_query_associations)

  }

  #n_omitted <- nrow(query_genes_df) - nrow(tcga_gene_stats)
  log4r_info(logger, paste0("Found n = ",n_query_associations,
                           " prognostic cancer associations within the target set"))

  return(prognostic_query_associations)

}

km_cshl_survival_genes <-

  function(qgenes,
           qsource = "symbol",
           projectsurvivaldb = NULL,
           logger = NULL,
           ...){

    stopifnot(!is.null(logger))
    dot_args <- list(...)
    log4r_info(logger,
      paste0("Project Survival_KM_CSHL: retrieval of ",
             "genetic determinants of cancer survival - ",
            dot_args$genetic_feature
      ))
    stopifnot(is.character(qgenes))
    stopifnot(!is.null(projectsurvivaldb))
    validate_db_df(projectsurvivaldb,
                   dbtype = "survival_km_cshl")

    target_genes <- data.frame("symbol" = qgenes, stringsAsFactors = F)

    km_survival_targets <- data.frame()

    targets_survival <- as.data.frame(
      projectsurvivaldb %>%
        dplyr::inner_join(target_genes, by = c("symbol"))
    )

    symlevels <- levels(projectsurvivaldb$symbol)

    if(nrow(targets_survival) > 0){

      km_survival_targets <- as.data.frame(
        targets_survival %>%
          dplyr::select(.data$symbol,
                        .data$tcga_cohort,
                        .data$z_score) %>%
        dplyr::mutate(
          symbol = factor(.data$symbol, levels = symlevels))
      )


    }

    return(km_survival_targets)

  }

