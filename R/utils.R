validate_query_genes <- function(qgenes,
                               q_id_type = "symbol",
                               qtype = "target",
                               ignore_id_err = F,
                               genedb = NULL,
                               transcript_xref_db = NULL,
                               logger = NULL){

  stopifnot(q_id_type == "symbol" |
              q_id_type == "entrezgene" |
              q_id_type == "ensembl_mrna" |
              q_id_type == "ensembl_protein" |
              q_id_type == "refseq_protein" |
              q_id_type == "refseq_mrna" |
              q_id_type == "uniprot_acc" |
              q_id_type == "ensembl_gene")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(logger))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(transcript_xref_db, dbtype = "transcript_xref")

  alias2entrez <-
    transcript_xref_db %>%
    dplyr::filter(property == "alias") %>%
    dplyr::rename(alias = value) %>%
    dplyr::select(-property)

  target_genes <- data.frame(
    'qid' = unique(qgenes), stringsAsFactors = F)
  gdb <- genedb %>%
    dplyr::select(.data$symbol,
                  .data$entrezgene,
                  .data$name,
                  .data$ensembl_gene_id) %>%
    dplyr::distinct()

  queryset <- list()
  queryset[['found']] <- data.frame()
  queryset[['not_found']] <- data.frame()
  queryset[['all']] <- data.frame()
  queryset[['match_status']] <- "perfect_go"

  qtype_id <- 'symbol'
  if(q_id_type == 'entrezgene' |
     q_id_type == 'symbol' |
     q_id_type == 'ensembl_gene'){

    qtype_id <- q_id_type
    if(q_id_type == 'ensembl_gene'){
      qtype_id <- 'ensembl_gene_id'
    }

    target_genes <- target_genes %>%
      dplyr::left_join(gdb,
                       by = c("qid" = q_id_type)) %>%
      dplyr::mutate(!!rlang::sym(q_id_type) := .data$qid) %>%
      dplyr::distinct()

  }else{

    qtype_id <- q_id_type
    if(q_id_type == "refseq_protein"){
      qtype_id <- "refseq_peptide"
    }
    if(q_id_type == "ensembl_mrna"){
      qtype_id <- "ensembl_transcript_id"
    }
    if(q_id_type == "ensembl_protein"){
      qtype_id <- "ensembl_protein_id"
    }


    assertable::assert_colnames(
      transcript_xref_db,
      c("entrezgene","property","value"),
      only_colnames = T,
      quiet = T
    )

    gene_xref_map <-
      transcript_xref_db %>%
      dplyr::filter(property == qtype_id) %>%
      dplyr::rename(!!rlang::sym(qtype_id) := value) %>%
      dplyr::select(entrezgene, rlang::sym(qtype_id))

    target_genes <- as.data.frame(
      target_genes %>%

        ## map query to entrezgene
        dplyr::left_join(
          gene_xref_map,
          by = c("qid" = qtype_id)) %>%
        dplyr::mutate(!!rlang::sym(qtype_id) := .data$qid) %>%
        ## append other gene annotations
        dplyr::left_join(gdb, by = c("entrezgene")) %>%
        dplyr::distinct() %>%
        dplyr::group_by(
          .data$symbol,
          .data$entrezgene,
          .data$ensembl_gene_id,
          .data$name) %>%
        dplyr::summarise(
          !!rlang::sym(qtype_id) := paste(
            !!rlang::sym(qtype_id), collapse=","),
          qid = paste(.data$qid, collapse=","),
          .groups = "drop")
    )

  }


  # if(q_id_type == 'entrezgene'){
  #   target_genes <- target_genes %>%
  #     dplyr::left_join(gdb,
  #                      by = c("qid" = "entrezgene")) %>%
  #     dplyr::mutate(entrezgene = .data$qid) %>%
  #     dplyr::distinct()
  # }
  # if(q_id_type == 'symbol'){
  #   target_genes <- target_genes %>%
  #     dplyr::left_join(gdb,
  #                      by = c("qid" = "symbol")) %>%
  #     dplyr::mutate(symbol = .data$qid) %>%
  #     dplyr::distinct()
  #
  # }
  # if(q_id_type == 'uniprot_acc'){
  #   target_genes <- as.data.frame(target_genes %>%
  #     dplyr::left_join(uniprot_xref,
  #                      by = c("qid" = "uniprot_acc")) %>%
  #     dplyr::mutate(uniprot_acc = .data$qid) %>%
  #     dplyr::left_join(gdb, by = c("symbol")) %>%
  #     dplyr::distinct() %>%
  #     dplyr::group_by(.data$symbol, .data$entrezgene,
  #                     .data$ensembl_gene_id,
  #                     .data$name) %>%
  #     dplyr::summarise(uniprot_acc = paste(.data$uniprot_acc, collapse=","),
  #                      qid = paste(.data$qid, collapse=","),
  #                      .groups = "drop")
  #   )
  #
  # }
  # if(q_id_type == 'refseq_mrna'){
  #   target_genes <- as.data.frame(target_genes %>%
  #     dplyr::left_join(refseq_mrna_xref,
  #                      by = c("qid" = "refseq_mrna")) %>%
  #     dplyr::mutate(refseq_mrna = .data$qid) %>%
  #     dplyr::left_join(gdb, by = c("symbol")) %>%
  #     dplyr::distinct() %>%
  #     dplyr::group_by(.data$symbol,
  #                     .data$entrezgene,
  #                     .data$ensembl_gene_id,
  #                     .data$name) %>%
  #     dplyr::summarise(refseq_mrna =
  #                        paste(.data$refseq_mrna, collapse=","),
  #                      qid = paste(.data$qid, collapse=","),
  #                      .groups = "drop")
  #   )
  # }
  # if(q_id_type == 'ensembl_mrna'){
  #   target_genes <- as.data.frame(target_genes %>%
  #     dplyr::left_join(ensembl_mrna_xref,
  #                      by = c("qid" = "ensembl_transcript_id")) %>%
  #     dplyr::mutate(ensembl_transcript_id = .data$qid) %>%
  #     dplyr::left_join(gdb, by = c("symbol")) %>%
  #     dplyr::distinct() %>%
  #     dplyr::group_by(.data$symbol, .data$entrezgene,
  #                     .data$ensembl_gene_id,
  #                     .data$name) %>%
  #     dplyr::summarise(ensembl_transcript_id =
  #                        paste(.data$ensembl_transcript_id, collapse=","),
  #                      qid = paste(.data$qid, collapse=","),
  #                      .groups = "drop")
  #   )
  # }
  #
  # if(q_id_type == 'refseq_protein'){
  #   target_genes <- as.data.frame(
  #     target_genes %>%
  #       dplyr::left_join(refseq_protein_xref,
  #                        by = c("qid" = "refseq_peptide")) %>%
  #       dplyr::mutate(refseq_peptide = .data$qid) %>%
  #       dplyr::left_join(gdb, by = c("symbol")) %>%
  #       dplyr::distinct() %>%
  #       dplyr::group_by(.data$symbol, .data$entrezgene,
  #                       .data$ensembl_gene_id,
  #                       .data$name) %>%
  #       dplyr::summarise(refseq_peptide =
  #                          paste(.data$refseq_peptide, collapse=","),
  #                        qid = paste(.data$qid, collapse=","),
  #                        .groups = "drop")
  #   )
  # }
  #
  # if(q_id_type == 'ensembl_protein'){
  #   target_genes <- as.data.frame(
  #     target_genes %>%
  #       dplyr::left_join(ensembl_protein_xref,
  #                        by = c("qid" = "ensembl_protein_id")) %>%
  #       dplyr::mutate(ensembl_protein_id = .data$qid) %>%
  #       dplyr::left_join(gdb, by = c("symbol")) %>%
  #       dplyr::distinct() %>%
  #       dplyr::group_by(.data$symbol, .data$entrezgene,
  #                       .data$ensembl_gene_id,
  #                       .data$name) %>%
  #       dplyr::summarise(ensembl_protein_id =
  #                          paste(.data$ensembl_protein_id, collapse=","),
  #                        qid = paste(.data$qid, collapse=","),
  #                        .groups = "drop")
  #   )
  # }
  #
  #
  # if(q_id_type == 'ensembl_gene'){
  #   target_genes <- target_genes %>%
  #     dplyr::left_join(gdb, by = c("qid" = "ensembl_gene_id")) %>%
  #     dplyr::mutate(ensembl_gene_id = .data$qid) %>%
  #     dplyr::distinct()
  # }

  queryset[['found']] <- target_genes %>%
    dplyr::filter(!is.na(.data$symbol) &
                    !is.na(.data$entrezgene) &
                    !is.na(.data$ensembl_gene_id)) %>%
    dplyr::mutate(alias = F)

  queryset[['not_found']] <- target_genes %>%
    dplyr::filter(is.na(.data$symbol) |
                    is.na(.data$entrezgene) |
                    is.na(.data$ensembl_gene_id))

  if(nrow(queryset[['not_found']]) > 0){
    if(ignore_id_err == T){
      queryset[['match_status']] <- "imperfect_go"

      if(q_id_type == 'symbol'){
        log4r_info(logger, paste0("WARNING: query gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid,collapse=", ")))
        log4r_info(logger, paste0("Trying to map query identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid,collapse=", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], .data$qid),
            alias2entrez,
            by = c("qid" = "alias"))

        ## Check that alias is not an alias for existing query entries (found)
        ## anti_join against found entries

        if(nrow(query_as_alias) > 0){

          if(nrow(queryset[['found']]) > 0){
            query_as_alias <- query_as_alias %>%
              dplyr::anti_join(queryset[['found']],
                               by = "entrezgene") %>%
              dplyr::distinct()
          }
          if(nrow(query_as_alias) > 0){

            query_as_alias <- query_as_alias %>%
              dplyr::left_join(gdb,
                               by = "entrezgene") %>%
                               #by = c("symbol" = "symbol")) %>%
              dplyr::distinct() %>%
              dplyr::mutate(alias = T)

            log4r_info(logger,
              paste0("Mapped query identifiers as gene aliases ",
                     paste0(queryset[['not_found']]$qid, collapse=", ")," ---> ",
                     paste0(query_as_alias$symbol, collapse=", ")))

            queryset[['found']] <-
              dplyr::bind_rows(queryset[['found']], query_as_alias)
            queryset[['not_found']] <- queryset[['not_found']] %>%
              dplyr::anti_join(query_as_alias, by = "qid")

            if(nrow(queryset[['not_found']]) > 0){
              log4r_warn(logger,
                paste0("Query gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid,collapse=", "),
                       " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
            }else{
              queryset[['match_status']] <- "perfect_go"
            }
          }
        }else{
          log4r_warn(logger,
            paste0("Query gene identifiers NOT found: ",
                   paste0(queryset[['not_found']]$qid,collapse=", "),
                   " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
        }

      }else{
        log4r_warn(logger, paste0("Query gene identifiers NOT found: ",
                                 paste0(queryset[['not_found']]$qid,collapse=", ")))
      }

      ## Indicate that processing should stop when encountering invalid query identifiers
    }else{
      queryset[['match_status']] <- "imperfect_stop"

      if(q_id_type == 'symbol'){
        log4r_warn(logger, paste0("WARNING: query gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid,collapse=", ")))
        log4r_info(logger, paste0("Trying to map query identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid,collapse=", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], .data$qid),
            alias2entrez,
            by = c("qid" = "alias"))

        if(nrow(query_as_alias) > 0){
          query_as_alias <- query_as_alias %>%
            dplyr::left_join(gdb,
                             by = "entrezgene") %>%
                             #by = c("symbol" = "symbol")) %>%
            dplyr::distinct() %>%
          dplyr::mutate(alias = T)


          log4r_info(logger,
            paste0("Mapped query identifiers as gene aliases ",
                   paste0(queryset[['not_found']]$qid,collapse=", ")," ---> ",
                   paste0(query_as_alias$qid,collapse=", ")))

          queryset[['found']] <-
            dplyr::bind_rows(queryset[['found']], query_as_alias)
          queryset[['not_found']] <- queryset[['not_found']] %>%
            dplyr::anti_join(query_as_alias, by = "qid")

          if(nrow(queryset[['not_found']]) > 0){
            log4r_info(logger,
              paste0("ERROR: Query gene identifiers NOT found: ",
                     paste0(queryset[['not_found']]$qid,collapse=", "),
                     " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
          }else{
            queryset[['match_status']] <- "perfect_go"

          }
        }else{
          log4r_info(logger, paste0("ERROR: query gene identifiers NOT found: ",
                         paste0(queryset[['not_found']]$qid, collapse=", "),
                         " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))

        }
      }
      else{
        log4r_info(logger, paste0("ERROR: query gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid, collapse=", "),
                       " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
      }
    }
  }
  if(nrow(queryset[['found']]) == length(qgenes)){
    log4r_info(logger, paste0('SUCCESS: Identified all genes (n = ',
                             nrow(queryset[['found']]),') in ',qtype,' set'))
  }else{
    if(nrow(queryset[['found']]) == 0){
      log4r_info(logger, paste0(
        "ERROR: NO query gene identifiers found: ",
        paste0(target_genes$qid,collapse=", "),
        " - wrong query_id_type (",q_id_type,")?","\n"))
      queryset[['match_status']] <- "imperfect_stop"
    }else{
      log4r_info(logger,
        paste0('Identified n = ',
               nrow(queryset[['found']]),' entries in ',
               'target set (n = ',
               nrow(queryset[['not_found']]),' invalid entries)'))
    }

  }

  if(nrow(queryset[['found']]) > 0){
    queryset[['found']] <- queryset[['found']] %>%
      dplyr::mutate(status = 'found') %>%
      dplyr::mutate(
        status = dplyr::if_else(
          .data$status == "found" & .data$alias == T,
          as.character("found_as_alias"),
          as.character(.data$status))
      ) %>%
      dplyr::arrange(dplyr::desc(.data$status), .data$symbol) %>%
      dplyr::mutate(
        genename = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                          .data$entrezgene,"' target='_blank'>",.data$name,"</a>")
      )
  }
  if(nrow(queryset[['not_found']]) > 0){
    queryset[['not_found']] <- queryset[['not_found']] %>%
      dplyr::mutate(status = 'not_found') %>%
      dplyr::mutate(genename = NA) %>%
      dplyr::arrange(.data$qid)
  }

  if(nrow(queryset[['found']]) > 0 | nrow(queryset[['not_found']]) > 0){

    queryset[['all']] <- as.data.frame(
      queryset[['not_found']] %>%
      dplyr::bind_rows(queryset[['found']]) %>%
      dplyr::rename(query_id = .data$qid) %>%
      dplyr::select(.data$query_id,
                    .data$status,
                    .data$symbol,
                    .data$genename) %>%
      dplyr::rowwise() %>%
        dplyr::mutate(
          symbol = dplyr::if_else(.data$status == "not_found",
                                  as.character(NA),
                                  as.character(.data$symbol)))
    )
  }

  queryset[['found']]$qid <- NULL
  queryset[['not_found']]$qid <- NULL


  return(queryset)

}

#' Function that validates a particular db object (data.frame) in oncoEnrichR
#'
#' @param df data.frame with annotation data for oncoEnrichR
#' @param dbtype type of oncoEnrichR datasource
#'
#' @keywords internal
#'

validate_db_df <- function(df, dbtype = "genedb"){

  val <- assertthat::validate_that(
    is.data.frame(df)
  )
  if(!is.logical(val)){
    message(val)
  }

  dbtypes <- c("genedb",
               "protein_complex",
               "dorothea",
               "ligand_receptor",
               "transcript_xref",
               "comppidb",
               "ppi_nodes",
               "fitness_scores",
               "target_priority_scores",
               "survival_km_cshl",
               "ppi_edges")
  if(!(dbtype %in% dbtypes)){
    stop(
      paste0("dbtype '",dbtype,
             "' not recognized, possible values are: ",
             paste(sort(dbtypes), collapse=", ")))
  }
  if(dbtype == "genedb"){
    cols <- c('symbol',
              'entrezgene',
              'oncogene',
              'tumor_suppressor',
              'cancer_driver',
              'citation_links_oncogene',
              'citation_links_tsgene',
              'citation_links_cdriver',
              'ensembl_gene_id',
              'name',
              'gene_summary',
              'gene_biotype',
              'AB_tractability_category',
              'AB_tractability_support',
              'SM_tractability_category',
              'SM_tractability_support',
              'genename',
              'targeted_cancer_drugs_lp',
              'targeted_cancer_drugs_ep',
              'num_go_terms',
              'cancer_max_rank',
              'unknown_function_rank',
              'has_gene_summary')
  }
  if(dbtype == "protein_complex"){
    cols <- c('complex_id',
              'complex_name',
              'purification_method',
              'complex_comment',
              'disease_comment',
              'sources',
              'confidence',
              'complex_literature',
              'complex_literature_support')
  }
  if(dbtype == "ligand_receptor"){
    cols <- c('interaction_id',
              'interaction_name',
              'annotation',
              'pathway_name',
              'interaction_members',
              'ligand',
              'receptor',
              'agonist',
              'antagonist',
              'co_A_receptor',
              'co_I_receptor',
              'literature_support')
  }
  if(dbtype == "dorothea"){
    cols <- c('regulator',
              'target',
              'interaction_sources',
              'confidence_level',
              'mode_of_regulation',
              'tf_target_literature_support',
              'tf_target_literature')
  }
  if(dbtype == "transcript_xref"){
    cols <- c('entrezgene',
              'property',
              'value')
  }
  if(dbtype == "survival_km_cshl"){
    cols <- c('symbol',
              'tcga_cohort',
              'z_score')
  }

  if(dbtype == "target_priority_scores"){
    cols <- c('symbol',
              'gene_id',
              'priority_score',
              'tumor_type')
  }

  if(dbtype == "fitness_scores"){
    cols <- c('symbol',
              'model_name',
              'loss_of_fitness',
              'model_id',
              'model_type',
              'tissue',
              'entrezgene',
              'sample_site',
              'gene_id_project_score')
  }

  if(dbtype == "comppidb"){
    cols <- c('uniprot_acc',
              'go_id',
              'go_term',
              'confidence',
              'annotation_source',
              'annotation_type')
  }

  if(dbtype == "ppi_nodes"){
    cols <- c('symbol',
              'entrezgene',
              'genename',
              'query_node',
              'cancer_driver',
              'id',
              'tumor_suppressor',
              'oncogene')
  }

  if(dbtype == "ppi_edges"){
    cols <- c('preferredName_A',
              'preferredName_B',
              'entrezgene_a',
              'entrezgene_b',
              'oncogene_A',
              'oncogene_B',
              'tsgene_A',
              'tsgene_B',
              'cdriver_A',
              'cdriver_B',
              'query_node_A',
              'query_node_B',
              'weight',
              'from',
              'to',
              'fscore',
              'tscore',
              'score',
              'ascore',
              'pscore',
              'nscore',
              'dscore',
              'escore',
              'interaction_symbol')

  }
  assertable::assert_colnames(df,
                              colnames = cols,
                              only_colnames = F,
                              quiet = T)


}

add_excel_sheet <- function(
  report = NULL,
  workbook = NULL,
  analysis_output = "disease_association",
  tableStyle = "TableStyleMedium15"){

  invisible(assertthat::assert_that(!is.null(report)))
  invisible(assertthat::assert_that(!is.null(report$data)))
  invisible(assertthat::assert_that(!is.null(workbook)))

  target_df <- data.frame()

  if(analysis_output == "query"){
    if(is.data.frame(report$data$query$target)){
      if(NROW(report$data$query$target) > 0){
        target_df <- report$data$query$target %>%
          dplyr::mutate(
            genename =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$genename)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "unknown_function"){
    if(is.data.frame(report$data$unknown_function$hits_df)){
      if(NROW(report$data$unknown_function$hits_df) > 0){
        target_df <- report$data$unknown_function$hits_df %>%
          dplyr::mutate(
            annotation_source = "GO (MsigDB 7.4)/NCBI Gene/UniProt (2021_02)",
            version = NA) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          ) %>%
          dplyr::mutate(
            gene_summary =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$gene_summary)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "disease_association"){
    if(is.data.frame(report$data$disease$target)){
      if(NROW(report$data$disease$target) > 0){
        target_df <- report$data$disease$target %>%
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) %>%
          dplyr::select(-c(.data$cancer_association_links,
                           .data$disease_association_links,
                           .data$cancergene_evidence)) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          ) %>%
          dplyr::mutate(
            gene_summary =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$gene_summary)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "cancer_hallmark"){
    if(is.data.frame(report$data$cancer_hallmark$target)){
      if(NROW(report$data$cancer_hallmark$target) > 0){
        target_df <- report$data$cancer_hallmark$target %>%
          dplyr::mutate(
            symbol =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
              )
          ) %>%
          dplyr::select(-.data$literature_support)
      }
    }
  }


  if(analysis_output == "drug_known"){
    if(is.data.frame(report$data$drug$target_drugs)){
      if(NROW(report$data$drug$target_drugs)){
        target_df <- report$data$drug$target_drugs %>%
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          ) %>%
          dplyr::mutate(
            targeted_cancer_drugs_lp =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(.data$targeted_cancer_drugs_lp)
                  )
                ),
                " , ",
                ", "
              )
          ) %>%
          dplyr::mutate(
            targeted_cancer_drugs_ep =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(.data$targeted_cancer_drugs_ep)
                  )
                ),
                " , ",
                ", "
              )
          )
      }
    }
  }

  if (analysis_output == "drug_tractability") {
    if (is.data.frame(report$data$drug$tractability_ab) &
        is.data.frame(report$data$drug$tractability_sm)) {
      if (NROW(report$data$drug$tractability_ab) |
          NROW(report$data$drug$tractability_sm)) {

        target_df <- data.frame()

        if (NROW(report$data$drug$tractability_sm)){
          df <- report$data$drug$tractability_sm %>%
            dplyr::mutate(
              annotation_source = report$config$resources$opentargets$name,
              version = report$config$resources$opentargets$version) %>%
            dplyr::rename(tractability_category = .data$SM_tractability_category,
                          tractability_support = .data$SM_tractability_support) %>%
            dplyr::mutate(tractability_drugtype = "Small molecules/compounds") %>%
            dplyr::select(.data$annotation_source, .data$version,
                          .data$tractability_drugtype,
                          dplyr::everything()) %>%
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
            ) %>%
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(.data$tractability_support)
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }

        if (NROW(report$data$drug$tractability_ab)){
          df <- report$data$drug$tractability_ab %>%
            dplyr::mutate(
              annotation_source = report$config$resources$opentargets$name,
              version = report$config$resources$opentargets$version) %>%
            dplyr::rename(tractability_category = .data$AB_tractability_category,
                          tractability_support = .data$AB_tractability_support) %>%
            dplyr::mutate(tractability_drugtype = "Antibody") %>%
            dplyr::select(.data$annotation_source, .data$version,
                          .data$tractability_drugtype,
                          dplyr::everything()) %>%
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
            ) %>%
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(.data$tractability_support)
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }
    }
  }

  if(analysis_output == "ligand_receptor"){
    for(c in c('secreted_signaling','cell_cell_contact',
               'ecm_receptor')){

      if(is.data.frame(report$data$ligand_receptor[[c]])){
        if(NROW(report$data$ligand_receptor[[c]]) > 0){

          df <- report$data$ligand_receptor[[c]] %>%
            dplyr::mutate(
              annotation_source = report$config$resources[['cellchatdb']]$name,
              version = report$config$resources[['cellchatdb']]$version) %>%
            dplyr::select(.data$annotation_source, .data$version,
                          dplyr::everything()) %>%
            dplyr::mutate(
              literature_support =
                stringr::str_trim(
                  textclean::replace_html(.data$literature_support)
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }

    }
  }


  if(analysis_output == "protein_complex"){

    for(c in c('omnipath','humap2')){

      if(is.data.frame(report$data$protein_complex[[c]])){
        if(NROW(report$data$protein_complex[[c]]) > 0){


          df <- report$data$protein_complex[[c]] %>%
            dplyr::mutate(
              annotation_source = report$config$resources[[c]]$name,
              version = report$config$resources[[c]]$version) %>%
            dplyr::select(.data$annotation_source, .data$version,
                          dplyr::everything()) %>%
            dplyr::mutate(
              complex_genes =
                stringr::str_replace_all(
                  stringr::str_squish(
                    stringr::str_trim(
                      textclean::replace_html(.data$complex_genes)
                    )
                  ),
                  " , ",
                  ", "
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }
    }

    if(NROW(target_df) > 0 &
       "literature" %in% colnames(target_df)){
      target_df <- target_df %>%
        dplyr::mutate(
          literature =
            stringr::str_trim(
              textclean::replace_html(.data$literature)
            )
        ) %>%
        dplyr::mutate(
          complex_name =
            stringr::str_trim(
              textclean::replace_html(.data$complex_name)
            )
        )
    }

  }


  if(analysis_output == "prognostic_association"){
    if(is.data.frame(report$data$cancer_prognosis$hpa$assocs)){
      if(NROW(report$data$cancer_prognosis$hpa$assocs) > 0){
        target_df <- report$data$cancer_prognosis$hpa$assocs %>%
          dplyr::mutate(
            annotation_source = report$config$resources$hpa$name,
            version = report$config$resources$hpa$version) %>%
          dplyr::select(.data$annotation_source,
                        .data$version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            tumor_types =
              stringr::str_trim(
                textclean::replace_html(.data$tumor_types)
              )
          )
      }
    }
  }

  if(analysis_output == "survival_association"){
    for(t in c('cna','mut','exp')){
      if(is.data.frame(report$data$cancer_prognosis$km_cshl$assocs[[t]])){
        if(NROW(report$data$cancer_prognosis$km_cshl$assocs[[t]]) > 0){
          df <- report$data$cancer_prognosis$km_cshl$assocs[[t]] %>%
            dplyr::mutate(
              annotation_source = "Genetic determinants of survival in cancer (Smith et al., bioRxiv, 2021)",
              version = "v2") %>%
            dplyr::mutate(feature_type = t) %>%
            dplyr::arrange(.data$feature_type, .data$z_score) %>%
            dplyr::select(.data$annotation_source, .data$version,
                          dplyr::everything())

          target_df <- target_df %>%
            dplyr::bind_rows(df)

        }
      }
    }
  }

  if(analysis_output == "coexpression"){

    ## co-expression
    if(is.data.frame(report$data$tcga$co_expression)){
      if(NROW(report$data$tcga$co_expression) > 0){
        target_df <- report$data$tcga$co_expression %>%
          dplyr::mutate(
            annotation_source = report$config$resources$tcga$name,
            version = report$config$resources$tcga$version) %>%
          dplyr::rename(tcga_cohort = .data$tumor) %>%
          dplyr::select(-.data$entrezgene) %>%
          dplyr::rename(partner_oncogene = .data$oncogene,
                        partner_tumor_suppressor = .data$tumor_suppressor,
                        partner_cancer_driver = .data$cancer_driver) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything())
      }
    }
  }

  if(analysis_output == "regulatory"){

    ## regulatory interactions
    for(c in c('pancancer','global')){

      if(is.data.frame(report$data$regulatory$interactions[[c]])){
        if(NROW(report$data$regulatory$interactions[[c]]) > 0){
          df <- report$data$regulatory$interactions[[c]] %>%
            dplyr::mutate(
              dorothea_collection = c,
              annotation_source = report$config$resources$dorothea$name,
              version = report$config$resources$dorothea$version) %>%
            dplyr::select(.data$annotation_source, .data$version,
                          .data$dorothea_collection,
                          dplyr::everything()) %>%
            dplyr::mutate(
              target_name =
                stringr::str_trim(
                  textclean::replace_html(.data$target_name)
                )
            ) %>%
            dplyr::mutate(
              regulator_name =
                stringr::str_trim(
                  textclean::replace_html(.data$regulator_name)
                )
            ) %>%
            dplyr::mutate(
              literature_support =
                stringr::str_trim(
                  textclean::replace_html(.data$literature_support)
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }
    }

  }

  if(analysis_output == "aberration"){

    ## cna aberrations
    for(t in c('cna_ampl','cna_homdel')){
      if(is.data.frame(report$data$tcga$aberration$table[[t]])){
        if(NROW(report$data$tcga$aberration$table[[t]]) > 0){
          df <-
            report$data$tcga$aberration$table[[t]] %>%
            dplyr::mutate(
              annotation_source = report$config$resources$tcga$name,
              version = report$config$resources$tcga$version) %>%
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$gene)
                )
            ) %>%
            dplyr::select(-.data$gene) %>%
            dplyr::select(.data$annotation_source, .data$version,
                          .data$symbol, dplyr::everything())

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }
    }

  }

  if(analysis_output == "subcellcomp"){
    if(is.data.frame(report$data$subcellcomp$all)){
      if(NROW(report$data$subcellcomp$all) > 0){
        target_df <- report$data$subcellcomp$all %>%
          dplyr::mutate(
            annotation_source = report$config$resources$comppi$name,
            version = report$config$resources$comppi$version) %>%
          dplyr::select(-c(.data$ggcompartment)) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            compartment =
              stringr::str_trim(
                textclean::replace_html(.data$compartment)
              )
          ) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          )
      }
    }
  }

  if(analysis_output == "enrichment"){
    enrichment_df <- data.frame()
    for(e in c('go','wikipathway','kegg','msigdb','netpath')){
      if(NROW(report$data$enrichment[[e]]) > 0){
        enrichment_df <- enrichment_df %>%
          dplyr::bind_rows(
            report$data$enrichment[[e]] %>%
              dplyr::mutate(
                annotation_source = report$config$resources[[e]]$name,
                version = report$config$resources[[e]]$version) %>%
              dplyr::rename(category = .data$db,
                            entrezgene = .data$gene_id)
          ) %>%
          dplyr::select(-c(.data$gene_symbol_link,
                           .data$description_link)) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        .data$category, .data$description,
                        dplyr::everything())
      }
    }
    target_df <- enrichment_df
  }



  if(analysis_output == "crispr_ps_fitness"){
    if(is.data.frame(report$data$crispr_ps$fitness_scores$targets)){
      if(NROW(report$data$crispr_ps$fitness_scores$targets) > 0){
        target_df <- report$data$crispr_ps$fitness_scores$targets %>%
          dplyr::mutate(
            annotation_source = report$config$resources$projectscore$name,
            version = report$config$resources$projectscore$version) %>%
          dplyr::select(-c(.data$symbol_link_ps, .data$cmp_link)) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything())
      }
    }
  }

  if(analysis_output == "crispr_ps_prioritized"){
    if(is.data.frame(report$data$crispr_ps$target_priority_scores$targets)){
      if(NROW(report$data$crispr_ps$target_priority_scores$targets) > 0){
        target_df <- report$data$crispr_ps$target_priority_scores$targets %>%
          dplyr::mutate(
            annotation_source = report$config$resources$projectscore$name,
            version = report$config$resources$projectscore$version) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything())
      }
    }
  }

  if(analysis_output == "ppi"){

    target_df <- data.frame()
    if(is.data.frame(report$data$ppi$complete_network$edges)){
      if(NROW(report$data$ppi$complete_network$edges) > 0){
        target_df <- report$data$ppi$complete_network$edges %>%
          dplyr::rename(symbol_A = .data$preferredName_A,
                        symbol_B = .data$preferredName_B,
                        is_target_A = .data$query_node_A,
                        is_target_B = .data$query_node_B,
                        string_score = .data$score,
                        string_nscore = .data$nscore,
                        string_fscore = .data$fscore,
                        string_pscore = .data$pscore,
                        string_ascore = .data$ascore,
                        string_escore = .data$escore,
                        string_dscore = .data$dscore,
                        string_tscore = .data$tscore
                        ) %>%
          dplyr::select(.data$symbol_A, .data$symbol_B,
                        .data$is_target_A, .data$is_target_B,
                        .data$string_score,
                        .data$string_nscore,
                        .data$string_fscore,
                        .data$string_pscore,
                        .data$string_ascore,
                        .data$string_escore,
                        .data$string_dscore,
                        .data$string_tscore) %>%
          dplyr::mutate(
            annotation_source = report$config$resources$string$name,
            version = report$config$resources$string$version) %>%
          dplyr::select(.data$annotation_source, .data$version,
                        dplyr::everything())
      }
    }

  }

  if(analysis_output == "cell_tissue"){

    target_df <- data.frame()
    for(e in c("tissue_enrichment","scRNA_enrichment")){
      if(is.data.frame(report$data$cell_tissue[[e]]$per_gene)){
        if(NROW(report$data$cell_tissue[[e]]$per_gene) > 0){
          if(e == "tissue_enrichment"){
            df <-
              report$data$cell_tissue[[e]]$per_gene %>%
              dplyr::mutate(
                annotation_source = report$config$resources$gtex$name,
                version = report$config$resources$gtex$version,
                category = stringr::str_replace(e,"_enrichment","")) %>%
              dplyr::select(.data$annotation_source,
                            .data$version,
                            .data$category,
                            dplyr::everything()) %>%
              dplyr::rename(tissue_or_celltype = .data$tissue)

          }else{
            df <-
              report$data$cell_tissue[[e]]$per_gene %>%
              dplyr::mutate(
                annotation_source = report$config$resources$hpa$name,
                version = report$config$resources$hpa$version,
                category = stringr::str_replace(e,"_enrichment","")) %>%
              dplyr::select(.data$annotation_source, .data$version,
                            .data$category, dplyr::everything()) %>%
              dplyr::rename(tissue_or_celltype = .data$cell_type)
          }
          target_df <- target_df %>%
            dplyr::bind_rows(df) %>%
            dplyr::mutate(
              genename =
                stringr::str_trim(
                  textclean::replace_html(.data$genename)
                )
            )
        }
      }
    }
  }

  if(nrow(target_df) > 0){

    analysis_output = stringr::str_replace(
      analysis_output,"_ps_","_")

    openxlsx::addWorksheet(workbook,
                           sheetName = toupper(analysis_output))

    ## set automatic column widths
    openxlsx::setColWidths(workbook,
                           sheet = toupper(analysis_output),
                           cols = 1:ncol(target_df),
                           widths = "auto")

    ## write with default Excel Table style
    openxlsx::writeDataTable(workbook,
                             sheet = toupper(analysis_output),
                             x = target_df,
                             startRow = 1,
                             startCol = 1,
                             colNames = TRUE,
                             tableStyle = tableStyle)
  }

  return(workbook)
}


log4r_layout <-
  function(level, ...) {
    paste0(format(Sys.time()), " - ",
           level, " - ", ..., "\n",
           collapse = "")
  }

log4r_info <- function(log4r_logger, msg) {
  log4r::info(log4r_logger, msg)
}

log4r_debug <- function(log4r_logger, msg) {
  log4r::debug(log4r_logger, msg)
}

log4r_warn <- function(log4r_logger, msg) {
  log4r::warn(log4r_logger, msg)
}

log4r_err <- function(log4r_logger, msg) {
  log4r::error(log4r_logger, msg)
}

file_is_writable <- function(path) {
  assertthat::is.string(path) &&
    file.exists(path) &&
    assertthat::is.writeable(path)
}


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' RCurl url.exists
#'
#' @name url_exists
#' @rdname url_exists
#' @keywords internal
#' @importFrom RCurl url.exists
NULL


#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

#' GSEAbase imports
#'
#' @name gseabase
#' @keywords internal
#' @importFrom GSEABase GeneSet SymbolIdentifier ENSEMBLIdentifier EntrezIdentifier
NULL

utils::globalVariables(c("."))
