validate_query_genes <- function(qgenes,
                               q_id_type = "symbol",
                               qtype = "target",
                               ignore_id_err = F,
                               genedb = NULL,
                               ensembl_mrna_xref = NULL,
                               refseq_mrna_xref = NULL,
                               refseq_protein_xref = NULL,
                               ensembl_protein_xref = NULL,
                               uniprot_xref = NULL){

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
  stopifnot(!is.null(uniprot_xref))
  stopifnot(!is.null(refseq_mrna_xref))
  stopifnot(!is.null(ensembl_mrna_xref))
  stopifnot(!is.null(refseq_protein_xref))
  stopifnot(!is.null(ensembl_protein_xref))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR:::validate_db_df(uniprot_xref, dbtype = "uniprot_xref")
  oncoEnrichR:::validate_db_df(refseq_mrna_xref, dbtype = "refseq_mrna_xref")
  oncoEnrichR:::validate_db_df(ensembl_mrna_xref, dbtype = "ensembl_mrna_xref")
  oncoEnrichR:::validate_db_df(refseq_protein_xref, dbtype = "refseq_protein_xref")
  oncoEnrichR:::validate_db_df(ensembl_protein_xref, dbtype = "ensembl_protein_xref")


  target_genes <- data.frame('qid' = unique(qgenes), stringsAsFactors = F)
  gdb <- genedb %>%
    dplyr::select(symbol, entrezgene,
                  name, ensembl_gene_id) %>%
    dplyr::distinct()

  queryset <- list()
  queryset[['found']] <- data.frame()
  queryset[['not_found']] <- data.frame()
  queryset[['all']] <- data.frame()
  queryset[['match_status']] <- "perfect_go"

  if(q_id_type == 'entrezgene'){
    target_genes <- target_genes %>%
      dplyr::left_join(gdb,
                       by = c("qid" = "entrezgene")) %>%
      dplyr::mutate(entrezgene = qid) %>%
      dplyr::distinct()
  }
  if(q_id_type == 'symbol'){
    target_genes <- target_genes %>%
      dplyr::left_join(gdb,
                       by = c("qid" = "symbol")) %>%
      dplyr::mutate(symbol = qid) %>%
      dplyr::distinct()

  }
  if(q_id_type == 'uniprot_acc'){
    target_genes <- as.data.frame(target_genes %>%
      dplyr::left_join(uniprot_xref,
                       by = c("qid" = "uniprot_acc")) %>%
      dplyr::mutate(uniprot_acc = qid) %>%
      dplyr::left_join(gdb, by = c("symbol")) %>%
      dplyr::distinct() %>%
      dplyr::group_by(symbol, entrezgene,
                      ensembl_gene_id,
                      name) %>%
      dplyr::summarise(uniprot_acc = paste(uniprot_acc, collapse=","),
                       qid = paste(qid, collapse=","),
                       .groups = "drop")
    )

  }
  if(q_id_type == 'refseq_mrna'){
    target_genes <- as.data.frame(target_genes %>%
      dplyr::left_join(refseq_mrna_xref,
                       by = c("qid" = "refseq_mrna")) %>%
      dplyr::mutate(refseq_mrna = qid) %>%
      dplyr::left_join(gdb, by = c("symbol")) %>%
      dplyr::distinct() %>%
      dplyr::group_by(symbol, entrezgene,
                      ensembl_gene_id,
                      name) %>%
      dplyr::summarise(refseq_mrna =
                         paste(refseq_mrna, collapse=","),
                       qid = paste(qid, collapse=","),
                       .groups = "drop")
    )
  }
  if(q_id_type == 'ensembl_mrna'){
    target_genes <- as.data.frame(target_genes %>%
      dplyr::left_join(ensembl_mrna_xref,
                       by = c("qid" = "ensembl_transcript_id")) %>%
      dplyr::mutate(ensembl_transcript_id = qid) %>%
      dplyr::left_join(gdb, by = c("symbol")) %>%
      dplyr::distinct() %>%
      dplyr::group_by(symbol, entrezgene,
                      ensembl_gene_id,
                      name) %>%
      dplyr::summarise(ensembl_transcript_id =
                         paste(ensembl_transcript_id, collapse=","),
                       qid = paste(qid, collapse=","),
                       .groups = "drop")
    )
  }

  if(q_id_type == 'refseq_protein'){
    target_genes <- as.data.frame(
      target_genes %>%
        dplyr::left_join(refseq_protein_xref,
                         by = c("qid" = "refseq_peptide")) %>%
        dplyr::mutate(refseq_peptide = qid) %>%
        dplyr::left_join(gdb, by = c("symbol")) %>%
        dplyr::distinct() %>%
        dplyr::group_by(symbol, entrezgene,
                        ensembl_gene_id,
                        name) %>%
        dplyr::summarise(refseq_peptide =
                           paste(refseq_peptide, collapse=","),
                         qid = paste(qid, collapse=","),
                         .groups = "drop")
    )
  }

  if(q_id_type == 'ensembl_protein'){
    target_genes <- as.data.frame(
      target_genes %>%
        dplyr::left_join(ensembl_protein_xref,
                         by = c("qid" = "ensembl_protein_id")) %>%
        dplyr::mutate(ensembl_protein_id = qid) %>%
        dplyr::left_join(gdb, by = c("symbol")) %>%
        dplyr::distinct() %>%
        dplyr::group_by(symbol, entrezgene,
                        ensembl_gene_id,
                        name) %>%
        dplyr::summarise(ensembl_protein_id =
                           paste(ensembl_protein_id, collapse=","),
                         qid = paste(qid, collapse=","),
                         .groups = "drop")
    )
  }


  if(q_id_type == 'ensembl_gene'){
    target_genes <- target_genes %>%
      dplyr::left_join(gdb, by = c("qid" = "ensembl_gene_id")) %>%
      dplyr::mutate(ensembl_gene_id = qid) %>%
      dplyr::distinct()
  }

  queryset[['found']] <- target_genes %>%
    dplyr::filter(!is.na(symbol) &
                    !is.na(entrezgene) &
                    !is.na(ensembl_gene_id)) %>%
    dplyr::mutate(alias = F)

  queryset[['not_found']] <- target_genes %>%
    dplyr::filter(is.na(symbol) |
                    is.na(entrezgene) |
                    is.na(ensembl_gene_id))

  if(nrow(queryset[['not_found']]) > 0){
    if(ignore_id_err == T){
      queryset[['match_status']] <- "imperfect_go"

      if(q_id_type == 'symbol'){
        rlogging::message(paste0("WARNING: query gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid,collapse=", ")))
        rlogging::message(paste0("Trying to map query identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid,collapse=", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], qid),
            oncoEnrichR::genedb[['alias2primary']],
            by = c("qid" = "alias"))

        ## Check that alias is not an alias for existing query entries (found)
        ## anti_join against found entries

        if(nrow(query_as_alias) > 0){

          if(nrow(queryset[['found']]) > 0){
            query_as_alias <- query_as_alias %>%
              dplyr::anti_join(queryset[['found']],
                               by = "symbol") %>%
              dplyr::distinct()
          }
          if(nrow(query_as_alias) > 0){

            query_as_alias <- query_as_alias %>%
              dplyr::left_join(gdb,
                               by = c("symbol" = "symbol")) %>%
              dplyr::distinct() %>%
              dplyr::mutate(alias = T)

            rlogging::message(
              paste0("Mapped query identifiers as gene aliases ",
                     paste0(queryset[['not_found']]$qid,collapse=", ")," ---> ",
                     paste0(query_as_alias$qid,collapse=", ")))

            queryset[['found']] <-
              dplyr::bind_rows(queryset[['found']], query_as_alias)
            queryset[['not_found']] <- queryset[['not_found']] %>%
              dplyr::anti_join(query_as_alias, by = "qid")

            if(nrow(queryset[['not_found']]) > 0){
              rlogging::message(
                paste0("WARNING: query gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid,collapse=", "),
                       " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
            }else{
              queryset[['match_status']] <- "perfect_go"
            }
          }
        }else{
          rlogging::message(
            paste0("WARNING: query gene identifiers NOT found: ",
                   paste0(queryset[['not_found']]$qid,collapse=", "),
                   " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
        }

      }else{
        rlogging::message(paste0("WARNING: query gene identifiers NOT found: ",
                                 paste0(queryset[['not_found']]$qid,collapse=", ")))
      }

      ## Indicate that processing should stop when encountering invalid query identifiers
    }else{
      queryset[['match_status']] <- "imperfect_stop"

      if(q_id_type == 'symbol'){
        rlogging::message(paste0("WARNING: query gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid,collapse=", ")))
        rlogging::message(paste0("Trying to map query identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid,collapse=", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], qid),
            oncoEnrichR::genedb[['alias2primary']],
            by = c("qid" = "alias"))

        if(nrow(query_as_alias) > 0){
          query_as_alias <- query_as_alias %>%
            dplyr::left_join(gdb,
                             by = c("symbol" = "symbol")) %>%
            dplyr::distinct() %>%
          dplyr::mutate(alias = T)


          rlogging::message(
            paste0("Mapped query identifiers as gene aliases ",
                   paste0(queryset[['not_found']]$qid,collapse=", ")," ---> ",
                   paste0(query_as_alias$qid,collapse=", ")))

          queryset[['found']] <-
            dplyr::bind_rows(queryset[['found']], query_as_alias)
          queryset[['not_found']] <- queryset[['not_found']] %>%
            dplyr::anti_join(query_as_alias, by = "qid")

          if(nrow(queryset[['not_found']]) > 0){
            rlogging::message(
              paste0("ERROR: query gene identifiers NOT found: ",
                     paste0(queryset[['not_found']]$qid,collapse=", "),
                     " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
          }else{
            queryset[['match_status']] <- "perfect_go"

          }
        }else{
          rlogging::message(paste0("ERROR: query gene identifiers NOT found: ",
                         paste0(queryset[['not_found']]$qid, collapse=", "),
                         " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))

        }
      }
      else{
        rlogging::message(paste0("ERROR: query gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid, collapse=", "),
                       " (make sure that primary identifiers/symbols are used, not aliases or synonyms)"))
      }
    }
  }
  if(nrow(queryset[['found']]) == length(qgenes)){
    rlogging::message(paste0('SUCCESS: Identified all genes (n = ',
                             nrow(queryset[['found']]),') in ',qtype,' set'))
  }else{
    if(nrow(queryset[['found']]) == 0){
      rlogging::message(paste0(
        "ERROR: NO query gene identifiers found: ",
        paste0(target_genes$qid,collapse=", "),
        " - wrong query_id_type (",q_id_type,")?"),"\n")
      queryset[['match_status']] <- "imperfect_stop"
    }else{
      rlogging::message(
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
          status == "found" & alias == T,
          as.character("found_as_alias"),
          as.character(status))
      ) %>%
      dplyr::arrange(desc(status), symbol) %>%
      dplyr::mutate(
        genename = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                          entrezgene,"' target='_blank'>",name,"</a>")
      )
  }
  if(nrow(queryset[['not_found']]) > 0){
    queryset[['not_found']] <- queryset[['not_found']] %>%
      dplyr::mutate(status = 'not_found') %>%
      dplyr::mutate(genename = NA) %>%
      dplyr::arrange(qid)
  }

  if(nrow(queryset[['found']]) > 0 | nrow(queryset[['not_found']]) > 0){

    queryset[['all']] <- as.data.frame(
      queryset[['not_found']] %>%
      dplyr::bind_rows(queryset[['found']]) %>%
      dplyr::rename(query_id = qid) %>%
      dplyr::select(query_id, status, symbol,
                    genename) %>%
      dplyr::rowwise() %>%
        dplyr::mutate(
          symbol = dplyr::if_else(status == "not_found",
                                  as.character(NA),
                                  as.character(symbol)))
    )
  }

  queryset[['found']]$qid <- NULL
  queryset[['not_found']]$qid <- NULL


  return(queryset)

}

validate_db_df <- function(df, dbtype = "genedb"){

  val <- assertthat::validate_that(
    is.data.frame(df)
  )
  if(!is.logical(val)){
    message(val)
  }

  dbtypes <- c("genedb",
               "corum",
               "uniprot_xref",
               "refseq_mrna_xref",
               "ensembl_mrna_xref",
               "refseq_protein_xref",
               "ensembl_protein_xref",
               "comppidb",
               "ppi_nodes",
               "fitness_scores",
               "target_priority_scores",
               "ppi_edges",
               "pdf")
  if(!(dbtype %in% dbtypes)){
    rlogging::stop(
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
              'ensembl_gene_id',
              'name',
              'gene_summary',
              'gencode_gene_biotype',
              'ot_tractability_compound',
              'genename',
              'targeted_cancer_drugs_lp',
              'targeted_cancer_drugs_ep')
  }
  if(dbtype == "corum"){
    cols <- c('complex_id',
              'complex_name',
              'protein_complex_purification_method',
              'complex_comment',
              'disease_comment',
              'citation','citation_link')
  }
  ## poorly defined genes (pdf)
  if(dbtype == "pdf"){
    cols <- c('symbol',
              'genename',
              'num_go_terms',
              'unknown_function_rank',
              'gene_summary',
              'has_gene_summary')
  }
  if(dbtype == "uniprot_xref"){
    cols <- c('symbol','uniprot_acc')
  }
  if(dbtype == "refseq_mrna_xref"){
    cols <- c('symbol','refseq_mrna')
  }
  if(dbtype == "ensembl_mrna_xref"){
    cols <- c('symbol','ensembl_transcript_id')
  }
  if(dbtype == "ensembl_protein_xref"){
    cols <- c('symbol','ensembl_protein_id')
  }
  if(dbtype == "refseq_protein_xref"){
    cols <- c('symbol','refseq_peptide')
  }

  if(dbtype == "target_priority_scores"){
    cols <- c('symbol',
              'gene_id',
              'priority_score',
              'tumor_type')
  }

  if(dbtype == "fitness_scores"){
    cols <- c('symbol','model_name',
              'loss_of_fitness','model_id',
              'model_type','tissue',
              'entrezgene',
              'sample_site',
              'gene_id_project_score')
  }

  if(dbtype == "comppidb"){
    cols <- c('symbol',
              'go_id',
              'go_term',
              'annotation_source',
              'annotation_type')
  }

  if(dbtype == "ppi_nodes"){
    cols <- c('symbol',
              'entrezgene',
              'genename',
              'name',
              'gencode_gene_biotype',
              'ot_tractability_compound',
              'query_node',
              'cancer_driver',
              'id',
              'tumor_suppressor','oncogene')
  }

  if(dbtype == "ppi_edges"){
    cols <- c('preferredName_A','preferredName_B',
              'entrezgene_a','entrezgene_b',
              'oncogene_A','oncogene_B',
              'tsgene_A','tsgene_B',
              'cdriver_A','cdriver_B',
              'query_node_A','query_node_B',
              'weight','from','to','fscore',
              'tscore','score','ascore',
              'pscore','nscore','dscore','escore',
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
  analysis_output = "disease",
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
                  textclean::replace_html(genename)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "disease"){
    if(is.data.frame(report$data$disease$target)){
      if(NROW(report$data$disease$target) > 0){
        target_df <- report$data$disease$target %>%
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) %>%
          dplyr::select(-c(cancer_association_links,
                           disease_association_links,
                           cancergene_evidence)) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(genename)
              )
          ) %>%
          dplyr::mutate(
            gene_summary =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(gene_summary)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "protein_complex"){
    if(is.data.frame(report$data$protein_complex$complex)){
      if(NROW(report$data$protein_complex$complex) > 0){
        target_df <- report$data$protein_complex$complex %>%
          dplyr::mutate(
            annotation_source = report$config$resources$corum$name,
            version = report$config$resources$corum$version) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            complex_genes =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(complex_genes)
                  )
                ),
                " , ",
                ", "
              )
          ) %>%
          dplyr::mutate(
            citation =
              stringr::str_trim(
                textclean::replace_html(citation)
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
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(genename)
              )
          ) %>%
          dplyr::mutate(
            gene_summary =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(gene_summary)
                )
              )
          )
      }
    }
  }

  if(analysis_output == "cancer_prognosis"){
    if(is.data.frame(report$data$cancer_prognosis$assocs)){
      if(NROW(report$data$cancer_prognosis$assocs) > 0){
        target_df <- report$data$cancer_prognosis$assocs %>%
          dplyr::mutate(
            annotation_source = report$config$resources$hpa$name,
            version = report$config$resources$hpa$version) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            tumor_types =
              stringr::str_trim(
                textclean::replace_html(tumor_types)
              )
          )
      }
    }
  }

  if(analysis_output == "tcga_coexpression"){

    ## co-expression
    if(is.data.frame(report$data$tcga$co_expression)){
      if(NROW(report$data$tcga$co_expression) > 0){
        target_df <- report$data$tcga$co_expression %>%
          dplyr::mutate(
            annotation_source = report$config$resources$tcga$name,
            version = report$config$resources$tcga$version) %>%
          dplyr::rename(tcga_cohort = tumor) %>%
          dplyr::select(-entrezgene) %>%
          dplyr::rename(partner_oncogene = oncogene,
                        partner_tumor_suppressor = tumor_suppressor,
                        partner_cancer_driver = cancer_driver) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything())
      }
    }
  }

  if(analysis_output == "tcga_aberration"){

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
                  textclean::replace_html(gene)
                )
            ) %>%
            dplyr::select(-gene) %>%
            dplyr::select(annotation_source, version,
                          symbol, dplyr::everything())

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
          dplyr::select(-c(ggcompartment)) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            compartment =
              stringr::str_trim(
                textclean::replace_html(compartment)
              )
          ) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(genename)
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
              dplyr::rename(category = db,
                            entrezgene = gene_id)
          ) %>%
          dplyr::select(-c(gene_symbol_link,
                           description_link)) %>%
          dplyr::select(annotation_source, version,
                        category, description,
                        dplyr::everything())
      }
    }
    target_df <- enrichment_df
  }

  if(analysis_output == "drug"){
    if(is.data.frame(report$data$drug$target_drugs)){
      if(NROW(report$data$drug$target_drugs)){
        target_df <- report$data$drug$target_drugs %>%
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) %>%
          dplyr::select(annotation_source, version,
                        dplyr::everything()) %>%
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(genename)
              )
          ) %>%
          dplyr::mutate(
            targeted_cancer_drugs_lp =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(targeted_cancer_drugs_lp)
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
                    textclean::replace_html(targeted_cancer_drugs_ep)
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
            dplyr::rename(tractability_category = SM_tractability_category,
                          tractability_support = SM_tractability_support) %>%
            dplyr::mutate(tractability_drugtype = "Small molecules/compounds") %>%
            dplyr::select(annotation_source, version,
                          tractability_drugtype,
                          dplyr::everything()) %>%
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(symbol)
                )
            ) %>%
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(tractability_support)
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
            dplyr::rename(tractability_category = AB_tractability_category,
                          tractability_support = AB_tractability_support) %>%
            dplyr::mutate(tractability_drugtype = "Antibody") %>%
            dplyr::select(annotation_source, version,
                          tractability_drugtype,
                          dplyr::everything()) %>%
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(symbol)
                )
            ) %>%
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(tractability_support)
                )
            )

          target_df <- target_df %>%
            dplyr::bind_rows(df)
        }
      }
    }
  }

  if(analysis_output == "crispr_ps_fitness"){
    if(is.data.frame(report$data$crispr_ps$fitness_scores$targets)){
      if(NROW(report$data$crispr_ps$fitness_scores$targets) > 0){
        target_df <- report$data$crispr_ps$fitness_scores$targets %>%
          dplyr::mutate(
            annotation_source = report$config$resources$projectscore$name,
            version = report$config$resources$projectscore$version) %>%
          dplyr::select(-c(symbol_link_ps, cmp_link)) %>%
          dplyr::select(annotation_source, version,
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
          dplyr::select(annotation_source, version,
                        dplyr::everything())
      }
    }
  }

  if(analysis_output == "ppi"){

    target_df <- data.frame()
    if(is.data.frame(report$data$ppi$complete_network$edges)){
      if(NROW(report$data$ppi$complete_network$edges) > 0){
        target_df <- report$data$ppi$complete_network$edges %>%
          dplyr::rename(symbol_A = preferredName_A,
                        symbol_B = preferredName_B,
                        is_target_A = query_node_A,
                        is_target_B = query_node_B,
                        string_score = score,
                        string_nscore = nscore,
                        string_fscore = fscore,
                        string_pscore = pscore,
                        string_ascore = ascore,
                        string_escore = escore,
                        string_dscore = dscore,
                        string_tscore = tscore
                        ) %>%
          dplyr::select(symbol_A, symbol_B,
                        is_target_A, is_target_B,
                        string_score,
                        string_nscore,
                        string_fscore,
                        string_pscore,
                        string_ascore,
                        string_escore,
                        string_dscore,
                        string_tscore) %>%
          dplyr::mutate(
            annotation_source = report$config$resources$string$name,
            version = report$config$resources$string$version) %>%
          dplyr::select(annotation_source, version,
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
              dplyr::select(annotation_source, version,
                            category, dplyr::everything()) %>%
              dplyr::rename(tissue_or_celltype = tissue)

          }else{
            df <-
              report$data$cell_tissue[[e]]$per_gene %>%
              dplyr::mutate(
                annotation_source = report$config$resources$hpa$name,
                version = report$config$resources$hpa$version,
                category = stringr::str_replace(e,"_enrichment","")) %>%
              dplyr::select(annotation_source, version,
                            category, dplyr::everything()) %>%
              dplyr::rename(tissue_or_celltype = cell_type)
          }
          target_df <- target_df %>%
            dplyr::bind_rows(df) %>%
            dplyr::mutate(
              genename =
                stringr::str_trim(
                  textclean::replace_html(genename)
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

