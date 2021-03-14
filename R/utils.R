validate_query_genes <- function(qgenes,
                               q_id_type = "symbol",
                               qtype = "target",
                               ignore_id_err = F,
                               genedb = NULL,
                               uniprot_acc = NULL){

  stopifnot(q_id_type == "symbol" | q_id_type == "entrezgene" |
            q_id_type == "uniprot_acc" | q_id_type == "ensembl_gene_id")
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(uniprot_acc))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  oncoEnrichR:::validate_db_df(uniprot_acc, dbtype = "uniprot_acc")

  target_genes <- data.frame('qid' = qgenes, stringsAsFactors = F)
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
    target_genes <- target_genes %>%
      dplyr::left_join(uniprot_acc,
                       by = c("qid" = "uniprot_acc")) %>%
      dplyr::mutate(uniprot_acc = qid) %>%
      dplyr::distinct()
  }
  if(q_id_type == 'ensembl_gene_id'){
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
               'query set (n = ',
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

  # invisible(assertthat::assert_that(
  #   is.data.frame(df),
  #   msg = paste0("Argument 'df' must by of ",
  #                "type data.frame, not ", class(df))))
  dbtypes <- c("genedb","corum","uniprot_acc",
               "comppidb","ppi_nodes",
               "projectscoredb",
               "ppi_edges","pdf")
  if(!(dbtype %in% dbtypes)){
    rlogging::stop(
      paste0("dbtype '",dbtype,
             "' not recognized, possible values are: ",
             paste(sort(dbtypes), collapse=", ")))
  }
  if(dbtype == "genedb"){
    cols <- c('symbol','entrezgene','oncogene',
              'tumor_suppressor','cancer_driver',
              'ensembl_gene_id','name',
              'gene_summary',
              'gencode_gene_biotype',
              'ot_tractability_compound',
              'genename','targeted_cancer_drugs_lp',
              'targeted_cancer_drugs_ep')
  }
  if(dbtype == "corum"){
    cols <- c('complex_id','complex_name',
              'protein_complex_purification_method',
              'complex_comment','disease_comment',
              'citation','citation_link')
  }
  ## poorly defined genes (pdf)
  if(dbtype == "pdf"){
    cols <- c('symbol','genename',
              'num_go_terms',
              'unknown_function_rank','gene_summary',
              'has_gene_summary')
  }
  if(dbtype == "uniprot_acc"){
    cols <- c('symbol','entrezgene','uniprot_acc')
  }

  if(dbtype == "comppidb"){
    cols <- c('symbol','entrezgene','uniprot_acc',
              'go_id','go_term',
              'go_ontology','annotation_source',
              'annotation_type')
  }

  if(dbtype == "ppi_nodes"){
    cols <- c('symbol','entrezgene','genename',
              'name','gencode_gene_biotype',
              'ot_tractability_compound',
              'query_node','cancer_driver',
              'id','tumor_suppressor','oncogene')
  }

  if(dbtype == "projectscoredb"){
    cols <- c('symbol','model_name',
              'loss_of_fitness','model_id',
              'model_link_cmp','synonyms',
              'model_type','tissue',
              'cancer_type','cancer_type_detail',
              'sample_site','gender',
              'ethnicity','symbol_link_ps',
              'gene_id_project_score')
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
    if(is.data.frame(report$data$disease$target$target)){
      if(NROW(report$data$disease$target$target) > 0){
        target_df <- report$data$disease$target$target %>%
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
            annotation_source = "GO (MsigDB 7.2)/NCBI Gene/UniProt (2020_06)",
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
          dplyr::select(-c(colour, ggcompartment)) %>%
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
    if(is.data.frame(report$data$drug$target$target)){
      if(NROW(report$data$drug$target$target)){
        target_df <- report$data$drug$target$target %>%
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

  if(analysis_output == "loss_of_fitness"){
    if(is.data.frame(report$data$loss_of_fitness$hits_df)){
      if(NROW(report$data$loss_of_fitness$hits_df) > 0){
        target_df <- report$data$loss_of_fitness$hits_df %>%
          dplyr::mutate(
            annotation_source = report$config$resources$projectscore$name,
            version = report$config$resources$projectscore$version) %>%
          dplyr::select(-c(symbol_link_ps, cmp_link)) %>%
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

