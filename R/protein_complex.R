
annotate_protein_complex <- function(query_entrez,
                                     genedb = NULL,
                                     complex_db = NULL,
                                     transcript_xref_db = NULL,
                                     logger = NULL){

  stopifnot(!is.null(logger))
  log4r_info(logger, "OmniPath: retrieval of protein complex information for target set - multiple underlying sources")
  stopifnot(is.integer(query_entrez))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(complex_db))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(complex_db$db, dbtype = "protein_complex")
  validate_db_df(transcript_xref_db, dbtype = "transcript_xref")

  uniprot_xref <- transcript_xref_db %>%
    dplyr::filter(property == "uniprot_acc") %>%
    dplyr::rename(uniprot_acc = value) %>%
    dplyr::select(entrezgene, uniprot_acc)

  all_protein_complexes <- as.data.frame(
    complex_db$up_xref %>%
      dplyr::inner_join(
        uniprot_xref,
        by = "uniprot_acc") %>%
      dplyr::inner_join(
        dplyr::select(genedb, .data$entrezgene, .data$symbol),
        by = "entrezgene") %>%
      dplyr::arrange(.data$complex_id, .data$symbol) %>%
      dplyr::mutate(
        genelink = paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                          .data$entrezgene,
                          "' target='_blank'>", .data$symbol,"</a>")) %>%
      dplyr::group_by(.data$complex_id) %>%
      dplyr::summarise(
        complex_genes = paste(unique(.data$genelink), collapse=", "),
        .groups = "drop")) %>%
    dplyr::left_join(complex_db$db, by = "complex_id") %>%
    dplyr::rename(annotation_source = .data$sources)

  target_complex_genes <- list()


  target_complex_genes[['omnipath']] <-
    data.frame("entrezgene" = query_entrez, stringsAsFactors = F) %>%
    dplyr::left_join(uniprot_xref, by = "entrezgene") %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(.data$uniprot_acc))

  target_complex_genes[['humap2']] <-
    target_complex_genes[['omnipath']]

  target_complex_overlap <- list()

  for(class in c('omnipath', 'humap2')){
    target_complex_overlap[[class]] <- data.frame()

    if(nrow(target_complex_genes[[class]]) > 0){

      target_complex_overlap[[class]] <- as.data.frame(
        complex_db$up_xref %>%
          dplyr::inner_join(target_complex_genes[[class]], by = "uniprot_acc")
      )

      if(nrow(target_complex_overlap[[class]]) > 0){

        target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
          dplyr::left_join(dplyr::select(genedb, entrezgene, symbol),
                           by = "entrezgene")

        if(class == "humap2"){
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
            dplyr::filter(stringr::str_detect(.data$complex_id,"HuMAP2"))
        }else{
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
            dplyr::filter(!stringr::str_detect(.data$complex_id,"HuMAP2"))
        }

        if(nrow(target_complex_overlap[[class]]) > 0){
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
            dplyr::group_by(.data$complex_id) %>%
            dplyr::summarise(target_genes = paste(unique(sort(.data$symbol)),
                                                  collapse = ", ")) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(num_target_members =
                            stringr::str_count(.data$target_genes,",") + 1) %>%
            dplyr::arrange(dplyr::desc(.data$num_target_members)) %>%
            dplyr::left_join(all_protein_complexes, by = "complex_id") %>%
            dplyr::arrange(dplyr::desc(.data$num_target_members),
                           .data$confidence) %>%
            dplyr::select(-.data$num_target_members) %>%
            dplyr::select(.data$complex_name,
                          .data$target_genes,
                          .data$complex_literature_support,
                          .data$complex_genes,
                          .data$annotation_source,
                          .data$disease_comment,
                          .data$complex_comment,
                          .data$confidence,
                          .data$purification_method) %>%
            dplyr::rename(literature = .data$complex_literature_support)


          if(class == "humap2"){
            target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
              dplyr::select(.data$complex_name,
                            .data$target_genes,
                            .data$complex_genes,
                            .data$annotation_source,
                            .data$confidence) %>%
              dplyr::mutate(complex_name = paste0(
                "<a href='http://humap2.proteincomplexes.org/displayComplexes?complex_key=",
                .data$complex_name,
                "' target='_blank'>", .data$complex_name, "</a>"
              ))

          }
        }
      }
    }
  }

  return(target_complex_overlap)

}
