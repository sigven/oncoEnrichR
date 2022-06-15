
annotate_protein_complex <- function(query_entrez,
                                     genedb = NULL,
                                     complex_db = NULL,
                                     complex_up_xref = NULL,
                                     transcript_xref = NULL,
                                     otdb_gene_rank = NULL,
                                     logger = NULL){

  stopifnot(!is.null(logger))
  log4r_info(logger, "OmniPath: retrieval of protein complex information for target set - multiple underlying sources")
  stopifnot(is.integer(query_entrez))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(complex_db))
  stopifnot(!is.null(complex_up_xref))
  stopifnot(!is.null(transcript_xref))
  stopifnot(!is.null(otdb_gene_rank))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(complex_db, dbtype = "protein_complex")
  validate_db_df(transcript_xref, dbtype = "transcript_xref")

  otdb_gene_rank <- otdb_gene_rank %>%
    dplyr::select(.data$ensembl_gene_id,
                  .data$global_assoc_rank) %>%
    dplyr::distinct()

  uniprot_xref <- transcript_xref %>%
    dplyr::filter(.data$property == "uniprot_acc") %>%
    dplyr::rename(uniprot_acc = .data$value) %>%
    dplyr::select(.data$entrezgene, .data$uniprot_acc)

  all_protein_complexes <- as.data.frame(
    complex_up_xref %>%
      dplyr::inner_join(
        uniprot_xref,
        by = "uniprot_acc") %>%
      dplyr::inner_join(
        dplyr::select(genedb,
                      .data$entrezgene,
                      .data$symbol,
                      .data$ensembl_gene_id),
        by = "entrezgene") %>%
      dplyr::left_join(
        dplyr::select(otdb_gene_rank,
                      .data$ensembl_gene_id,
                      .data$global_assoc_rank),
        by = "ensembl_gene_id") %>%
      dplyr::mutate(global_assoc_rank = dplyr::if_else(
        is.na(.data$global_assoc_rank),
        as.numeric(0),
        as.numeric(.data$global_assoc_rank)
      )) %>%
      dplyr::select(-.data$ensembl_gene_id) %>%
      dplyr::arrange(.data$complex_id, .data$symbol) %>%
      dplyr::mutate(
        genelink = paste0(
          "<a href ='http://www.ncbi.nlm.nih.gov/gene/",
          .data$entrezgene,
          "' target='_blank'>", .data$symbol,"</a>")) %>%
      dplyr::group_by(.data$complex_id) %>%
      dplyr::summarise(
        complex_genes = paste(
          unique(.data$genelink), collapse=", "),
        num_complex_members = dplyr::n(),
        complex_cancer_rank_score = round(mean(
          .data$global_assoc_rank), digits = 3),
        .groups = "drop")) %>%
    dplyr::left_join(complex_db, by = "complex_id") %>%
    dplyr::filter(
      .data$num_complex_members > 2) %>%
    dplyr::arrange(
      dplyr::desc(.data$complex_cancer_rank_score)) %>%
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
        complex_up_xref %>%
          dplyr::inner_join(target_complex_genes[[class]], by = "uniprot_acc")
      )

      if(nrow(target_complex_overlap[[class]]) > 0){

        target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
          dplyr::left_join(
            dplyr::select(genedb,
                          .data$entrezgene,
                          .data$symbol),
            by = "entrezgene")

        if(class == "humap2"){
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
            dplyr::filter(stringr::str_detect(.data$complex_id,"HuMAP2"))
        }else{
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
            dplyr::filter(!stringr::str_detect(.data$complex_id,"HuMAP2"))
        }

        if(nrow(target_complex_overlap[[class]]) > 0){
          target_complex_overlap[[class]] <- as.data.frame(
            target_complex_overlap[[class]] %>%
              dplyr::group_by(.data$complex_id) %>%
              dplyr::summarise(
                target_genes = paste(unique(sort(.data$symbol)),
                                     collapse = ", ")) %>%
              dplyr::ungroup() %>%
              dplyr::mutate(num_target_members =
                              stringr::str_count(.data$target_genes,",") + 1) %>%
              #dplyr::arrange(dplyr::desc(.data$num_target_members)) %>%
              dplyr::filter(!is.na(.data$complex_id)) %>%
              dplyr::inner_join(all_protein_complexes, by = "complex_id") %>%
              dplyr::arrange(dplyr::desc(.data$complex_cancer_rank_score),
                             .data$confidence) %>%
              #dplyr::select(-.data$num_target_members) %>%
              dplyr::select(.data$complex_name,
                            .data$target_genes,
                            .data$complex_literature_support,
                            .data$complex_genes,
                            .data$annotation_source,
                            .data$disease_comment,
                            .data$complex_cancer_rank_score,
                            .data$num_target_members,
                            .data$complex_comment,
                            .data$confidence,
                            .data$purification_method) %>%
              dplyr::rename(literature = .data$complex_literature_support)
          )


          if(class == "humap2"){
            target_complex_overlap[[class]] <- target_complex_overlap[[class]] %>%
              dplyr::mutate(complex_id = .data$complex_name) %>%
              dplyr::select(.data$complex_name,
                            .data$target_genes,
                            .data$complex_genes,
                            .data$complex_id,
                            .data$complex_cancer_rank_score,
                            .data$num_target_members,
                            .data$confidence) %>%

              dplyr::mutate(complex_id = paste0(
                "<a href='http://humap2.proteincomplexes.org/displayComplexes?complex_key=",
                .data$complex_id,
                "' target='_blank'>", .data$complex_id, "</a>"
              ))

          }
        }
      }
    }
  }

  return(target_complex_overlap)

}

annotate_protein_domain <- function(query_entrez,
                                    genedb = NULL,
                                    pfamdb = NULL,
                                    transcript_xref = NULL,
                                    logger = NULL){

  stopifnot(!is.null(logger))
  log4r_info(logger, "PFAM: retrieval of protein domain information for target set")
  stopifnot(is.integer(query_entrez))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(pfamdb))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(pfamdb, dbtype = "protein_domain")
  validate_db_df(transcript_xref, dbtype = "transcript_xref")

  uniprot_xref <- transcript_xref %>%
    dplyr::filter(.data$property == "uniprot_acc") %>%
    dplyr::rename(uniprot_acc = .data$value) %>%
    dplyr::select(.data$entrezgene, .data$uniprot_acc)

  all_protein_domains <- as.data.frame(
    pfamdb %>%
      dplyr::inner_join(
        uniprot_xref,
        by = "uniprot_acc")
  )

  target_protein_domains <-
    data.frame("entrezgene" = query_entrez, stringsAsFactors = F) %>%
    dplyr::left_join(all_protein_domains, by = "entrezgene") %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(.data$uniprot_acc))


  if(NROW(target_protein_domains) > 0){
    target_protein_domains <-
      as.data.frame(
        target_protein_domains %>%
          dplyr::inner_join(
            dplyr::select(
              genedb,
              .data$entrezgene,
              .data$symbol),
            by = "entrezgene") %>%
          dplyr::arrange(
            .data$pfam_id,
            .data$symbol) %>%
          dplyr::mutate(
            protein_domain =
              paste0(
                "<a href=\"http://pfam.xfam.org/family/",
                .data$pfam_id,
                "\" target='_blank'>",
                .data$pfam_long_name,
                "</a>")
          ) %>%
          dplyr::mutate(
            pfam_domain_gene =
              paste0(
                "<a href=\"http://pfam.xfam.org/protein/",
                .data$uniprot_acc,
                "\" target='_blank'>",
                .data$symbol,
                "</a>:",
                .data$domain_freq
              )
          ) %>%
          dplyr::arrange(
            .data$protein_domain,
            dplyr::desc(.data$domain_freq)
          ) %>%
          dplyr::group_by(
            .data$protein_domain,
            ) %>%
          dplyr::summarise(
            target_genes = paste(
              unique(.data$pfam_domain_gene),
              collapse = ", "
            ),
            target_set_frequency =
              sum(.data$domain_freq),
            .groups = "drop")
      ) %>%
      dplyr::mutate(
        source_db = "Pfam"
      ) %>%
      dplyr::arrange(
        dplyr::desc(.data$target_set_frequency))
  }

  return(target_protein_domains)

}

