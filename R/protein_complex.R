
annotate_protein_complex <- function(query_entrez,
                                     genedb = NULL,
                                     complex_db = NULL,
                                     complex_up_xref = NULL,
                                     transcript_xref = NULL,
                                     otdb_gene_rank = NULL) {

  lgr::lgr$info( "OmniPath: retrieval of protein complex information for target set - multiple underlying sources")
  stopifnot(is.integer(query_entrez))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(complex_db))
  stopifnot(!is.null(complex_up_xref))
  stopifnot(!is.null(transcript_xref))
  stopifnot(!is.null(otdb_gene_rank))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(complex_db, dbtype = "protein_complex")
  validate_db_df(transcript_xref, dbtype = "transcript_xref")

  otdb_gene_rank <- otdb_gene_rank |>
    dplyr::select(c("ensembl_gene_id",
                  "global_assoc_rank")) |>
    dplyr::distinct()

  uniprot_xref <- transcript_xref |>
    dplyr::filter(.data$property == "uniprot_acc") |>
    dplyr::rename(uniprot_acc = "value") |>
    dplyr::select(c("entrezgene", "uniprot_acc"))


  target_complex_overlap <- list()
  target_complex_overlap[['omnipath']] <- data.frame()
  target_complex_overlap[['humap2']] <- data.frame()

  all_protein_complexes <- as.data.frame(
    complex_up_xref |>
      dplyr::inner_join(
        uniprot_xref,
        by = "uniprot_acc", multiple = "all")
  )

  if(NROW(all_protein_complexes) == 0){
    return(target_complex_overlap)
  }

  all_protein_complexes <- as.data.frame(
    all_protein_complexes |>
      dplyr::inner_join(
        dplyr::select(genedb,
                      c("entrezgene",
                      "symbol",
                      "ensembl_gene_id")),
        by = "entrezgene", multiple = "all") |>
      dplyr::left_join(
        dplyr::select(otdb_gene_rank,
                      c("ensembl_gene_id",
                      "global_assoc_rank")),
        by = "ensembl_gene_id", multiple = "all") |>
      dplyr::mutate(global_assoc_rank = dplyr::if_else(
        is.na(.data$global_assoc_rank),
        as.numeric(0),
        as.numeric(.data$global_assoc_rank)
      )) |>
      dplyr::select(-c("ensembl_gene_id")) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$complex_id, .data$symbol) |>
      dplyr::mutate(
        genelink = paste0(
          "<a href ='http://www.ncbi.nlm.nih.gov/gene/",
          .data$entrezgene,
          "' target='_blank'>", .data$symbol,"</a>")) |>
      dplyr::group_by(.data$complex_id) |>
      dplyr::summarise(
        complex_genes = paste(
          unique(.data$genelink), collapse = ", "),
        num_complex_members = dplyr::n(),
        complex_cancer_rank_score = round(mean(
          .data$global_assoc_rank), digits = 3),
        .groups = "drop")) |>
    dplyr::left_join(complex_db, by = "complex_id", multiple = "all") |>
    dplyr::filter(
      .data$num_complex_members > 2) |>
    dplyr::arrange(
      dplyr::desc(.data$complex_cancer_rank_score)) |>
    dplyr::rename(annotation_source = "sources")

  target_complex_genes <- list()


  target_complex_genes[['omnipath']] <-
    data.frame("entrezgene" = query_entrez, stringsAsFactors = F) |>
    dplyr::left_join(uniprot_xref, by = "entrezgene", multiple = "all") |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$uniprot_acc))

  target_complex_genes[['humap2']] <-
    target_complex_genes[['omnipath']]


  for (class in c('omnipath', 'humap2')) {

    if (nrow(target_complex_genes[[class]]) > 0) {

      target_complex_overlap[[class]] <- as.data.frame(
        complex_up_xref |>
          dplyr::inner_join(target_complex_genes[[class]],
                            by = "uniprot_acc", multiple = "all")
      )

      if (nrow(target_complex_overlap[[class]]) > 0) {

        target_complex_overlap[[class]] <- target_complex_overlap[[class]] |>
          dplyr::left_join(
            dplyr::select(genedb,
                          c("entrezgene",
                          "symbol")),
            by = "entrezgene", multiple = "all")

        if (class == "humap2") {
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] |>
            dplyr::filter(stringr::str_detect(.data$complex_id,"HuMAP2"))
        } else {
          target_complex_overlap[[class]] <- target_complex_overlap[[class]] |>
            dplyr::filter(!stringr::str_detect(.data$complex_id,"HuMAP2"))
        }

        if (nrow(target_complex_overlap[[class]]) > 0) {
          target_complex_overlap[[class]] <- as.data.frame(
            target_complex_overlap[[class]] |>
              dplyr::group_by(.data$complex_id) |>
              dplyr::summarise(
                target_genes = paste(unique(sort(.data$symbol)),
                                     collapse = ", ")) |>
              dplyr::ungroup() |>
              dplyr::mutate(num_target_members =
                              stringr::str_count(.data$target_genes,",") + 1) |>
              #dplyr::arrange(dplyr::desc(.data$num_target_members)) |>
              dplyr::filter(!is.na(.data$complex_id)) |>
              dplyr::inner_join(all_protein_complexes,
                                by = "complex_id", multiple = "all") |>
              dplyr::arrange(dplyr::desc(.data$complex_cancer_rank_score),
                             .data$confidence) |>
              #dplyr::select(-.data$num_target_members) |>
              dplyr::select(c("complex_name",
                            "target_genes",
                            "complex_literature_support",
                            "complex_genes",
                            "annotation_source",
                            "disease_comment",
                            "complex_cancer_rank_score",
                            "num_target_members",
                            "complex_comment",
                            "confidence",
                            "purification_method")) |>
              dplyr::rename(literature = "complex_literature_support")
          )


          if (class == "humap2") {
            target_complex_overlap[[class]] <- target_complex_overlap[[class]] |>
              dplyr::mutate(complex_id = .data$complex_name) |>
              dplyr::select(c("complex_name",
                            "target_genes",
                            "complex_genes",
                            "complex_id",
                            "complex_cancer_rank_score",
                            "num_target_members",
                            "confidence")) |>

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
                                    transcript_xref = NULL) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( "PFAM: retrieval of protein domain information for target set")
  stopifnot(is.integer(query_entrez))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(pfamdb))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(pfamdb, dbtype = "protein_domain")
  validate_db_df(transcript_xref, dbtype = "transcript_xref")

  uniprot_xref <- transcript_xref |>
    dplyr::filter(.data$property == "uniprot_acc") |>
    dplyr::rename(uniprot_acc = "value") |>
    dplyr::select(c("entrezgene", "uniprot_acc"))

  all_protein_domains <- as.data.frame(
    pfamdb |>
      dplyr::inner_join(
        uniprot_xref,
        by = "uniprot_acc", multiple = "all") |>
      dplyr::distinct()
  )

  target_protein_domains <-
    data.frame("entrezgene" = query_entrez,
               stringsAsFactors = F) |>
    dplyr::left_join(
      all_protein_domains,
      by = "entrezgene", multiple = "all") |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$uniprot_acc))


  if (NROW(target_protein_domains) > 0) {
    target_protein_domains <-
      as.data.frame(
        target_protein_domains |>
          dplyr::inner_join(
            dplyr::select(
              genedb,
              c("entrezgene",
              "symbol")),
            by = "entrezgene", multiple = "all") |>
          dplyr::arrange(
            .data$pfam_id,
            .data$symbol) |>
          dplyr::mutate(
            protein_domain =
              paste0(
                "<a href=\"https://www.ebi.ac.uk/interpro/entry/pfam/",
                .data$pfam_id,
                "\" target='_blank'>",
                .data$pfam_long_name,
                "</a>")
          ) |>
          dplyr::mutate(
            pfam_domain_gene =
              paste0(
                "<a href=\"https://www.ebi.ac.uk/interpro/protein/UniProt/",
                .data$uniprot_acc,
                "\" target='_blank'>",
                .data$symbol,
                "</a>:",
                .data$domain_freq
              )
          ) |>
          dplyr::arrange(
            .data$protein_domain,
            dplyr::desc(.data$domain_freq)
          ) |>
          dplyr::group_by(
            .data$protein_domain,
            ) |>
          dplyr::summarise(
            target_genes = paste(
              unique(.data$pfam_domain_gene),
              collapse = ", "
            ),
            target_set_frequency =
              sum(.data$domain_freq),
            .groups = "drop")
      ) |>
      dplyr::mutate(
        source_db = "Pfam"
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$target_set_frequency))
  }

  return(target_protein_domains)

}

