get_go_enrichment <- function(query_entrez,
                              background_entrez = NULL,
                              bgset_description = "All protein-coding genes",
                              ontology = "MF",
                              genedb = NULL,
                              p_value_cutoff = 0.05,
                              q_value_cutoff = 0.2,
                              p_value_adjustment_method = "BH",
                              min_geneset_size = 10,
                              max_geneset_size = 500,
                              simplify = F,
                              pool = F) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info(
    paste0("Enrichment - GO: performing gene enrichment analysis of target set (subontology ",ontology,")"))
  lgr::lgr$info( paste0("Enrichment - GO: settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  lgr::lgr$info( paste0("Enrichment - GO: settings: p_value_adjustment_method = ",p_value_adjustment_method))
  lgr::lgr$info( paste0("Enrichment - GO: settings: minGSSize = ",min_geneset_size))
  lgr::lgr$info( paste0("Enrichment - GO: settings: maxGSSize = ",max_geneset_size))
  lgr::lgr$info( paste0("Enrichment - GO: settings: remove redundancy of enriched GO terms = ",simplify))
  lgr::lgr$info( paste0("Enrichment - GO: settings: Background geneset: '",bgset_description,"'"))
  lgr::lgr$info( paste0("Enrichment - GO: settings: Background geneset size = ",length(background_entrez)))


  stopifnot(p_value_adjustment_method %in%
              c("BH","BY","fdr","none","holm",
                "hochberg","hommel","bonferroni"))
  stopifnot(is.character(query_entrez))
  stopifnot(ontology == "MF" | ontology == "BP" |
              ontology == "CC" | ontology == "ALL")
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")
  stopifnot(p_value_cutoff > 0 & p_value_cutoff < 1)
  stopifnot(q_value_cutoff > 0 & q_value_cutoff < 1)
  stopifnot(!is.null(background_entrez))
  stopifnot(is.character(background_entrez))

  ego <-
    suppressMessages(
      clusterProfiler::enrichGO(
        gene          = query_entrez,
        OrgDb         = "org.Hs.eg.db",
        ont           = ontology,
        minGSSize     = min_geneset_size,
        maxGSSize     = max_geneset_size,
        pAdjustMethod = p_value_adjustment_method,
        pvalueCutoff  = p_value_cutoff,
        qvalueCutoff  = q_value_cutoff,
        universe      = background_entrez,
        readable      = F,
        pool = pool)
    )

  if (simplify == T) {
    ego <- clusterProfiler::simplify(ego, cutoff=0.8,
                                     by="p.adjust", select_fun=min)
  }

  df <- as.data.frame(utils::head(ego, 5000))
  rownames(df) <- NULL
  if (ontology == "ALL" & "ONTOLOGY" %in% colnames(df)) {
    df <- df |> dplyr::rename(db = "ONTOLOGY")
  } else {
    df <- df |> dplyr::mutate(db = ontology)
  }
  if (nrow(df) > 0) {
    df <- suppressWarnings(
      df |>
        dplyr::mutate(db = paste0("GO_", .data$db)) |>
        dplyr::rename(
          go_id = "ID",
          go_description = "Description",
          count = "Count",
          gene_ratio = "GeneRatio",
          background_ratio = "BgRatio",
          gene_id = "geneID") |>
        dplyr::mutate(
          go_description_link =
            paste0('<a href=\'http://amigo.geneontology.org/amigo/term/',
                   .data$go_id,'\' target=\'_blank\'>',
                   .data$go_description,'</a>')) |>
        tidyr::separate(.data$gene_ratio,
                        c('num_query_hits','num_query_all'),
                        sep = "/", remove = F, convert = T) |>
        dplyr::mutate(qvalue = as.numeric(.data$qvalue)) |>
        dplyr::mutate(pvalue = as.numeric(.data$pvalue)) |>
        dplyr::mutate(
          qvalue =
            dplyr::if_else(
              !is.na(.data$qvalue),
              as.numeric(formatC(.data$qvalue, format = "e",
                                 digits = 1)),
              as.numeric(NA))) |>
        dplyr::mutate(
          pvalue =
            dplyr::if_else(
              !is.na(.data$pvalue),
              as.numeric(formatC(.data$pvalue, format = "e",
                                 digits = 1)),
              as.numeric(NA))) |>
        tidyr::separate(.data$background_ratio,
                        c('num_background_hits','num_background_all'),
                        sep = "/", remove = F, convert = T) |>
        dplyr::mutate(
          enrichment_factor =
            round(as.numeric((.data$num_query_hits / .data$num_query_all) /
                               (.data$num_background_hits / .data$num_background_all)),
                  digits = 1)) |>
        dplyr::select(-c("num_query_hits",
                         "num_query_all",
                         "num_background_hits",
                         "num_background_all"))
    )

    gene2id <- NULL
    if (!is.null(genedb)) {
      gene2id <- as.data.frame(
        df |>
          dplyr::select(c("go_id", "gene_id")) |>
          tidyr::separate_rows("gene_id", sep = "/") |>
          dplyr::mutate(
            gene_id = as.integer(.data$gene_id)) |>
          dplyr::left_join(
            dplyr::select(genedb,
                          c("entrezgene",
                          "symbol")),
            by=c("gene_id" = "entrezgene"),
            multiple = "all") |>
          dplyr::mutate(
            entrez_url =
              paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                     .data$gene_id, "' target=\'blank_\'>",
                     .data$symbol, "</a>")) |>
          dplyr::group_by(.data$go_id) |>
          dplyr::summarise(
            gene_symbol_link =
              paste(unique(.data$entrez_url), collapse = ", "),
            gene_symbol = paste(unique(.data$symbol), collapse = ", "))
      )
    }
    if (!is.null(gene2id)) {
      df <- df |>
        dplyr::left_join(
          dplyr::select(gene2id,
                        c("go_id",
                        "gene_symbol_link",
                        "gene_symbol")),
          by="go_id", multiple = "all") |>
        dplyr::rename(exact_source = "go_id",
                      description = "go_description",
                      description_link = "go_description_link") |>
        dplyr::mutate(standard_name = .data$exact_source)
    }

    if (!is.null(df)) {
      df$setting_p_value_cutoff <- p_value_cutoff
      df$setting_q_value_cutoff <- q_value_cutoff
      df$setting_p_value_adj_method <- p_value_adjustment_method
      df$setting_min_geneset_size <- min_geneset_size
      df$setting_max_geneset_size <- max_geneset_size
    }
  } else {
    df <- NULL
  }

  return(df)

}


get_universal_enrichment <- function(query_entrez,
                                     background_entrez = NULL,
                                     bgset_description = "All protein-coding genes",
                                     genedb = NULL,
                                     p_value_cutoff = 0.05,
                                     p_value_adjustment_method = "BH",
                                     min_geneset_size = 10,
                                     max_geneset_size = 500,
                                     q_value_cutoff = 0.2,
                                     TERM2GENE = NULL,
                                     TERM2NAME = NULL,
                                     TERM2SOURCE = NULL,
                                     dbsource = "") {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( paste0("Enrichment - ",dbsource,": performing gene enrichment analysis of target set"))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: p_value_adjustment_method = ",p_value_adjustment_method))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: minGSSize = ",min_geneset_size))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: maxGSSize = ",max_geneset_size))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: Background geneset: '",bgset_description,"'"))
  lgr::lgr$info( paste0("Enrichment - ",dbsource,": settings: Background geneset size = ",length(background_entrez)))

  stopifnot(is.character(query_entrez))
  stopifnot(!is.null(background_entrez))
  stopifnot(is.character(background_entrez))
  stopifnot(p_value_adjustment_method %in%
              c("BH","BY","fdr","none","holm",
                "hochberg","hommel","bonferroni"))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")
  stopifnot(p_value_cutoff > 0 & p_value_cutoff < 1)
  stopifnot(q_value_cutoff > 0 & q_value_cutoff < 1)
  stopifnot(!is.null(TERM2GENE))
  stopifnot(!is.null(TERM2NAME))
  stopifnot(is.data.frame(TERM2GENE))
  stopifnot(is.data.frame(TERM2NAME))
  stopifnot(identical(colnames(TERM2NAME), c("standard_name","name")))
  stopifnot(identical(colnames(TERM2GENE), c("standard_name","entrez_gene")))
  stopifnot(is.character(TERM2GENE$entrez_gene))

  stopifnot(NROW(
    dplyr::inner_join(TERM2NAME, TERM2GENE, by = "standard_name", multiple = "all")
  ) > 0)

  enr <-
    suppressMessages(
      clusterProfiler::enricher(gene          = query_entrez,
                                universe      = background_entrez,
                                pAdjustMethod = p_value_adjustment_method,
                                minGSSize     = min_geneset_size,
                                maxGSSize     = max_geneset_size,
                                pvalueCutoff  = p_value_cutoff,
                                qvalueCutoff  = q_value_cutoff,
                                TERM2GENE     = TERM2GENE,
                                TERM2NAME     = TERM2NAME)
    )


  df <- as.data.frame(utils::head(enr,5000))
  rownames(df) <- NULL
  if (nrow(df) > 0) {
    df <- suppressWarnings(
      df |>
        dplyr::rename(
          standard_name = "ID",
          description = "Description",
          count = "Count",
          gene_ratio = "GeneRatio",
          background_ratio = "BgRatio",
          gene_id = "geneID")
    )

    if (dbsource == "WikiPathways") {
      df <- df |>
        dplyr::mutate(
          description_link = paste0(
            "<a href=\"https://www.wikipathways.org/index.php/Pathway:",
            .data$standard_name,"\" target='_blank'>",
            .data$description,"</a>"
          )) |>
        dplyr::mutate(
          exact_source = "https://wikipathways.org",
          external_url = "https://wikipathways.org",
          db = dbsource)
    }
    else if (dbsource == "KEGG") {
      df <- df |>
        dplyr::mutate(
          description_link = paste0(
            "<a href=\"https://www.genome.jp/kegg-bin/show_pathway?",
            stringr::str_replace(.data$standard_name,"hsa","map"),
            "\" target='_blank'>",
            .data$description,"</a>"
          )) |>
        dplyr::mutate(
          exact_source = "https://www.genome.jp/kegg/pathway.html",
          external_url = "https://www.genome.jp/kegg/pathway.html",
          db = dbsource)
    }
    else if (dbsource == "NetPath") {
      df <- df |>
        dplyr::mutate(
          description_link = paste0(
            "<a href=\"http://netpath.org/pathways?path_id=",
            .data$standard_name,"\" target='_blank'>",
            .data$description,"</a>"
          )) |>
        dplyr::mutate(
          exact_source = "http://netpath.org",
          external_url = "http://netpath.org",
          db = dbsource)
    }
    else{
      stopifnot(!is.null(TERM2SOURCE) | !is.data.frame(TERM2SOURCE))
      stopifnot("standard_name" %in% colnames(TERM2SOURCE))
      df <- df |>
        dplyr::left_join(TERM2SOURCE, by="standard_name", multiple = "all") |>
        dplyr::mutate(
          description_link =
            dplyr::if_else(
              !is.na(.data$external_url),
              paste0("<a href='",.data$external_url,
                     "' target='_blank'>",
                     .data$description,"</a>"),
              .data$description))


    }
    df <- suppressWarnings(
      df |>
        dplyr::mutate(qvalue = as.numeric(.data$qvalue)) |>
        dplyr::mutate(pvalue = as.numeric(.data$pvalue)) |>
        dplyr::mutate(
          qvalue =
            dplyr::if_else(
              !is.na(.data$qvalue),
              as.numeric(formatC(.data$qvalue, format = "e",
                                 digits = 1)),
              as.numeric(NA))) |>
        dplyr::mutate(
          pvalue =
            dplyr::if_else(
              !is.na(.data$pvalue),
              as.numeric(formatC(.data$pvalue, format = "e",
                                 digits = 1)),
              as.numeric(NA))) |>
        tidyr::separate(
          .data$gene_ratio,
          c('num_query_hits','num_query_all'),
          sep = "/", remove = F, convert = T) |>
        tidyr::separate(
          .data$background_ratio,
          c('num_background_hits','num_background_all'),
          sep = "/", remove = F, convert = T) |>
        dplyr::mutate(
          enrichment_factor =
            round(as.numeric((
              .data$num_query_hits / .data$num_query_all) /
                (.data$num_background_hits / .data$num_background_all)),
              digits = 1)) |>
        dplyr::select(-c("num_query_hits",
                         "num_query_all",
                         "num_background_hits",
                         "num_background_all")) |>
        dplyr::mutate(
          db = dplyr::if_else(is.na(.data$db) &
                                nchar(dbsource) > 0,
                              dbsource,
                              as.character(.data$db)))
    )

    gene2id <- NULL
    if (!is.null(genedb)) {
      gene2id <- as.data.frame(
        df |>
          dplyr::select(c("standard_name", "gene_id")) |>
          tidyr::separate_rows("gene_id", sep = "/") |>
          dplyr::mutate(gene_id = as.integer(.data$gene_id)) |>
          dplyr::left_join(
            dplyr::select(genedb, c("entrezgene", "symbol")),
            by = c("gene_id" = "entrezgene"), multiple = "all") |>
          dplyr::mutate(
            entrez_url =
              paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                     .data$gene_id,"' target=\'blank_\'>", .data$symbol,"</a>")) |>
          dplyr::group_by(.data$standard_name) |>
          dplyr::summarise(
            gene_symbol_link =
              paste(unique(.data$entrez_url), collapse = ", "),
            gene_symbol = paste(unique(.data$symbol), collapse = ", "))
      )
    }
    if (!is.null(gene2id)) {
      df <- df |>
        dplyr::left_join(
          dplyr::select(gene2id,
                        c("standard_name",
                        "gene_symbol_link",
                        "gene_symbol")),
          by = "standard_name", multiple = "all")
    }
    if (!is.null(df)) {
      df$setting_p_value_cutoff <- p_value_cutoff
      df$setting_q_value_cutoff <- q_value_cutoff
      df$setting_p_value_adj_method <- p_value_adjustment_method
      df$setting_min_geneset_size <- min_geneset_size
      df$setting_max_geneset_size <- max_geneset_size
    }
  } else {
    df <- NULL
  }
  return(df)
}
