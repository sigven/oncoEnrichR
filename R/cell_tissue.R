gene_tissue_cell_spec_cat <-
  function(qgenes,
           q_id_type = "symbol",
           resolution = "tissue",
           genedb = NULL,
           oeDB = NULL,
           logger = NULL){

    stopifnot(!is.null(logger))
    etype <- "tissue"
    edb <- "GTex"
    if(resolution != "tissue"){
      etype <- "cell type"
      edb <- "HPA"
    }
    log4r_info(logger,
               paste0("Retrieving ", etype,
                      " specificity (", edb,
                      ") category of target genes")
    )

    stopifnot(!is.null(genedb))
    stopifnot(!is.null(oeDB$tissuecelldb))
    validate_db_df(genedb, dbtype = "genedb")
    stopifnot(resolution == "tissue" | resolution == "single_cell")
    stopifnot(q_id_type == "symbol" | q_id_type == "entrezgene")
    query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
    if(q_id_type == 'entrezgene'){
      stopifnot(is.integer(qgenes))
      query_genes_df <-
        data.frame('entrezgene' = qgenes, stringsAsFactors = F) %>%
        dplyr::inner_join(genedb, by = "entrezgene") %>%
        dplyr::distinct()

    }else{
      stopifnot(is.character(qgenes))
      query_genes_df <- query_genes_df %>%
        dplyr::inner_join(genedb, by = "symbol") %>%
        dplyr::distinct()
    }
    specificity_categories <-
      c('Group enriched',
        'Low tissue specificity',
        'Not detected',
        'Mixed',
        'Tissue enhanced',
        'Tissue enriched')
    source = 'tissue'
    if(resolution == "single_cell"){
      specificity_categories <-
        c('Group enriched',
          'Low cell type specificity',
          'Not detected',
          'Mixed',
          'Cell type enhanced',
          'Cell type enriched')
      source = 'single_cell'
    }

    specificity_groups_target <- as.data.frame(
      oeDB$tissuecelldb[[source]][['te_df']] %>%
        dplyr::inner_join(
          dplyr::select(query_genes_df,
                        .data$symbol,
                        .data$ensembl_gene_id,
                        .data$name),
          by = "ensembl_gene_id")
    )


    if(resolution == "tissue"){
      specificity_groups_target <-
        specificity_groups_target %>%
        dplyr::mutate(
          genename = paste0(
            "<a href='https://gtexportal.org/home/gene/",
            .data$ensembl_gene_id,"' target='_blank'>",
            .data$name,"</a>")
        )
    }else{
      specificity_groups_target <-
        specificity_groups_target %>%
        dplyr::mutate(
          genename = paste0(
            "<a href='https://www.proteinatlas.org/",
            .data$ensembl_gene_id,
            "-",
            .data$symbol,
            "/celltype' target='_blank'>",
            .data$name,"</a>")
        )
    }
    specificity_groups_target <- as.data.frame(
      specificity_groups_target %>%
      dplyr::group_by(.data$category) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(tot = sum(.data$n)) %>%
      dplyr::mutate(pct = round((.data$n / .data$tot) * 100, digits = 2)) %>%
      dplyr::mutate(
        group = paste0("Target set (n = ", .data$tot,")"
        )
      )
    )

    tot <- unique(specificity_groups_target$tot)

    for(n in specificity_categories){
      df <- data.frame(
        'category' = n,
        'pct' = 0,
        'n' = 0,
        'tot' = 0,
        group = paste0("Target set (n = ",
                       formatC(tot, format="f",
                               big.mark = ",", digits=0),")"))
      if(nrow(dplyr::inner_join(
        df,
        specificity_groups_target,
        by = "category")) == 0){
        specificity_groups_target <-
          specificity_groups_target %>%
          dplyr::bind_rows(df)

      }
    }

    specificity_groups_all <- as.data.frame(
      oeDB$tissuecelldb[[source]][['te_df']] %>%
        dplyr::group_by(.data$category) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(tot = sum(.data$n)) %>%
        dplyr::mutate(pct = round((.data$n / .data$tot) * 100, digits = 2)) %>%
        dplyr::mutate(group =
                        paste0("All HPA proteins (n = ",
                               formatC(.data$tot, format="f",
                                       big.mark = ",",
                                       digits=0),")"
                        )
        )
    )

    category_df <- dplyr::bind_rows(
      specificity_groups_all,
      specificity_groups_target
    )

    exp_dist_df <- data.frame()
    exp_dist_df <- oeDB$tissuecelldb[[resolution]]$expr_df
    exp_dist_df$ensembl_gene_id <-
      rownames(exp_dist_df)

    if(resolution == "tissue"){
      exp_dist_df <- as.data.frame(
        exp_dist_df %>%
          tidyr::pivot_longer(cols = !tidyr::starts_with("ensembl"),
                              names_to = "tissue",
                              values_to = "nTPM") %>%
          dplyr::inner_join(
            dplyr::select(query_genes_df,
                          .data$symbol,
                          .data$ensembl_gene_id),
            by = "ensembl_gene_id") %>%
          dplyr::mutate(exp = round(log2(.data$nTPM), digits = 3)) %>%
          dplyr::mutate(exp_measure = "log2(nTPM)")
      )
    }else{
      exp_dist_df <- as.data.frame(
        exp_dist_df %>%
          tidyr::pivot_longer(cols = !tidyr::starts_with("ensembl"),
                              names_to = "cell_type",
                              values_to = "nTPM") %>%
          dplyr::inner_join(
            dplyr::select(query_genes_df,
                          .data$symbol,
                          .data$ensembl_gene_id),
            by = "ensembl_gene_id") %>%
          dplyr::mutate(exp = round(log2(.data$nTPM), digits = 3)) %>%
          dplyr::mutate(exp_measure = "log2(nTPM)")
      )
    }

    return(list('category_df' = category_df,
                'exp_dist_df' = exp_dist_df))

}


gene_tissue_cell_enrichment <-
  function(qgenes_entrez,
           background_entrez = NULL,
           genedb = NULL,
           oeDB = NULL,
           resolution = "tissue",
           logger = logger){

    stopifnot(!is.null(logger))
    etype <- "tissues"
    edb <- "GTex"
    if(resolution != "tissue"){
      etype <- "cell types"
      edb <- "HPA"
    }
    log4r_info(logger,
      paste0("Estimating enrichment of ", etype,
             " (", edb,
             ") in target set with TissueEnrich"))

    stopifnot(!is.null(qgenes_entrez) &
                is.integer(qgenes_entrez))
    stopifnot(!is.null(genedb) |
                !is.data.frame(genedb))
    stopifnot(!is.null(oeDB$tissuecelldb))
    stopifnot(resolution == "tissue" | resolution == "single_cell")
    assertable::assert_colnames(
      genedb, c("ensembl_gene_id",
                "entrezgene",
                "gene_biotype",
                "cancer_max_rank"),
      only_colnames = F, quiet = T)

    df <- data.frame('entrezgene' = as.integer(qgenes_entrez),
                     stringsAsFactors = F) %>%
      dplyr::left_join(
        dplyr::select(genedb,
                      .data$entrezgene,
                      .data$ensembl_gene_id,
                      .data$symbol,
                      .data$name,
                      .data$cancer_max_rank),
        by = "entrezgene")

    if(resolution == "tissue"){
      df <- df %>%
        dplyr::mutate(
          genename = paste0(
            "<a href='https://gtexportal.org/home/gene/",
            .data$ensembl_gene_id,"' target='_blank'>",
            .data$name,"</a>")
        )
    }else{
      df <- df %>%
        dplyr::mutate(
          genename = paste0(
            "<a href='https://www.proteinatlas.org/",
            .data$ensembl_gene_id,
            "-",
            .data$symbol,
            "/celltype' target='_blank'>",
            .data$name,"</a>")
        )
    }



    bg <- oeDB$tissuecelldb[[resolution]][['te_df']]

    q <- bg %>%
      dplyr::select(.data$ensembl_gene_id) %>%
      dplyr::inner_join(df, by = "ensembl_gene_id") %>%
      dplyr::distinct()
    query_ensembl <- q$ensembl_gene_id

    specificities_per_gene <- bg %>%
      dplyr::filter(!is.na(.data$ensembl_gene_id)) %>%
      dplyr::inner_join(df, by = "ensembl_gene_id")

    if(nrow(specificities_per_gene) > 0){

      if(resolution == "tissue"){
        specificities_per_gene <- specificities_per_gene %>%
          dplyr::select(.data$symbol,
                        .data$genename,
                        .data$category,
                        .data$tissue,
                        .data$cancer_max_rank) %>%
          dplyr::arrange(.data$category,
                         dplyr::desc(.data$cancer_max_rank))
      }else{
        specificities_per_gene <- specificities_per_gene %>%
          dplyr::select(.data$symbol,
                        .data$genename,
                        .data$category,
                        .data$cell_type,
                        .data$cancer_max_rank) %>%
          dplyr::arrange(.data$category,
                         dplyr::desc(.data$cancer_max_rank))
      }
    }

    background_ensembl <- bg$ensembl_gene_id
    if(!is.null(background_entrez)){
      df <-
        data.frame('entrezgene' = as.integer(background_entrez),
                   stringsAsFactors = F) %>%
        dplyr::inner_join(
          dplyr::select(genedb, ensembl_gene_id, entrezgene),
          by = "entrezgene"
        )
      bg <- bg %>%
        dplyr::inner_join(df, by = "ensembl_gene_id") %>%
        dplyr::distinct()
      background_ensembl <- bg$ensembl_gene_id

      background_ensembl <- unique(
        c(background_ensembl, query_ensembl)
      )
    }

    gene_id_TE <- GSEABase::ENSEMBLIdentifier()
    gs_query <- GSEABase::GeneSet(
      geneIds = query_ensembl,
      organism = "Homo Sapiens",
      geneIdType = gene_id_TE)

    gs_background <- GSEABase::GeneSet(
      geneIds = background_ensembl,
      organism = "Homo Sapiens",
      geneIdType = gene_id_TE)

    te_output <- NULL

    ## get pre-defined tissue-specificities (SummarisedExperiement) - created with
    ## tissueEnrich::teGeneRetrieval on expression data sets
    se <- oeDB$tissuecelldb[[resolution]][['te_SE']]

    ## perform tissue enrichment analysis for query dataset (gs_query),
    ## using gs_background as the background dataset, and the annotated
    ## tissue/cell-type-specific genes in se as basis for detection
    ## of enriched tissues/cell types
    suppressMessages(
      te_output <- TissueEnrich::teEnrichmentCustom(
        inputGenes = gs_query,
        tissueSpecificGenes = se,
        backgroundGenes = gs_background)
    )

    enrichment_df <- data.frame()

    if(!is.null(te_output)){
      if(!is.null(te_output[[1]])){
        se_enrich_output <- te_output[[1]]
        enrichment_df <- stats::setNames(
          data.frame(SummarizedExperiment::assay(se_enrich_output),
                     row.names = SummarizedExperiment::rowData(se_enrich_output)[,1]),
          SummarizedExperiment::colData(se_enrich_output)[,1])

        enrichment_df$Tissue <- rownames(enrichment_df)
        rownames(enrichment_df) <- NULL
        enrichment_df <- enrichment_df %>%
          dplyr::filter(.data$Tissue != "All") %>%
          dplyr::rename(fold_change = .data$fold.change,
                        tissue = .data$Tissue,
                        tissue_specific_genes = .data$Tissue.Specific.Genes,
                        log10_pvalue = .data$Log10PValue) %>%
          dplyr::select(.data$tissue,
                        .data$tissue_specific_genes,
                        .data$fold_change,
                        .data$log10_pvalue) %>%
          dplyr::arrange(dplyr::desc(.data$fold_change))

        if(resolution == "single_cell"){
          enrichment_df <- enrichment_df %>%
            dplyr::rename(cell_type = .data$tissue,
                          celltype_specific_genes =
                            .data$tissue_specific_genes)
        }

      }
    }

    return(list('per_gene' = specificities_per_gene,
                'per_type' = enrichment_df))

  }

