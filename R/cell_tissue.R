gene_tissue_cell_spec_cat <-
  function(qgenes,
           qsource = "symbol",
           resolution = "tissue",
           genedb = NULL){

    etype <- "tissue"
    edb <- "GTex"
    if(resolution != "tissue"){
      etype <- "cell type"
      edb <- "HPA"
    }
    rlogging::message(
      paste0("Retrieving ", etype,
             " specificity (", edb,
             ") category of target genes")
    )

  stopifnot(!is.null(genedb))
  oncoEnrichR:::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(resolution == "tissue" | resolution == "single_cell")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == 'entrezgene'){
    query_genes_df <-
      data.frame('entrezgene' = qgenes, stringsAsFactors = F) %>%
      dplyr::inner_join(genedb, by = "entrezgene") %>%
      dplyr::distinct()

  }else{
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
    oncoEnrichR::tissue_cell_expr[[source]][['te_df']] %>%
      dplyr::inner_join(
        dplyr::select(query_genes_df, symbol,
                      ensembl_gene_id, genename),
        by = "ensembl_gene_id") %>%
      dplyr::group_by(category) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(tot = sum(n)) %>%
      dplyr::mutate(pct = round((n / tot) * 100, digits = 2)) %>%
      dplyr::mutate(
        group = paste0("Target set (n = ", tot,")"
          )
      )
  )

  tot <- unique(specificity_groups_target$tot)

  for(n in specificity_categories){
    df <- data.frame(
      'category' = n, 'pct' = 0,
      'n' = 0, 'tot' = 0,
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
    oncoEnrichR::tissue_cell_expr[[source]][['te_df']] %>%
      dplyr::group_by(category) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(tot = sum(n)) %>%
      dplyr::mutate(pct = round((n / tot) * 100, digits = 2)) %>%
      dplyr::mutate(group =
                      paste0("All HPA proteins (n = ",
                             formatC(tot, format="f",
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
  if(resolution == "tissue"){
    exp_dist_df <- as.data.frame(
      oncoEnrichR::tissue_cell_expr[[resolution]]$expr_df %>%
        dplyr::mutate(ensembl_gene_id = rownames(.)) %>%
        tidyr::pivot_longer(cols = !starts_with("ensembl"),
                     names_to = "tissue",
                     values_to = "NX") %>%
        dplyr::inner_join(
          dplyr::select(query_genes_df, symbol,
                        ensembl_gene_id),
          by = "ensembl_gene_id") %>%
        dplyr::mutate(exp = round(log2(NX), digits = 3)) %>%
        dplyr::mutate(exp_measure = "log2(NX)")
      )
  }else{
    exp_dist_df <- as.data.frame(
      oncoEnrichR::tissue_cell_expr[[resolution]]$expr_df %>%
        dplyr::mutate(ensembl_gene_id = rownames(.)) %>%
        tidyr::pivot_longer(cols = !starts_with("ensembl"),
                     names_to = "cell_type",
                     values_to = "TPM") %>%
        dplyr::inner_join(
          dplyr::select(query_genes_df, symbol,
                        ensembl_gene_id),
          by = "ensembl_gene_id") %>%
        dplyr::mutate(exp = round(log2(TPM), digits = 3)) %>%
        dplyr::mutate(exp_measure = "log2(TPM)")
    )
  }

  return(list('category_df' = category_df,
              'exp_dist_df' = exp_dist_df))

}


gene_tissue_cell_enrichment <-
  function(qgenes_entrez,
           background_entrez = NULL,
           genedb = NULL,
           resolution = "tissue"){

    etype <- "tissues"
    edb <- "GTex"
    if(resolution != "tissue"){
      etype <- "cell types"
      edb <- "HPA"
    }
    rlogging::message(
      paste0("Estimating enrichment of ", etype,
             " (", edb,
             ") in target set with TissueEnrich"))

    stopifnot(!is.null(qgenes_entrez) &
                is.character(qgenes_entrez))
    stopifnot(!is.null(genedb) |
                !is.data.frame(genedb))
    stopifnot(resolution == "tissue" | resolution == "single_cell")
    assertable::assert_colnames(
      genedb, c("ensembl_gene_id","entrezgene",
                "gencode_gene_biotype"),
      only_colnames = F, quiet = T)

    df <- data.frame('entrezgene' = as.character(qgenes_entrez),
                     stringsAsFactors = F)

    bg <- oncoEnrichR::tissue_cell_expr[[resolution]][['te_df']]

    q <- bg %>%
      dplyr::select(ensembl_gene_id, entrezgene) %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::inner_join(df, by = "entrezgene") %>%
      dplyr::distinct()
    query_ensembl <- q$ensembl_gene_id

    specificities_per_gene <- bg %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::inner_join(df, by = "entrezgene")
    if(nrow(specificities_per_gene) > 0){

      if(resolution == "tissue"){
        specificities_per_gene <- specificities_per_gene %>%
          dplyr::select(symbol, genename, category, tissue) %>%
          dplyr::arrange(category)
      }else{
        specificities_per_gene <- specificities_per_gene %>%
          dplyr::select(symbol, genename, category, cell_type) %>%
          dplyr::arrange(category)
      }
    }

    background_ensembl <- bg$ensembl_gene_id
    if(!is.null(background_entrez)){
      df <-
        data.frame('entrezgene' = as.character(background_entrez),
                   stringsAsFactors = F)
      bg <- bg %>%
        dplyr::inner_join(df, by = "entrezgene") %>%
        dplyr::distinct()
      background_ensembl <- bg$ensembl_gene_id
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
    se <- oncoEnrichR::tissue_cell_expr[[resolution]][['te_SE']]

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
        enrichment_df <- setNames(
          data.frame(assay(se_enrich_output),
                     row.names = rowData(se_enrich_output)[,1]),
          colData(se_enrich_output)[,1]) %>%
          dplyr::mutate(Tissue = rownames(.)) %>%
          magrittr::set_rownames(NULL) %>%
          dplyr::filter(Tissue != "All") %>%
          dplyr::rename(fold_change = fold.change,
                        tissue = Tissue,
                        tissue_specific_genes = Tissue.Specific.Genes,
                        log10_pvalue = Log10PValue) %>%
          dplyr::select(tissue, tissue_specific_genes,
                        fold_change, log10_pvalue) %>%
          dplyr::arrange(desc(fold_change))

        if(resolution == "single_cell"){
          enrichment_df <- enrichment_df %>%
            dplyr::rename(cell_type = tissue,
                          celltype_specific_genes =
                            tissue_specific_genes)
        }

      }
    }

    return(list('per_gene' = specificities_per_gene,
                'per_type' = enrichment_df))

  }

