
test_that("GO enrichment/ORA functionality (clusterProfiler) - testing ", {
  expect_error(oncoEnrichR:::get_go_enrichment(
    query_entrez = "200"))
  expect_error(oncoEnrichR:::get_go_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez,
    p_value_adjustment_method = "UNKNOWN"))
  expect_error(oncoEnrichR:::get_go_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez,
    p_value_cutoff = 2))
  expect_error(oncoEnrichR:::get_go_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez,
    q_value_cutoff = 2))
  expect_error(oncoEnrichR:::get_go_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez,
    ontology = "UNKNOWN"))

   expect_error(
     oncoEnrichR:::get_go_enrichment(
       query_entrez = as.integer(
         c(3845, 4609)),
       genedb = oedb$genedb$all,
       background_entrez = background_full_entrez,
     )
   )

   expect_true(
     is.data.frame(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         background_entrez = background_full_entrez,

       )
     )
   )

   expect_gte(
     NROW(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         background_entrez = background_full_entrez,

       )
     ),
     as.integer(0)
   )


   expect_identical(
     colnames(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         background_entrez = background_full_entrez,
       )
     ),
     c("exact_source", "description",
       "gene_ratio", "background_ratio",
       "rich_factor","fold_enrichment","z_score",
       "pvalue", "p.adjust","qvalue",
       "gene_id","count","db",
       "description_link",
       "enrichment_factor",
       "gene_symbol_link",
       "gene_symbol",
       "standard_name",
       "setting_p_value_cutoff",
       "setting_q_value_cutoff",
       "setting_p_value_adj_method",
       "setting_min_geneset_size",
       "setting_max_geneset_size")
   )

   expect_identical(
     oncoEnrichR:::get_go_enrichment(
       query_entrez = c("3845", "4609"),
       background_entrez = utils::head(
         background_sample_entrez, 20),
       genedb = oedb$genedb$all
     ),
     NULL
   )

})

test_that("Universal enrichment/ORA functionality (clusterProfiler) - testing ", {
  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200"))
  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    background_entrez = background_full_entrez))

  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    p_value_adjustment_method = "UNKNOWN",
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez))

  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    p_value_cutoff =  2,
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez))

  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    q_value_cutoff =  2,
    genedb = oedb$genedb$all,
    background_entrez = background_full_entrez))

  for(c in c("WikiPathways","NetPath","KEGG")){

    pway_db <- tolower(c)

    expect_gte(
      NROW(
        oncoEnrichR:::get_universal_enrichment(
          query_entrez = as.character(
            myc_data$entrezgene),
          background_entrez = background_full_entrez,
          genedb = oedb$genedb$all,
          dbsource = c,
          TERM2GENE = oedb$pathwaydb[[pway_db]]$TERM2GENE,
          TERM2NAME = oedb$pathwaydb[[pway_db]]$TERM2NAME)
      ),
      1
    )

    expect_identical(
      colnames(
        oncoEnrichR:::get_universal_enrichment(
          query_entrez = as.character(
            myc_data$entrezgene),
          genedb = oedb$genedb$all,
          background_entrez = background_full_entrez,
          dbsource = c,
          TERM2GENE = oedb$pathwaydb[[pway_db]]$TERM2GENE,
          TERM2NAME = oedb$pathwaydb[[pway_db]]$TERM2NAME)
      ),
      c("standard_name","description","gene_ratio","background_ratio",
        "rich_factor","fold_enrichment","z_score",
        "pvalue","p.adjust","qvalue","gene_id","count","description_link",
        "exact_source",
        "external_url","url","db","enrichment_factor",
        "gene_symbol_link","gene_symbol","setting_p_value_cutoff",
        "setting_q_value_cutoff","setting_p_value_adj_method",
        "setting_min_geneset_size","setting_max_geneset_size")
    )

  }

  expect_gte(
    NROW(
      oncoEnrichR:::get_universal_enrichment(
        query_entrez = as.character(
          myc_data$entrezgene),
        background_entrez = background_full_entrez,
        genedb = oedb$genedb$all,
        dbsource = "MSIGDB/HALLMARK",
        TERM2GENE = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2GENE,
        TERM2NAME = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2NAME,
        TERM2SOURCE = oedb$pathwaydb$msigdb$TERM2SOURCE)
    ),
    1
  )

  expect_identical(
    colnames(
      oncoEnrichR:::get_universal_enrichment(
        query_entrez = as.character(
          myc_data$entrezgene),
        genedb = oedb$genedb$all,
        background_entrez = background_full_entrez,
        dbsource = "MSIGDB/HALLMARK",
        TERM2GENE = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2GENE,
        TERM2NAME = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2NAME,
        TERM2SOURCE = oedb$pathwaydb$msigdb$TERM2SOURCE)
    ),
    c("standard_name","description","gene_ratio","background_ratio",
      "rich_factor","fold_enrichment","z_score",
      "pvalue","p.adjust","qvalue","gene_id","count","exact_source",
      "external_url","url","db","description_link","enrichment_factor",
      "gene_symbol_link","gene_symbol","setting_p_value_cutoff",
      "setting_q_value_cutoff","setting_p_value_adj_method",
      "setting_min_geneset_size","setting_max_geneset_size")
  )

  expect_identical(
    oncoEnrichR:::get_universal_enrichment(
      query_entrez = as.character(
        myc_data$entrezgene),
      genedb = oedb$genedb$all,
      background_entrez = background_sample_entrez,
      dbsource = "NetPath",
      TERM2GENE = oedb$pathwaydb$netpath$TERM2GENE,
      TERM2NAME = oedb$pathwaydb$netpath$TERM2NAME),
    NULL
  )

})
