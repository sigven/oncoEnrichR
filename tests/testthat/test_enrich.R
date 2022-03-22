
test_that("GO enrichment/ORA functionality (clusterProfiler) - testing ", {
   expect_error(oncoEnrichR:::get_go_enrichment(query_entrez = NULL))
   expect_error(oncoEnrichR:::get_go_enrichment(genedb = NULL))
   expect_error(oncoEnrichR:::get_go_enrichment(ontology = "UNKNOWN"))
   expect_error(oncoEnrichR:::get_go_enrichment(logger = NULL))
   expect_error(oncoEnrichR:::get_go_enrichment(p_value_adjustment_method = "UNKNOWN"))
   expect_error(oncoEnrichR:::get_go_enrichment(p_value_cutoff = 2))
   expect_error(oncoEnrichR:::get_go_enrichment(q_value_cutoff = 2))


   expect_error(
     oncoEnrichR:::get_go_enrichment(
       query_entrez = as.integer(
         c(3845, 4609)),
       genedb = oedb$genedb$all,
       logger = log4r_logger,
       ontology = "MF"
     )
   )

   expect_true(
     is.data.frame(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         logger = log4r_logger,
         ontology = "MF"
       )
     )
   )

   expect_gte(
     NROW(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         logger = log4r_logger,
         ontology = "MF"
       )
     ),
     as.integer(0)
   )


   expect_identical(
     colnames(
       oncoEnrichR:::get_go_enrichment(
         query_entrez = c("3845", "4609"),
         genedb = oedb$genedb$all,
         logger = log4r_logger,
         ontology = "MF"
       )
     ),
     c("exact_source", "description",
       "gene_ratio", "background_ratio",
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

   expect_error(
     oncoEnrichR:::get_go_enrichment(
       query_entrez = c("3845", "4609"),
       background_entrez = background_sample_set$entrezgene,
       genedb = oedb$genedb$all,
       logger = log4r_logger,
       ontology = "MF"
     )
   )

   expect_identical(
     oncoEnrichR:::get_go_enrichment(
       query_entrez = c("3845", "4609"),
       background_entrez = as.character(
         head(background_sample_set, 20)),
       genedb = oedb$genedb$all,
       logger = log4r_logger,
       ontology = "MF"
     ),
     NULL
   )

})

test_that("Universal enrichment/ORA functionality (clusterProfiler) - testing ", {
  expect_error(oncoEnrichR:::get_universal_enrichment(logger = NULL))
  expect_error(oncoEnrichR:::get_universal_enrichment(
    p_value_adjustment_method = "UNKNOWN",
    query_entrez = "200",
    logger = log4r_logger))
  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    genedb = NULL,
    logger = log4r_logger))
  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    p_value_cutoff = 2,
    logger = log4r_logger))
  expect_error(oncoEnrichR:::get_universal_enrichment(
    query_entrez = "200",
    genedb = oedb$genedb$all,
    q_value_cutoff = 2,
    logger = log4r_logger))

  for(c in c("WikiPathways","NetPath","KEGG")){

    pway_db <- tolower(c)

    expect_gte(
      NROW(
        oncoEnrichR:::get_universal_enrichment(
          query_entrez = as.character(
            myc_data$entrezgene),
          genedb = oedb$genedb$all,
          logger = log4r_logger,
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
          logger = log4r_logger,
          dbsource = c,
          TERM2GENE = oedb$pathwaydb[[pway_db]]$TERM2GENE,
          TERM2NAME = oedb$pathwaydb[[pway_db]]$TERM2NAME)
      ),
      c("standard_name","description","gene_ratio","background_ratio",
        "pvalue","p.adjust","qvalue","gene_id","count","description_link",
        "exact_source",
        "external_url","db","enrichment_factor",
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
        genedb = oedb$genedb$all,
        logger = log4r_logger,
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
        logger = log4r_logger,
        dbsource = "MSIGDB/HALLMARK",
        TERM2GENE = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2GENE,
        TERM2NAME = oedb$pathwaydb$msigdb$COLLECTION$H$ALL$TERM2NAME,
        TERM2SOURCE = oedb$pathwaydb$msigdb$TERM2SOURCE)
    ),
    c("standard_name","description","gene_ratio","background_ratio",
      "pvalue","p.adjust","qvalue","gene_id","count","exact_source",
      "external_url","db","description_link","enrichment_factor",
      "gene_symbol_link","gene_symbol","setting_p_value_cutoff",
      "setting_q_value_cutoff","setting_p_value_adj_method",
      "setting_min_geneset_size","setting_max_geneset_size")
  )

  expect_identical(
    oncoEnrichR:::get_universal_enrichment(
      query_entrez = as.character(
        myc_data$entrezgene),
      genedb = oedb$genedb$all,
      background_entrez = as.character(background_sample_set$entrezgene),
      logger = log4r_logger,
      dbsource = "NetPath",
      TERM2GENE = oedb$pathwaydb$netpath$TERM2GENE,
      TERM2NAME = oedb$pathwaydb$netpath$TERM2NAME),
    NULL
  )

})
