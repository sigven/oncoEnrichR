
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
     oncoEnrichR::cp_output_cols
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
      oncoEnrichR::cp_output_cols
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
    oncoEnrichR::cp_output_cols
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
