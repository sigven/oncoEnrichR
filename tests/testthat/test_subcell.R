
test_that("Subcellular compartment annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400))))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  # expect_error(oncoEnrichR:::annotate_subcellular_compartments(
  #   query_entrez = as.integer(c(300,400)),
  #   genedb = oedb$genedb$all,
  #   compartments = oedb$subcelldb))

  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = c("200","300"),
    genedb = oedb$genedb$all,
    compartments = oedb$subcelldb))


  expect_identical(
    typeof(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(c(1956,3845)),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)
      ),
    "list"
  )

  expect_identical(
    names(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)
    ),
    c("all","grouped","comp_density")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)$comp_density
    ),
    c("id","name","n_comp","genes","proportion","bin","fill")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)$grouped
    ),
    c("compartment","targets","targetlinks","n")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)$all
    ),
    c("symbol","genename",
      "compartment",
      "cancer_max_rank",
      "maximum_confidence",
      "supporting_channels",
      "supporting_channels_confidence",
      "supporting_sources",
      "n_supporting_channels",
      "id","name","subcellular_location_id")
  )


  expect_gt(
    (NROW(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb)$all
    ) -
      NROW(
        oncoEnrichR:::annotate_subcellular_compartments(
          query_entrez = as.integer(1956),
          compartments_min_confidence = 4,
          genedb = oedb$genedb$all,
          compartments = oedb$subcelldb)$all
      )),
    0
  )

})
