
test_that("Subcellular compartment annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400))))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all,
    compartments = oedb$subcelldb$compartments))

  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = c("200","300"),
    genedb = oedb$genedb$all,
    compartments = oedb$subcelldb$compartments,
    go_gganatogram_map = oedb$subcelldb$go_gganatogram_map))

  expect_identical(
    typeof(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)
      ),
    "list"
  )

  expect_identical(
    names(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)
    ),
    c("all","grouped","anatogram")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)$anatogram
    ),
    c("organ","type","colour","value")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)$grouped
    ),
    c("compartment","targets","targetlinks","n")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)$all
    ),
    c("symbol","genename","compartment","minimum_confidence",
      "supporting_channels","supporting_channels_confidence",
      "supporting_sources",
      "n_supporting_channels",
      "ggcompartment")
  )


  expect_gt(
    (NROW(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        compartments = oedb$subcelldb$compartments,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)$all
    ) -
      NROW(
        oncoEnrichR:::annotate_subcellular_compartments(
          query_entrez = as.integer(1956),
          compartments_min_confidence = 4,
          genedb = oedb$genedb$all,
          compartments = oedb$subcelldb$compartments,
          go_gganatogram_map = oedb$subcelldb$go_gganatogram_map)$all
      )),
    0
  )

})
