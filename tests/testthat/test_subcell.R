
test_that("Protein complex annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400))))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = as.integer(c(300,400)),
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all,
    comppidb = oedb$subcelldb$comppidb))

  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = c("200","300"),
    genedb = oedb$genedb$all,
    comppidb = oedb$subcelldb$comppidb,
    go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
    transcript_xref = oedb$genedb$transcript_xref))

  expect_identical(
    typeof(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)
      ),
    "list"
  )

  expect_identical(
    names(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)
    ),
    c("all","grouped","anatogram")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)$anatogram
    ),
    c("organ","type","colour","value")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)$grouped
    ),
    c("compartment","targets","targetlinks","n")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)$all
    ),
    c("symbol","genename","compartment","uniprot_acc",
      "annotation_source","annotation_type","confidence",
      "ggcompartment")
  )


  expect_gt(
    (NROW(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref)$all
    ) -
      NROW(
        oncoEnrichR:::annotate_subcellular_compartments(
          query_entrez = as.integer(1956),
          minimum_confidence = 2,
          genedb = oedb$genedb$all,
          comppidb = oedb$subcelldb$comppidb,
          go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
          transcript_xref = oedb$genedb$transcript_xref)$all
      )),
    0
  )

})
