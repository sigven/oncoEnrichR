
test_that("TF-target annotations - testing ", {
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC")))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    collection = "UNKNOWN_COLLECTION",
    tf_target_interactions = oedb$tftargetdb,
    genedb = oedb$genedb$all))

  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = as.integer(200),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb,
    collection = "global"))

  expect_gt(NROW(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC","TP53"),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb,
    collection = "pancancer")), 0)

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_confidence = "D",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer")
      ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_confidence = "C",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer")
    )
  )

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_confidence = "B",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer")
    ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_confidence = "A",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer")
    )
  )

  expect_gt(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("SMARCA2"),
        genedb = oedb$genedb$all,
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer")
      ),
    0)

  expect_identical(
    type(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("SMARCA2"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb,
          collection = "pancancer")
      )),
    "list"
  )

  expect_identical(
    type(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("MYC","TP53"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb,
          collection = "pancancer")
      )),
    "list"
  )

  expect_identical(
    names(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("MYC","TP53"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb,
          collection = "pancancer")
      )),
    c("nodes", "edges")
  )

})
