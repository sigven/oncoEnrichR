
test_that("TF-target annotations - testing ", {
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC")))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    genedb = oedb$genedb$all))
  #expect_error(oncoEnrichR:::annotate_tf_targets(
  #  qgenes = c("MYC"),
  #  tf_target_interactions = oedb$tftargetdb,
  #  genedb = oedb$genedb$all))

  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = as.integer(200),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb))

  expect_gt(NROW(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC","TP53"),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb)), 0)

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_resources = 1,
        tf_target_interactions  = oedb$tftargetdb)
      ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_resources = 2,
        tf_target_interactions  = oedb$tftargetdb)
    )
  )

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_resources = 2,
        tf_target_interactions  = oedb$tftargetdb)
    ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        regulatory_min_resources = 3,
        tf_target_interactions  = oedb$tftargetdb)
    )
  )

  expect_gt(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("SMARCA2"),
        genedb = oedb$genedb$all,
        tf_target_interactions  = oedb$tftargetdb)
      ),
    0)

  expect_identical(
    typeof(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("SMARCA2"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb)
      )),
    "list"
  )

  expect_identical(
    typeof(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("MYC","TP53"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb)
      )),
    "list"
  )

  expect_identical(
    names(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("MYC","TP53"),
          genedb = oedb$genedb$all,
          tf_target_interactions = oedb$tftargetdb)
      )),
    c("nodes", "edges")
  )

})
