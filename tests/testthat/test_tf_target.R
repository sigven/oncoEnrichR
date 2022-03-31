library(magrittr)

test_that("TF-target annotations - testing ", {
  expect_error(oncoEnrichR:::annotate_tf_targets())
  expect_error(oncoEnrichR:::annotate_tf_targets(logger = log4r_logger))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    logger = log4r_logger))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    genedb = oedb$genedb$all,
    logger = log4r_logger))
  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC"),
    collection = "UNKNOWN_COLLECTION",
    tf_target_interactions = oedb$tftargetdb,
    genedb = oedb$genedb$all,
    logger = log4r_logger))

  expect_error(oncoEnrichR:::annotate_tf_targets(
    qgenes = as.integer(200),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb,
    collection = "global",
    logger = log4r_logger))

  expect_gt(NROW(oncoEnrichR:::annotate_tf_targets(
    qgenes = c("MYC","TP53"),
    genedb = oedb$genedb$all,
    tf_target_interactions  = oedb$tftargetdb,
    collection = "pancancer",
    logger = log4r_logger)), 0)

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        min_confidence_reg_interaction = "D",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer",
        logger = log4r_logger)
      ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        min_confidence_reg_interaction = "C",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer",
        logger = log4r_logger)
    )
  )

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        min_confidence_reg_interaction = "B",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer",
        logger = log4r_logger)
    ),
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("TP53"),
        genedb = oedb$genedb$all,
        min_confidence_reg_interaction = "A",
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer",
        logger = log4r_logger)
    )
  )

  expect_gt(
    NROW(
      oncoEnrichR:::annotate_tf_targets(
        qgenes = c("SMARCA2"),
        genedb = oedb$genedb$all,
        tf_target_interactions  = oedb$tftargetdb,
        collection = "pancancer",
        logger = log4r_logger)
      ),
    0)

  expect_identical(
    type(
      oncoEnrichR:::retrieve_tf_target_network(
        oncoEnrichR:::annotate_tf_targets(
          qgenes = c("SMARCA2"),
          genedb = oedb$genedb$all,
          tf_target_interactions  = oedb$tftargetdb,
          collection = "pancancer",
          logger = log4r_logger)
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
          collection = "pancancer",
          logger = log4r_logger)
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
          collection = "pancancer",
          logger = log4r_logger)
      )),
    c("nodes", "edges")
  )

})
