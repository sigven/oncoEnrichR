
test_that("Gene fitness scores (Project Score) - testing ", {
  expect_error(oncoEnrichR:::get_fitness_lof_scores(
    qgenes = c("EGFR","KRAS")))

  expect_identical(
    names(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        genedb = oedb$genedb$all,
        cellmodeldb = oedb$cellmodeldb
      )
    ),
    c("targets",
      "n_targets",
      "n_essential",
      "common_essential")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        genedb = oedb$genedb$all,
        cellmodeldb = oedb$cellmodeldb
      )$targets
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        genedb = oedb$genedb$all,
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    c("symbol",
      "symbol_link_ps",
      "essential_gene",
      "model_name",
      "scaled_BF",
      "tissue",
      "model_link_ps",
      "cancer_type",
      "sample_site",
      "tissue_status",
      "n_models_dependent_tissue",
      "n_models_dependent_total")
  )

  expect_gt(
    NROW(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        genedb = oedb$genedb$all,
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("AQP4"),
        genedb = oedb$genedb$all,
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    as.integer(0)
  )

})


test_that("Target priority scores (Project Score) - testing ", {
  expect_error(oncoEnrichR:::get_target_priority_scores(
    qgenes = c("EGFR","KRAS")))

  expect_identical(
    names(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        cellmodeldb = oedb$cellmodeldb
      )
    ),
    c("targets","n_pri_targets")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        cellmodeldb = oedb$cellmodeldb
      )$targets
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    c("symbol","tumor_type",
      "priority_score")
  )

  expect_gt(
    NROW(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("AQP4"),
        cellmodeldb = oedb$cellmodeldb
      )$targets
    ),
    as.integer(0)
  )

})


