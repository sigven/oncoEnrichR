library(magrittr)

test_that("Gene fitness scores (Project Score) - testing ", {
  expect_error(oncoEnrichR:::get_fitness_lof_scores())
  expect_error(oncoEnrichR:::get_fitness_lof_scores(logger = log4r_logger))
  expect_error(oncoEnrichR:::get_fitness_lof_scores(
    qgenes = c("EGFR","KRAS"),
    logger = log4r_logger))

  expect_identical(
    names(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )
    ),
    c("targets","n_targets")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    c("symbol", "symbol_link_ps",
      "tissue", "n_gene_tissue",
      "cell_lines",
      "cmp_link", "n_gene")
  )

  expect_gt(
    NROW(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_fitness_lof_scores(
        qgenes = c("AQP4"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    as.integer(0)
  )

})


test_that("Target priority scores (Project Score) - testing ", {
  expect_error(oncoEnrichR:::get_target_priority_scores())
  expect_error(oncoEnrichR:::get_target_priority_scores(logger = log4r_logger))
  expect_error(oncoEnrichR:::get_target_priority_scores(
    qgenes = c("EGFR","KRAS"),
    logger = log4r_logger))

  expect_identical(
    names(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )
    ),
    c("targets","n_pri_targets")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    c("symbol","tumor_type",
      "priority_score")
  )

  expect_gt(
    NROW(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("EGFR","KRAS"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_target_priority_scores(
        qgenes = c("AQP4"),
        projectscoredb = oedb$projectscoredb,
        logger = log4r_logger
      )$targets
    ),
    as.integer(0)
  )

})


