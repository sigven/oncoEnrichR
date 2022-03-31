library(magrittr)

load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))

# log4r_logger <- log4r::logger(
#   threshold = "WARN",
#   appenders = log4r::console_appender(oncoEnrichR:::log4r_layout))

test_that("Ligand-receptor interactions - testing ", {
  expect_error(oncoEnrichR:::annotate_ligand_receptor_interactions())
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      logger = log4r_logger
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF"),
      logger = log4r_logger
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF"),
      genedb = oedb$genedb$all,
      logger = log4r_logger
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF"),
      genedb = oedb$genedb$all,
      ligand_receptor_db = oedb$ligandreceptordb$db,
      logger = log4r_logger
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = as.integer(c(200,300)),
      genedb = oedb$genedb$all,
      ligand_receptor_db = oedb$ligandreceptordb$db,
      ligand_receptor_xref = oedb$ligandreceptordb$xref,
      logger = log4r_logger)
  )

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_ligand_receptor_interactions(
        qgenes = c("EGFR", "EGF"),
        genedb = oedb$genedb$all,
        ligand_receptor_db = oedb$ligandreceptordb$db,
        ligand_receptor_xref = oedb$ligandreceptordb$xref,
        logger = log4r_logger)$secreted_signaling
      ),
    1
  )

})
