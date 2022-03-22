library(magrittr)

# load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))
#
# log4r_logger <- log4r::logger(
#   threshold = "WARN",
#   appenders = log4r::console_appender(oncoEnrichR:::log4r_layout))

test_that("Protein complex annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_subcellular_compartments())
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(logger = log4r_logger))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    logger = log4r_logger,
    query_entrez = as.integer(c(300,400))))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    logger = log4r_logger,
    query_entrez = as.integer(c(300,400)),
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    logger = log4r_logger,
    query_entrez = as.integer(c(300,400)),
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    logger = log4r_logger,
    query_entrez = as.integer(c(300,400)),
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all,
    comppidb = oedb$subcelldb$comppidb))

  expect_error(oncoEnrichR:::annotate_subcellular_compartments(
    query_entrez = c("200","300"),
    genedb = oedb$genedb$all,
    comppidb = oedb$subcelldb$comppidb,
    go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
    transcript_xref = oedb$genedb$transcript_xref,
    logger = log4r_logger))

  expect_identical(
    typeof(
      oncoEnrichR:::annotate_subcellular_compartments(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        comppidb = oedb$subcelldb$comppidb,
        go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)
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
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)
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
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$anatogram
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
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$grouped
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
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$all
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
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$all
    ) -
      NROW(
        oncoEnrichR:::annotate_subcellular_compartments(
          query_entrez = as.integer(1956),
          minimum_confidence = 2,
          genedb = oedb$genedb$all,
          comppidb = oedb$subcelldb$comppidb,
          go_gganatogram_map = oedb$subcelldb$go_gganatogram_map,
          transcript_xref = oedb$genedb$transcript_xref,
          logger = log4r_logger)$all
      )),
    0
  )

})
