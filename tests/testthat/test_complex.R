library(magrittr)

# load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))
#
# log4r_logger <- log4r::logger(
#   threshold = "WARN",
#   appenders = log4r::console_appender(oncoEnrichR:::log4r_layout))

test_that("Protein complex annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_protein_complex())
  expect_error(oncoEnrichR:::annotate_protein_complex(logger = log4r_logger))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    complex_db = oedb$genedb$proteincomplexdb$db,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    complex_db = oedb$genedb$proteincomplexdb$db,
    complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))

  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    complex_db = oedb$genedb$proteincomplexdb$db,
    complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all,
    query_entrez = "1042"))

  expect_error(oncoEnrichR:::annotate_protein_complex(
    logger = log4r_logger,
    complex_db = oedb$genedb$proteincomplexdb$db,
    complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
    transcript_xref = oedb$genedb$all,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))

  ## Get protein complexes related to EGFR (entrezgene = 1956)
  ## check that output is list
  expect_identical(
    typeof(
      oncoEnrichR:::annotate_protein_complex(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        complex_db = oedb$genedb$proteincomplexdb$db,
        complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)
    ),
    "list"
  )

  expect_identical(
    colnames(
      oncoEnrichR:::annotate_protein_complex(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        complex_db = oedb$genedb$proteincomplexdb$db,
        complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$omnipath
    ),
    c("complex_name", "target_genes", "literature", "complex_genes",
      "annotation_source", "disease_comment", "complex_comment",
      "confidence", "purification_method")
  )

  ## check that EGFR returns more than one hit from the OmniPath resource
  expect_gt(
    NROW(
      oncoEnrichR:::annotate_protein_complex(
        query_entrez = as.integer(1956),
        genedb = oedb$genedb$all,
        complex_db = oedb$genedb$proteincomplexdb$db,
        complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger)$omnipath
    ),
    0
  )

  expect_gt(
    NROW(
      oncoEnrichR:::annotate_protein_complex(
        query_entrez = as.integer(29844),
        genedb = oedb$genedb$all,
        complex_db = oedb$genedb$proteincomplexdb$db,
        complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
        transcript_xref = oedb$genedb$transcript_xref,
        logger = log4r_logger
      )$humap2
    ),
    as.integer(0)
  )

})
