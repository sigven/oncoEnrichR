
test_that("Query gene validation routine - testing ", {

  expect_error(oncoEnrichR:::validate_query_genes(
    qgenes = c("EGFR","KRAS")))
  expect_error(oncoEnrichR:::validate_query_genes(
    genedb = oedb$genedb$all,
    qgenes = c("EGFR","KRAS")))
  expect_error(oncoEnrichR:::validate_query_genes(
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all,
    q_id_type = "UNKNOWN",
    qgenes = c("EGFR","KRAS")))

  expect_identical(
    names(
      oncoEnrichR:::validate_query_genes(
        qgenes = c("EGFR","KRAS"),
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "symbol"
      )
    ),
    c("found", "not_found", "all", "match_status")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::validate_query_genes(
        qgenes = c("EGFR","KRAS"),
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "symbol"
      )$all
    ),
    c("query_id", "status", "symbol", "genename")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::validate_query_genes(
        qgenes = c("EGFR","KRAS"),
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "symbol"
      )$not_found
    ),
    c("entrezgene", "name",
      "symbol")
  )

  expect_identical(
    colnames(
      oncoEnrichR:::validate_query_genes(
        qgenes = c("EGFR","KRAS"),
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "symbol"
      )$found
    ),
    c("entrezgene", "name",
      "symbol", "alias", "status", "genename")
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "NP_005219",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "refseq_protein"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "NM_005228",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "refseq_transcript_id"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "1956",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "entrezgene"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "ENST00000275493",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "ensembl_mrna"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "ENSP00000275493",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "ensembl_protein"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "ENSG00000146648",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "ensembl_gene"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::validate_query_genes(
        qgenes = "P00533",
        transcript_xref = oedb$genedb$transcript_xref,
        genedb = oedb$genedb$all,
        q_id_type = "uniprot_acc"
      )$found
    ),
    as.integer(1)
  )

  expect_identical(
    NROW(
      dplyr::filter(
        oncoEnrichR:::validate_query_genes(
          qgenes = "NISBD2",
          transcript_xref = oedb$genedb$transcript_xref,
          genedb = oedb$genedb$all,
          q_id_type = "symbol"
        )$found,
        status == "found_as_alias")
    ),
    as.integer(1)
  )

  expect_output(
    oncoEnrichR:::validate_query_genes(
      qgenes = "NISBD2",
      transcript_xref = oedb$genedb$transcript_xref,
      genedb = oedb$genedb$all,
      q_id_type = "symbol"
    ), regexp = "WARNING: target gene identifiers NOT found as primary symbols"
  )

  expect_output(
    oncoEnrichR:::validate_query_genes(
      qgenes = c("NISBD2","UNKNOWN"),
      ignore_id_err = F,
      transcript_xref = oedb$genedb$transcript_xref,
      genedb = oedb$genedb$all,
      q_id_type = "symbol"
    ), regexp = "ERROR: target gene identifiers NOT found: UNKNOWN"
  )

  expect_output(
    oncoEnrichR:::validate_query_genes(
      qgenes = c("UNKNOWN"),
      ignore_id_err = F,
      transcript_xref = oedb$genedb$transcript_xref,
      genedb = oedb$genedb$all,
      q_id_type = "symbol"
    ), regexp = "ERROR: NO target gene identifiers found"
  )

})

