
test_that("Protein complex annotations - testing input parameters ", {
  expect_error(oncoEnrichR:::annotate_protein_complex())
  expect_error(oncoEnrichR:::annotate_protein_complex(
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    complex_db = oedb$genedb$proteincomplexdb$db,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))
  expect_error(oncoEnrichR:::annotate_protein_complex(
    complex_db = oedb$genedb$proteincomplexdb$db,
    complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
    genedb = oedb$genedb$all,
    query_entrez = as.integer(1042)))

  expect_error(oncoEnrichR:::annotate_protein_complex(
    complex_db = oedb$genedb$proteincomplexdb$db,
    complex_up_xref = oedb$genedb$proteincomplexdb$up_xref,
    transcript_xref = oedb$genedb$transcript_xref,
    genedb = oedb$genedb$all,
    query_entrez = "1042"))

  expect_error(oncoEnrichR:::annotate_protein_complex(
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
        otdb_gene_rank = oedb$otdb$gene_rank)
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
        otdb_gene_rank = oedb$otdb$gene_rank)$omnipath
    ),
    c("complex_name", "target_genes",
      "literature", "complex_genes",
      "annotation_source", "disease_comment",
      "complex_cancer_rank_score", "num_target_members",
      "complex_comment",
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
        otdb_gene_rank = oedb$otdb$gene_rank)$omnipath
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
        otdb_gene_rank = oedb$otdb$gene_rank
      )$humap2
    ),
    as.integer(0)
  )

})
