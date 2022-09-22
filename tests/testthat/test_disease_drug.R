
test_that("Target drug associations", {
  expect_error(oncoEnrichR:::target_drug_associations())
  expect_error(oncoEnrichR:::target_drug_associations(
    qgenes = c("BRAF","EGFR")
  ))

  expect_error(oncoEnrichR:::target_drug_associations(
    genedb = oedb$genedb$all,
    qgenes = c(1042,1044)
  ))

  expect_identical(
    names(oncoEnrichR:::target_drug_associations(
    qgenes = c("EGFR","BRAF"),
    genedb = oedb$genedb$all)),
    c("target_drugs","tractability_ab","tractability_sm")
  )

  expect_identical(
    NROW(oncoEnrichR:::target_drug_associations(
      qgenes = c("MITF"),
      genedb = oedb$genedb$all)$target_drugs),
    as.integer(0)
  )

  expect_identical(
    names(oncoEnrichR:::target_drug_associations(
      qgenes = c("EGFR","BRAF"),
      genedb = oedb$genedb$all)),
    c("target_drugs","tractability_ab","tractability_sm")
  )

  expect_identical(
    colnames(oncoEnrichR:::target_drug_associations(
      qgenes = c("EGFR","BRAF"),
      genedb = oedb$genedb$all)$target_drugs),
    c("symbol","genename",
      "drugs_late_phase",
      "drugs_early_phase",
      "approved_drugs")
  )

  expect_identical(
    colnames(oncoEnrichR:::target_drug_associations(
      qgenes = c("EGFR","BRAF"),
      genedb = oedb$genedb$all)$tractability_ab),
    c("symbol",
      "AB_tractability_category",
      "AB_tractability_support")
  )

  expect_identical(
    colnames(oncoEnrichR:::target_drug_associations(
      qgenes = c("EGFR","BRAF"),
      genedb = oedb$genedb$all)$tractability_sm),
    c("symbol",
      "SM_tractability_category",
      "SM_tractability_support")
  )

})

test_that("TEST - Genes of unknown function", {
  expect_error(oncoEnrichR:::get_genes_unknown_function())
  expect_error(oncoEnrichR:::get_genes_unknown_function(
    qgenes = 'C1orf43',
  ))

  expect_identical(
    NROW(oncoEnrichR:::get_genes_unknown_function(
      qgenes = 'C1orf43',
      genedb = oedb$genedb$all)),
    as.integer(1)
  )
  expect_identical(
    NROW(oncoEnrichR:::get_genes_unknown_function(
      qgenes = 'BRAF',
      genedb = oedb$genedb$all)),
    as.integer(0)
  )
  expect_identical(
    colnames(
      oncoEnrichR:::get_genes_unknown_function(
        qgenes = 'C1orf43',
        genedb = oedb$genedb$all)
    ),
    c("symbol","genename","num_go_terms",
      "gene_summary","unknown_function_rank",
      "has_gene_summary")
  )
})


test_that("TEST - Target disease associations", {
  expect_error(oncoEnrichR:::target_disease_associations())
  expect_error(oncoEnrichR:::target_disease_associations(
    qgenes = c("EGFR","BRAF")))
  expect_error(oncoEnrichR:::target_disease_associations(
    qgenes = c("EGFR","BRAF"),
    otdb_all = oedb$otdb$all))
  expect_error(oncoEnrichR:::target_disease_associations(
    qgenes = c("EGFR","BRAF"),
    otdb_all = oedb$otdb$all,
    otdb_gene_rank = oedb$otdb$gene_rank))
  expect_error(oncoEnrichR:::target_disease_associations(
    qgenes = c(300,400),
    otdb_all = oedb$otdb$all,
    otdb_gene_rank = oedb$otdb$gene_rank,
    genedb = oedb$genedb$all))

  expect_identical(
    names(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)
    ),
    c("target", "assoc_pr_gene", "ttype_matrix")
  )

  expect_gt(
    nchar(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        show_top_diseases_only = F,
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$target$cancer_associations
    ),
    nchar(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        show_top_diseases_only = T,
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$target$cancer_associations
    )
  )

  expect_identical(
    NROW(
      oncoEnrichR:::target_disease_associations(
        qgenes = "C1orf43",
        show_top_diseases_only = T,
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$assoc_pr_gene$cancer
      ),
    as.integer(0)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::target_disease_associations(
        qgenes = "C1orf43",
        show_top_diseases_only = T,
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$assoc_pr_gene$other
    ),
    as.integer(0)
  )

  expect_identical(
    colnames(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$target
    ),
    c('symbol',
      'genename',
      'ensembl_gene_id',
      'oncogene',
      'tumor_suppressor',
      'cancer_driver',
      'targetset_cancer_rank',
      'global_cancer_rank',
      'cancergene_evidence',
      'cancer_associations',
      'cancer_association_links',
      'targetset_disease_rank',
      'disease_associations',
      'disease_association_links',
      'gene_summary')
  )

  expect_identical(
    typeof(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$assoc_pr_gene
    ),
    "list"
  )


  expect_identical(
    colnames(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$assoc_pr_gene$other
    ),
    c('symbol', 'targetset_disease_rank', 'ot_links','ot_diseases')
  )

  expect_identical(
    colnames(
      oncoEnrichR:::target_disease_associations(
        qgenes = "BRAF",
        genedb = oedb$genedb$all,
        otdb_all = oedb$otdb$all,
        otdb_gene_rank = oedb$otdb$gene_rank)$assoc_pr_gene$cancer
    ),
    c('symbol', 'targetset_cancer_rank', 'global_cancer_rank',
      'ot_cancer_links','ot_cancer_diseases')
  )

})

