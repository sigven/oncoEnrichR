
liver_genes <- c("ALB", "CYP3A4", "CYP2E1", "APOA1", "APOB",
                 "FGA", "FGB", "FGG", "HP", "TTR", "SERPINA1", "TF")

test_that("Tissue/cell-type expression specificity (HPA, tau index) - testing ", {

  expect_error(oncoEnrichR:::gene_tissue_celltype_specificity(
    qgenes = c("ALB","EGFR")))

  expect_error(oncoEnrichR:::gene_tissue_celltype_specificity(
    qgenes = c("ALB","EGFR"),
    genedb = oedb$genedb$all,
    tissuecelldb = oedb$tissuecelldb$tissue,
    resolution = "UNKNOWN"))

  expect_identical(
    names(
      oncoEnrichR:::gene_tissue_celltype_specificity(
        qgenes = c("ALB","EGFR"),
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )
    ),
    c("per_gene", "tau_distribution")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::gene_tissue_celltype_specificity(
        qgenes = c("ALB","EGFR"),
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_gene
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::gene_tissue_celltype_specificity(
        qgenes = c("ALB","EGFR"),
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_gene
    ),
    c("symbol",
      "genename",
      "cancer_max_rank",
      "tau",
      "top_context",
      "top_contexts",
      "top_expressions")
  )

  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_celltype_specificity(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_gene
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::gene_tissue_celltype_specificity(
        qgenes = c("NOTAREALSYMBOL123"),
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_gene
    ),
    as.integer(0)
  )

  ## ALB is a canonical liver-specific gene - regression guard for the
  ## tau calculation and context lookup
  alb_tissue <- oncoEnrichR:::gene_tissue_celltype_specificity(
    qgenes = c("ALB"),
    genedb = oedb$genedb$all,
    tissuecelldb = oedb$tissuecelldb$tissue,
    resolution = "tissue"
  )$per_gene

  expect_equal(alb_tissue$top_context, "liver")
  expect_gt(alb_tissue$tau, 0.9)

  alb_celltype <- oncoEnrichR:::gene_tissue_celltype_specificity(
    qgenes = c("ALB"),
    genedb = oedb$genedb$all,
    tissuecelldb = oedb$tissuecelldb$celltype,
    resolution = "celltype"
  )$per_gene

  expect_equal(alb_celltype$top_context, "hepatocytes")

})


test_that("Tissue/cell-type enrichment (HPA, hypergeometric test) - testing ", {

  expect_error(oncoEnrichR:::gene_tissue_celltype_enrichment(
    qgenes = liver_genes))

  expect_error(oncoEnrichR:::gene_tissue_celltype_enrichment(
    qgenes = liver_genes,
    genedb = oedb$genedb$all,
    tissuecelldb = oedb$tissuecelldb$tissue,
    resolution = "UNKNOWN"))

  expect_identical(
    names(
      oncoEnrichR:::gene_tissue_celltype_enrichment(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )
    ),
    c("per_gene", "per_type")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::gene_tissue_celltype_enrichment(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_type
    )
  )

  expect_identical(
    names(
      oncoEnrichR:::gene_tissue_celltype_enrichment(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_type
    ),
    c("tissue",
      "n_query_genes",
      "n_background_genes",
      "fold_enrichment",
      "p_adjusted",
      "p_value",
      "genes_contributing")
  )

  expect_identical(
    names(
      oncoEnrichR:::gene_tissue_celltype_enrichment(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$celltype,
        resolution = "celltype"
      )$per_type
    ),
    c("cell_type",
      "n_query_genes",
      "n_background_genes",
      "fold_enrichment",
      "p_adjusted",
      "p_value",
      "genes_contributing")
  )

  ## every context should be represented in per_type, regardless of
  ## whether any query genes are highly expressed there
  expect_equal(
    NROW(
      oncoEnrichR:::gene_tissue_celltype_enrichment(
        qgenes = liver_genes,
        genedb = oedb$genedb$all,
        tissuecelldb = oedb$tissuecelldb$tissue,
        resolution = "tissue"
      )$per_type
    ),
    NROW(oedb$tissuecelldb$tissue$enrichment_background_summary)
  )

  ## a well-established liver gene panel should come out as significantly
  ## enriched for liver tissue/hepatocytes - regression guard for the
  ## hypergeometric test
  liver_tissue_enrichment <-
    oncoEnrichR:::gene_tissue_celltype_enrichment(
      qgenes = liver_genes,
      genedb = oedb$genedb$all,
      tissuecelldb = oedb$tissuecelldb$tissue,
      resolution = "tissue"
    )$per_type |>
    dplyr::arrange(p_value)

  expect_equal(liver_tissue_enrichment$tissue[1], "liver")
  expect_lt(liver_tissue_enrichment$p_adjusted[1], 0.001)
  expect_gt(liver_tissue_enrichment$fold_enrichment[1], 1)

  liver_celltype_enrichment <-
    oncoEnrichR:::gene_tissue_celltype_enrichment(
      qgenes = liver_genes,
      genedb = oedb$genedb$all,
      tissuecelldb = oedb$tissuecelldb$celltype,
      resolution = "celltype"
    )$per_type |>
    dplyr::arrange(p_value)

  expect_equal(liver_celltype_enrichment$cell_type[1], "hepatocytes")
  expect_lt(liver_celltype_enrichment$p_adjusted[1], 0.001)

  ## restricting the background should shrink (or keep equal) the
  ## per-context background gene counts
  liver_tissue_enrichment_bg <-
    oncoEnrichR:::gene_tissue_celltype_enrichment(
      qgenes = liver_genes,
      genedb = oedb$genedb$all,
      tissuecelldb = oedb$tissuecelldb$tissue,
      resolution = "tissue",
      background_entrez = background_sample_entrez
    )$per_type

  expect_true(
    all(
      liver_tissue_enrichment_bg$n_background_genes <=
        liver_tissue_enrichment$n_background_genes[
          match(liver_tissue_enrichment_bg$tissue,
               liver_tissue_enrichment$tissue)]
    )
  )

})
