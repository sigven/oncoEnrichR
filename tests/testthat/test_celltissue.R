
test_that("Tissue/celltype enrichment categories ", {
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat())
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df))
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df))
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    genedb = oedb$genedb$all,
    qgenes = c("BRAF","EGFR","KRAS"),
    resolution = "UNKNOWN",
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df))
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    #logger = log4r_logger,
    genedb = oedb$genedb$all,
    qgenes = c("BRAF","EGFR","KRAS"),
    resolution = "tissue",
    q_id_type = "entrezgene",
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df))
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    qgenes = c("200","300"),
    q_id_type = "entrezgene",
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
    resolution = "tissue"))

  ## hpa enrichment df's not correctly provided
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    qgenes = as.integer(c(200, 300)),
    q_id_type = "entrezgene",
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$expr_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
    resolution = "tissue"))

  ## mismatch between resolution and hpa df's provided
  expect_error(oncoEnrichR:::gene_tissue_cell_spec_cat(
    qgenes = as.integer(c(200, 300)),
    q_id_type = "entrezgene",
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
    resolution = "single_cell"))

  expect_identical(
    NROW(
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        qgenes = as.integer(c(-1, -2)),
        q_id_type = "entrezgene",
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
        resolution = "tissue")$category_df
    ),
    as.integer(12)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        qgenes = as.integer(c(-1, -2)),
        q_id_type = "entrezgene",
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
        resolution = "tissue")$exp_dist_df
    ),
    as.integer(0)
  )

  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        qgenes = as.integer(c(1956, 673)),
        q_id_type = "entrezgene",
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_expr_db_df = oedb$tissuecelldb$tissue$expr_df,
        resolution = "tissue")$exp_dist_df
    ),
    as.integer(0)
  )

  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        qgenes = as.integer(c(1956, 673)),
        q_id_type = "entrezgene",
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_expr_db_df = oedb$tissuecelldb$single_cell$expr_df,
        resolution = "single_cell")$exp_dist_df
    ),
    as.integer(0)
  )

  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        qgenes = as.character(c("KRAS", "EGFR")),
        q_id_type = "symbol",
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_expr_db_df = oedb$tissuecelldb$single_cell$expr_df,
        resolution = "single_cell")$exp_dist_df
    ),
    as.integer(0)
  )



})

test_that("Tissue/celltype gene enrichment ", {
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment())
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment(
    genedb = oedb$genedb$all
  ))
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df
  ))
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE
  ))
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
    qgenes_entrez = as.character(c(1042))
  ))
  expect_error(oncoEnrichR:::gene_tissue_cell_enrichment(
    genedb = oedb$genedb$all,
    hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
    hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
    qgenes_entrez = as.integer(c(1042)),
    resolution = "UNKNOWN"
  ))

  expect_error(
    oncoEnrichR:::gene_tissue_cell_enrichment(
      qgenes_entrez = c("200","300"),
      genedb = oedb$genedb$all,
      hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
      hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
      resolution = "tissue")
  )

  ## hpa enrichment df's not correctly provided
  expect_error(
    oncoEnrichR:::gene_tissue_cell_enrichment(
      qgenes_entrez = as.integer(c(200, 300)),
      genedb = oedb$genedb$all,
      hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
      hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_df,
      resolution = "tissue"))

  ## mismatch between resolution and hpa df's provided
  expect_error(
    oncoEnrichR:::gene_tissue_cell_enrichment(
      qgenes_entrez = as.integer(c(200, 300)),
      genedb = oedb$genedb$all,
      hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
      hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
      resolution = "single_cell"))

  expect_identical(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(-1, -2)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
        resolution = "tissue")$per_gene),
    as.integer(0)
  )

  expect_identical(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(-1, -2)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
        resolution = "tissue")$per_type),
    as.integer(54))

  expect_identical(
    colnames(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(-1, -2)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
        resolution = "tissue")$per_type),
    c("tissue","tissue_specific_genes",
      "fold_change","log10_pvalue"))

  expect_identical(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(-1, -2)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$single_cell$te_SE,
        resolution = "single_cell")$per_type
    ),
    as.integer(79)
  )

  expect_identical(
    colnames(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(-1, -2)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$single_cell$te_SE,
        resolution = "single_cell")$per_type),
    c("cell_type","celltype_specific_genes",
      "fold_change","log10_pvalue"))


  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$tissue$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$tissue$te_SE,
        resolution = "tissue")$per_gene
    ),
    0
  )

  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$single_cell$te_SE,
        resolution = "single_cell")$per_gene
    ),
    0
  )


  expect_gt(
    NROW(
      oncoEnrichR:::gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(c(1956, 673)),
        background_entrez = as.integer(
          background_sample_entrez),
        genedb = oedb$genedb$all,
        hpa_enrichment_db_df = oedb$tissuecelldb$single_cell$te_df,
        hpa_enrichment_db_SE = oedb$tissuecelldb$single_cell$te_SE,
        resolution = "single_cell")$per_type
    ),
    0
  )

})
