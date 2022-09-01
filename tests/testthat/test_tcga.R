
test_that("TCGA oncoplot generation - testing ", {
  expect_error(oncoEnrichR:::tcga_oncoplot_genes(
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_oncoplot_genes(
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_oncoplot_genes(
    qgenes = c("KRAS","MYC"),
    site = "Unknown_Site",
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))

  expect_identical(
    class(
      oncoEnrichR:::tcga_oncoplot_genes(
        qgenes = c("MYC","KRAS"),
        tcgadb = oedb$tcgadb,
        site = "Lung",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    ),
    "data.frame"
  )

  expect_identical(
    class(
      oncoEnrichR:::tcga_oncoplot_genes(
        qgenes = as.integer(c(3845, 4609)),
        tcgadb = oedb$tcgadb,
        qsource = "entrezgene",
        site = "Lung",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    ),
    "data.frame"
  )


  expect_gte(
    NROW(
      oncoEnrichR:::tcga_oncoplot_genes(
        qgenes = as.integer(3845),
        qsource = "entrezgene",
        tcgadb = oedb$tcgadb,
        site = "Lung",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    ),
    as.integer(1)
  )


  expect_identical(
    colnames(
      oncoEnrichR:::tcga_oncoplot_genes(
        qgenes = c("MYC","KRAS"),
        tcgadb = oedb$tcgadb,
        site = "Lung",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    ),
    c("symbol", "variant_type", "primary_site", "percent_mutated",
      "samples_mutated", "cohort_size", "percentile")
  )


})

test_that("TCGA aberration table - testing ", {
  expect_error(oncoEnrichR:::tcga_aberration_table(
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_aberration_table(
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_aberration_table(
    qgenes = c("KRAS","MYC"),
    vtype = "UNKNOWN_TYPE",
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))

  ## No copy number data
  expect_true(
    is.data.frame(
      oncoEnrichR:::tcga_aberration_table(
        qgenes = as.integer(c(3845, 4609)),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        genedb = oedb$genedb$all
      )
    )
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::tcga_aberration_table(
        qgenes = c("MYC","KRAS"),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        qsource = "symbol",
        genedb = oedb$genedb$all
      )
    )
  )

  expect_gte(
    NROW(
      oncoEnrichR:::tcga_aberration_table(
        qgenes = c("MYC","KRAS"),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        qsource = "symbol",
        genedb = oedb$genedb$all
      )
    ), as.integer(0)
  )


  expect_identical(
    colnames(
      oncoEnrichR:::tcga_aberration_table(
        qgenes = as.integer(c(3845, 4609)),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        genedb = oedb$genedb$all
      )
    ),
    c("gene", "variant_type", "primary_site",
      "primary_diagnosis", "percent_mutated",
      "samples_mutated","cohort_size",
      "percentile")
  )

})


test_that("TCGA co-expression data - testing ", {
  expect_error(oncoEnrichR:::tcga_coexpression(
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_coexpression(
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_coexpression(
    qgenes = c(200,300),
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))

  expect_true(
    is.data.frame(
      oncoEnrichR:::tcga_coexpression(
        tcgadb = oedb$tcgadb,
        genedb = oedb$genedb$all,
        qgenes = c("MYC","KRAS")
        )
      )
  )
  expect_gte(
    NROW(
      oncoEnrichR:::tcga_coexpression(
        tcgadb = oedb$tcgadb,
        genedb = oedb$genedb$all,
        qgenes = c("MYC","KRAS")
      )
    ), 1
  )



})


test_that("TCGA aberration matrix - testing ", {

  expect_error(oncoEnrichR:::tcga_aberration_matrix(
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_aberration_matrix(
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_aberration_matrix(
    qgenes = c("KRAS","MYC"),
    vtype = "UNKNOWN_TYPE",
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))
  expect_error(oncoEnrichR:::tcga_aberration_matrix(
    qgenes = c("KRAS","MYC"),
    vtype = "cna_homdel",
    cstrata = "UNKNOWN_STRATA",
    tcgadb = oedb$tcgadb,
    genedb = oedb$genedb$all))

  ## No copy number data
  expect_identical(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("NATP","AAVS1"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_ampl",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    NULL
  )

  expect_output(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("NATP","AAVS1"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_ampl",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    regexp = "NOTE: NO genes in query set with TCGA aberration data"
  )

  expect_identical(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = as.integer(c(3845, 4609)),
      tcgadb = oedb$tcgadb,
      qsource = "entrezgene",
      vtype = "cna_ampl",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    NULL
  )

  ## Few query genes for heatmap
  expect_identical(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("MYC","EGFR"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_ampl",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    NULL
  )

  expect_output(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("MYC","EGFR"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_ampl",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    regexp = "NOTE: Limited number of genes"
  )

  ## Few query genes for heatmap
  expect_identical(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("CARD17","BIRC8","C3orf36"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_homdel",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    NULL
  )

  expect_output(
    oncoEnrichR:::tcga_aberration_matrix(
      qgenes = c("CARD17","BIRC8","C3orf36"),
      tcgadb = oedb$tcgadb,
      vtype = "cna_homdel",
      cstrata = "site",
      genedb = oedb$genedb$all
    ),
    regexp = "NOTE: NO genes in query set with TCGA aberration data"
  )

  expect_identical(
    typeof(
      oncoEnrichR:::tcga_aberration_matrix(
        qgenes = c("MYC","EGFR","MDM2"),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    ),
    "double"
  )

  expect_true(
    is.matrix(
      oncoEnrichR:::tcga_aberration_matrix(
        qgenes = c("MYC","EGFR","MDM2"),
        tcgadb = oedb$tcgadb,
        vtype = "cna_ampl",
        cstrata = "site",
        genedb = oedb$genedb$all
      )
    )
  )


})
