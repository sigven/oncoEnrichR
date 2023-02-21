
test_that("Prognostic associations - HPA - input check", {
  expect_error(oncoEnrichR:::hpa_prognostic_genes(
    qgenes = c(1042),
    genedb = oedb$genedb$all,
    hpadb = NULL))
  expect_error(oncoEnrichR:::hpa_prognostic_genes(
    qgenes = c(1042),
    genedb = oedb$genedb$all,
    hpadb = oedb$hpa,
    q_id_type = "WRONG_ID_TYPE"))
  expect_error(oncoEnrichR:::hpa_prognostic_genes(
    q_id_type = "entrezgene",
    genedb = oedb$genedb$all,
    hpadb = oedb$hpa,
    qgenes = c("200","300")))
  expect_error(oncoEnrichR:::hpa_prognostic_genes(
    q_id_type = "symbol", qgenes = c(300,400)))

  expect_gte(
    NROW(
      oncoEnrichR:::hpa_prognostic_genes(
        q_id_type = "entrezgene",
        qgenes = as.integer(3845), #KRAS
        hpadb = oedb$hpa,
        genedb = oedb$genedb$all)
      ),
    as.integer(1)
  )

  expect_equal(NROW(oncoEnrichR:::hpa_prognostic_genes(
    hpadb = oedb$hpa,
    genedb = oedb$genedb$all,
    qgenes = c("UNKNOWN_GENE"))), 0)

  expect_named(oncoEnrichR:::hpa_prognostic_genes(
    hpadb = oedb$hpa,
    genedb = oedb$genedb$all,
    qgenes = c("BRAF","KRAS")),
    c("symbol","primary_site","evidence_direction",
      "tumor_types", "p_value",
      "percentile_rank_site","percentile_rank_all",
      "log10_p_value"))


})

test_that("Prognostic associations - CSHL - input check", {
  expect_error(
    oncoEnrichR:::km_cshl_survival_genes(
      qgenes = c("BRAF"),
      projectsurvivaldb = NULL)
  )

  expect_error(
    oncoEnrichR:::km_cshl_survival_genes(
      qgenes = c(300,400)
      )
    )
  expect_named(
    oncoEnrichR:::km_cshl_survival_genes(
      survivaldb = oedb$survivaldb$cna,
      qgenes = c("BRAF","KRAS")),
    c("symbol","tcga_cohort","z_score")
    )
  expect_gt(
    NROW(
      oncoEnrichR:::km_cshl_survival_genes(
        survivaldb = oedb$survivaldb$cna,
        qgenes = c("BRAF","KRAS"))),
    0)

})

