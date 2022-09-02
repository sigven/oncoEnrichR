
ppi_settings <- list()
ppi_settings[["minimum_score"]] <- 900
ppi_settings[["visnetwork_shape"]] <- "dot"
ppi_settings[["visnetwork_shadow"]] <- T
ppi_settings[["show_drugs"]] <- T
ppi_settings[["add_nodes"]] <- 50
ppi_settings[["query_type"]] <- "network"

pc_genes <-
  oedb$genedb$all[oedb$genedb$all$gene_biotype == "protein-coding",]

test_that("Tissue/celltype enrichment categories ", {
  expect_error(oncoEnrichR:::get_ppi_network(
    settings = ppi_settings
  ))
  expect_error(oncoEnrichR:::get_ppi_network(
    settings = ppi_settings,
    genedb = oedb$genedb$all
  ))
  expect_error(oncoEnrichR:::get_ppi_network(
    settings = ppi_settings,
    genedb = oedb$genedb$all,
    cancerdrugdb = oedb$cancerdrugdb
  ))
  expect_error(oncoEnrichR:::get_ppi_network(
    settings = ppi_settings,
    genedb = oedb$genedb$all,
    cancerdrugdb = oedb$cancerdrugdb,
    qgenes = c("BRAF")
  ))

  expect_identical(
    names(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(25912),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )
    ),
    c("source","complete_network",
      "hubscores")
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(25912),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$complete_network$edges
    ),
    as.integer(0)
  )

  expect_equal(
    NROW(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(25912),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$complete_network$nodes
    ),
    as.integer(0)
  )

  ppi_settings$minimum_score <- 700

  expect_identical(
    names(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
        )
      ),
    c("source","complete_network",
      "community_network","hubscores")
  )

  expect_identical(
    names(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(dplyr::sample_n(pc_genes, 300)$entrezgene),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )
    ),
    c("source","complete_network",
      "community_network","hubscores")
  )

  expect_identical(
    names(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$complete_network
    ),
    c("nodes",
      "edges")
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$complete_network$nodes
    )
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$complete_network$edges
    )
  )

  expect_true(
    is.data.frame(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$hubscores
    )
  )

  expect_identical(
    colnames(
      oncoEnrichR:::get_ppi_network(
        qgenes = as.integer(c(1956, 673)),
        genedb = oedb$genedb$all,
        settings = ppi_settings,
        cancerdrugdb = oedb$cancerdrugdb
      )$hubscores
    ),
    c("symbol","name","hub_score")
  )

})
