

ppi_settings <- list()

ppi_settings[["string"]] <- list()
ppi_settings[["string"]][["minimum_score"]] <- 0.9
ppi_settings[["string"]][["visnetwork_shape"]] <- "dot"
ppi_settings[["string"]][["visnetwork_shadow"]] <- TRUE
ppi_settings[["string"]][["show_drugs"]] <- TRUE
ppi_settings[["string"]][["add_nodes"]] <- 30
ppi_settings[["string"]][["query_type"]] <- "network"
ppi_settings[["string"]][["network_type"]] <- "physical"

ppi_settings[["biogrid"]] <- list()
ppi_settings[["biogrid"]][['minimum_evidence']] <- 2
ppi_settings[["biogrid"]][["visnetwork_shape"]] <- "dot"
ppi_settings[["biogrid"]][["visnetwork_shadow"]] <- TRUE
ppi_settings[["biogrid"]][["show_drugs"]] <- TRUE
ppi_settings[["biogrid"]][["add_nodes"]] <- 30

pc_genes <-
  oedb$genedb$all[oedb$genedb$all$gene_biotype == "protein-coding",]

test_that("Protein-protein interaction network ", {
  expect_error(oncoEnrichR:::get_ppi_network(
    settings = ppi_settings[['string']]))

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

  ppi_settings$minimum_score <- 0.7

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
        qgenes = as.integer(
          dplyr::sample_n(pc_genes, 300)$entrezgene),
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
