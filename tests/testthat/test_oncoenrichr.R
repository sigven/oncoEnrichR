
test_that("oncoEnrichR - initialize report", {

  expect_identical(
    typeof(
      oncoEnrichR:::init_report(
        oedb,
        project_title = "Project title",
        project_owner = "Project owner",
        html_floating_toc = T,
        html_report_theme = "default",
        query_id_type = "symbol",
        ignore_id_err = TRUE,
        project_description = "Project description",
        ppi_string_min_score = 0.9,
        ppi_add_nodes = 30,
        ppi_show_drugs = T,
        bgset_description =
          "All annotated protein-coding genes",
        bgset_id_type = "symbol",
        enrichment_p_value_cutoff = 0.05,
        enrichment_p_value_adj = "BH",
        enrichment_q_value_cutoff = 0.2,
        enrichment_min_geneset_size = 10,
        enrichment_max_geneset_size = 500,
        enrichment_plot_num_terms = 20,
        enrichment_simplify_go = F,
        subcellcomp_min_confidence = 2,
        subcellcomp_channel_types = c("Experimental","Text mining", "Knowledge"),
        subcellcomp_show_cytosol = F,
        regulatory_min_confidence = "D",
        show_ppi = T,
        show_disease = T,
        show_top_diseases_only = T,
        show_cancer_hallmarks = T,
        show_drug = T,
        show_enrichment = T,
        show_aberration = T,
        show_coexpression = T,
        show_subcell_comp = T,
        show_fitness = T,
        show_cell_tissue = F,
        show_ligand_receptor = T,
        show_regulatory = T,
        show_unknown_function = T,
        show_prognostic = T,
        show_complex = T)
    ),
    "list"
  )

})


# test_that("oncoEnrichR - load database testing", {
#
#   expect_identical(
#     typeof(
#       oncoEnrichR:::load_db(cache_dir = "~/.oncoenrichr_cache", remote = F)
#     ),
#     "list"
#   )
#
# })

test_that("oncoEnrichR - generate report", {

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data[['symbol']],
      oeDB = NULL
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = NULL,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = as.integer(100),
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      bgset = as.integer(100),
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  ## single gene not allowed
  expect_output(
    oncoEnrichR::onco_enrich(
      query = "BRAF",
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      oeDB = oedb,
      enrichment_p_value_adj = "UNKNOWN_ADJUSTMENT"
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      oeDB = oedb,
      regulatory_min_confidence = "F"
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      oeDB = oedb,
      html_report_theme = "UNKNOWN_THEME"
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      oeDB = oedb,
      galaxy = T,
      show_enrichment = F,
      show_aberration = F,
      show_coexpression = F
    ),
    regexp = "NOTE: Running oncoEnrichR workflow in Galaxy mode"
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      query_id_type = "UNKNOWN_TYPE",
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      bgset_id_type = "UNKNOWN_TYPE",
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      ppi_string_min_score = 1.2,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      ppi_biogrid_min_evidence =  11,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      enrichment_plot_num_terms = 10.7,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      enrichment_plot_num_terms = 35,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      subcellcomp_min_confidence = 6,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      enrichment_p_value_cutoff = 1.2,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      enrichment_q_value_cutoff = 1.2,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      fitness_max_score = 1,
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_output(
    oncoEnrichR::onco_enrich(
      query = c("BRAF","TESTZ4"),
      oeDB = oedb
    ),
    regexp = "ERROR: "
  )

  expect_error(
    oncoEnrichR::onco_enrich(
      query = myc_data$symbol,
      ppi_add_nodes = 70,
      oeDB = oedb
    )
  )

  expect_identical(
    typeof(
      oncoEnrichR::onco_enrich(
        query = myc_data[['symbol']],
        oeDB = oedb,
        html_floating_toc = T,
        html_report_theme = "default",
        query_id_type = "symbol",
        show_cell_tissue = T,
        show_enrichment = T,
        project_title = "cMYC_BioID_screen",
        project_owner = "Raught et al.")
    ),
    "list"
  )

  # expect_identical(
  #   typeof(
  #     oncoEnrichR::onco_enrich(
  #       query = myc_data[['symbol']],
  #       bgset = background_sample_entrez,
  #       bgset_id_type = "entrezgene",
  #       bgset_description = "Sample background set",
  #       oeDB = oedb,
  #       html_floating_toc = T,
  #       html_report_theme = "default",
  #       query_id_type = "symbol",
  #       show_tcga_aberration = F,
  #       show_tcga_coexpression = F,
  #       show_enrichment = F,
  #       project_title = "cMYC_BioID_screen",
  #       project_owner = "Raught et al.")
  #   ),
  #   "list"
  # )

})

