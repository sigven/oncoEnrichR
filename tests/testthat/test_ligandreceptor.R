
# load(system.file("internal_db/oedb.rda", package = "oncoEnrichR"))
#
test_that("Ligand-receptor interactions - testing ", {
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF")
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF"),
      genedb = oedb$genedb$all
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = c("EGFR", "EGF"),
      genedb = oedb$genedb$all,
      ligand_receptor_db = oedb$ligandreceptordb$db
    )
  )
  expect_error(
    oncoEnrichR:::annotate_ligand_receptor_interactions(
      qgenes = as.integer(c(200,300)),
      genedb = oedb$genedb$all,
      ligand_receptor_db = oedb$ligandreceptordbdb,
      ligand_receptor_xref = oedb$ligandreceptordbxref)
  )

  expect_gte(
    NROW(
      oncoEnrichR:::annotate_ligand_receptor_interactions(
        qgenes = c("EGFR", "EGF"),
        genedb = oedb$genedb$all,
        ligand_receptor_db = oedb$ligandreceptordb$db,
        ligand_receptor_xref = oedb$ligandreceptordb$xref)$secreted_signaling
      ),
    1
  )

})
