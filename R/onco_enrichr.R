
init_report <- function(project_title = "Project title",
                        project_owner = "Project owner",
                        project_description = "Project description",
                        ppi_min_string_score = 900,
                        ppi_add_nodes = 50,
                        background_enrichment_description = "All protein-coding genes",
                        p_value_cutoff_enrichment = 0.05,
                        p_value_adjustment_method = "BH",
                        q_value_cutoff_enrichment = 0.2,
                        min_geneset_size = 10,
                        max_geneset_size = 500,
                        simplify_go = F,
                        show_ppi = T,
                        show_drugs_in_ppi = T,
                        show_disease = T,
                        show_drug = T,
                        show_gene_summary = F,
                        show_enrichment = T,
                        show_tcga_aberration = T,
                        show_tcga_coexpression = T,
                        show_subcell_comp = T,
                        show_crispr_lof = T,
                        #show_gtex_coexp = F,
                        show_complex = T) {

  ## report object
  rep <- list()

  ## two main elements
  # 1. data - contains all annotations and enrichment results
  # 2. config - contains all settings and underlying data sources
  for (e in c("data","config")) {
    rep[[e]] <- list()
  }

  ## release notes - software and database versions
  rep[["config"]][["resources"]] <-  oncoEnrichR::release_notes

  ## logicals indicating which sections/analyses of the report to include
  rep[["config"]][["show"]] <- list()
  rep[["config"]][["show"]][["ppi"]] <- show_ppi
  rep[["config"]][["show"]][["disease"]] <- show_disease
  rep[["config"]][["show"]][["drug"]] <- show_drug
  rep[["config"]][["show"]][["gene_summary"]] <- show_gene_summary
  rep[["config"]][["show"]][["enrichment"]] <- show_enrichment
  rep[["config"]][["show"]][["protein_complex"]] <- show_complex
  rep[["config"]][["show"]][["tcga_aberration"]] <- show_tcga_aberration
  rep[["config"]][["show"]][["tcga_coexpression"]] <- show_tcga_coexpression
  rep[["config"]][["show"]][["subcellcomp"]] <- show_subcell_comp
  rep[["config"]][["show"]][["loss_of_fitness"]] <- show_crispr_lof
  #rep[["config"]][["show"]][["gtex_coexp"]] <- show_gtex_coexp

  ## report metadata - project owner, background, project_title
  rep[["config"]][["project_title"]] <- project_title
  rep[["config"]][["project_description"]] <- project_description
  rep[["config"]][["project_owner"]] <- project_owner

  ## settings for showing disease associations
  rep[["config"]][["disease"]] <- list()
  rep[["config"]][["disease"]][["breaks"]] <- c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  rep[["config"]][["disease"]][["colors"]] <- c("#b8b8ba","#EFF3FF","#C6DBEF",
                                                "#9ECAE1","#6BAED6","#4292C6",
                                                "#2171B5","#084594")

  ## settings for crispr loss-of-fitness
  rep[["config"]][["loss_of_fitness"]] <- list()
  rep[["config"]][["loss_of_fitness"]][["plot_height"]] <- 10

  ## protein-protein interaction settings
  rep[["config"]][["ppi"]] <- list()
  rep[["config"]][["ppi"]][["stringdb"]] <- list()
  rep[["config"]][["ppi"]][["stringdb"]][["minimum_score"]] <- ppi_min_string_score
  rep[["config"]][["ppi"]][["stringdb"]][["visnetwork_shape"]] <- "dot"
  rep[["config"]][["ppi"]][["stringdb"]][["visnetwork_shadow"]] <- T
  rep[["config"]][["ppi"]][["stringdb"]][["show_drugs"]] <- show_drugs_in_ppi
  rep[["config"]][["ppi"]][["stringdb"]][["add_nodes"]] <- ppi_add_nodes
  rep[["config"]][["ppi"]][["stringdb"]][["query_type"]] <- "network"

  ## enrichment settings
  rep[["config"]][["enrichment"]] <- list()
  rep[["config"]][["enrichment"]][["p_value_cutoff"]] <- p_value_cutoff_enrichment
  rep[["config"]][["enrichment"]][["q_value_cutoff"]] <- q_value_cutoff_enrichment
  rep[["config"]][["enrichment"]][["p_adjust_method"]] <- p_value_adjustment_method
  rep[["config"]][["enrichment"]][["min_gs_size"]] <- min_geneset_size
  rep[["config"]][["enrichment"]][["max_gs_size"]] <- max_geneset_size
  rep[["config"]][["enrichment"]][["simplify_go"]] <- simplify_go
  rep[["config"]][["enrichment"]][["background_set"]] <- background_enrichment_description

  ## TCGA/GTex
  rep[["config"]][["co_expression_gtex"]] <- list()
  rep[["config"]][["co_expression_gtex"]][["plot_height"]] <- 5
  rep[["config"]][["tcga_aberration"]] <- list()
  rep[["config"]][["tcga_aberration"]][["plot_height"]] <- 14

  for (analysis in c("tcga","disease","ppi","tcga","gtex","enrichment","drug",
                    "protein_complex","subcellcomp","loss_of_fitness")) {
    rep[["data"]][[analysis]] <- list()
  }

  rep[["data"]][["drug"]][["target"]] <- data.frame()
  rep[["data"]][["disease"]][["target"]] <- data.frame()
  rep[["data"]][["disease"]][["target_assoc"]] <- data.frame()
  rep[["data"]][["disease"]][["assoc_pr_gene"]] <- list()
  rep[["data"]][["ppi"]][["complete_network"]] <- NULL
  rep[["data"]][["ppi"]][["hubscores"]] <- data.frame()
  rep[["data"]][["ppi"]][["community_network"]] <- NULL
  for (c in c("go","msigdb","wikipathwaydb","keggdb")) {
    rep[["data"]][["enrichment"]][[c]] <- data.frame()
  }

  rep[["data"]][["loss_of_fitness"]][["plot"]] <- NULL
  rep[["data"]][["loss_of_fitness"]][["df"]] <- data.frame()
  rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] <- 0

  rep[["data"]][["gtex"]][["co_expression"]] <- list()
  rep[["data"]][["gtex"]][["co_expression"]][["plots"]] <- NULL
  rep[["data"]][["gtex"]][["co_expression"]][["df"]] <- data.frame()
  rep[["data"]][["tcga"]][["co_expression"]] <- data.frame()
  rep[["data"]][["protein_complex"]][["complex"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["all"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["grouped"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["anatogram"]] <- data.frame()

  rep[["data"]][["tcga"]][["aberration"]] <- list()
  rep[["data"]][["tcga"]][["aberration"]][["table"]] <- list()
  rep[["data"]][["tcga"]][["aberration"]][["plot"]] <- list()

  for (v in c("cna_homdel","cna_ampl","snv_indel")) {
    if(v != "snv_indel"){
      rep[["data"]][["tcga"]][["aberration"]][["plot"]][[v]] <- NULL
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- data.frame()
    }
    else{
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- list()
      i <- 1
      while(i <= nrow(oncoEnrichR::maf_codes)){
        site <- oncoEnrichR::maf_codes[i,]$primary_site
        code <- oncoEnrichR::maf_codes[i,]$code
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]] <- list()
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]][['code']] <- code
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]][['top_mutated_genes']] <- data.frame()
        i <- i + 1
      }
    }
  }

  return(rep)
}

#' Function that interrogates and analyzes a list of human protein-coding genes for cancer relevance
#'
#' @param query character vector with gene/query identifiers
#' @param query_source character indicating source of query (one of "uniprot_acc", "symbol", "entrezgene", or "ensembl_gene_id")
#' @param ignore_unknown logical indicating if analysis should continue when uknown query identifiers are encountered
#' @param project_title project title (title of report)
#' @param project_owner name of project owner
#' @param project_description project background information
#' @param background_enrichment character vector with gene identifiers, used as reference/background for enrichment/over-representation analysis
#' @param background_enrichment_source character indicating source of background ("uniprot_acc","symbol","entrezgene","ensembl_gene_id")
#' @param background_enrichment_description character indicating type of background (e.g. "All lipid-binding proteins (n = 200)")
#' @param p_value_cutoff_enrichment cutoff p-value for enrichment/over-representation analysis
#' @param p_value_adjustment_method one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param q_value_cutoff_enrichment cutoff q-value for enrichment analysis
#' @param min_geneset_size minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param max_geneset_size maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param simplify_go remove highly similar GO terms in results from GO enrichment/over-representation analysis
#' @param ppi_add_nodes number of nodes to add to query set when computing the protein-protein interaction network (STRING)
#' @param ppi_score_threshold minimum score (0-1000) for retrieval of protein-protein interactions (STRING)
#' @param show_ppi logical indicating if report should contain protein-protein interaction data (STRING)
#' @param show_gene_summary logical indicating if report should fetch summary of gene function (from RefSeq/mygene.info, will increase processing time)
#' @param show_drugs_in_ppi logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform)
#' @param show_enrichment logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME etc.)
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of queryset with oncogenes/tumor suppressor genes
#' @param show_subcell_comp logical indicating if report should list subcellular compartment annotations (ComPPI)
#' @param show_crispr_lof logical indicating if report should list results from CRISPR/Cas9 loss-of-fitness screens (Project Score)
#' @param show_complex logical indicating if report should list proteins in known protein complexes (CORUM)
#' @export
#'
onco_enrich <- function(query,
                   query_source = "symbol",
                   ignore_unknown = FALSE,
                   project_title = "Project title",
                   project_owner = "Project owner",
                   project_description = "Project description",
                   background_enrichment = NULL,
                   background_enrichment_source = "symbol",
                   background_enrichment_description = "All protein-coding genes",
                   p_value_cutoff_enrichment = 0.05,
                   p_value_adjustment_method = "BH",
                   q_value_cutoff_enrichment = 0.2,
                   min_geneset_size = 10,
                   max_geneset_size = 500,
                   simplify_go = F,
                   ppi_add_nodes = 50,
                   ppi_score_threshold = 900,
                   show_gene_summary = F,
                   show_ppi = T,
                   show_drugs_in_ppi = F,
                   show_disease = T,
                   show_drug = T,
                   show_enrichment = T,
                   show_tcga_aberration = T,
                   show_tcga_coexpression = T,
                   show_subcell_comp = T,
                   show_crispr_lof = T,
                   #show_gtex_coexp = F,
                   show_complex = T) {
  #gtex_atlasassay_groups = c("g32","g9","g29","g10","g28","g44","g33","g50","g37","g38","g42","g35")) {
  stopifnot(is.character(query))
  stopifnot(p_value_adjustment_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  if (length(query) > 800 | length(query) < 20) {
    rlogging::message(paste0("ERROR: oncoEnrichR needs minimum 20 query identifiers, and accepts a maximum of 800. Query contained n = ",length(query), " identifiers"))
    return(NULL)
  }
  stopifnot(query_source == "symbol" | query_source == "entrezgene" |
              query_source == "uniprot_acc" | query_source == "ensembl_gene_id")
  stopifnot(ppi_score_threshold > 0 & ppi_score_threshold <= 1000)
  stopifnot(p_value_cutoff_enrichment > 0 & p_value_cutoff_enrichment < 1)
  stopifnot(q_value_cutoff_enrichment > 0 & q_value_cutoff_enrichment < 1)
  stopifnot(ppi_add_nodes <= 50)


  qgenes_match <-
    oncoEnrichR::verify_query_genes(query,
                                    qsource = query_source,
                                    ignore_unknown = ignore_unknown,
                                    genedb = oncoEnrichR::genedb,
                                    uniprot_acc = oncoEnrichR::uniprot_xref)


  background_entrez <- NULL
  background_genes_match <- NULL
  if (!is.null(background_enrichment)) {
    background_genes_match <-
      oncoEnrichR::verify_query_genes(background_enrichment,
                                      qsource = background_enrichment_source,
                                      genedb = oncoEnrichR::genedb,
                                      uniprot_acc = oncoEnrichR::uniprot_xref)
    if (background_genes_match[["success"]] == -1) {
      return(-1)
    }
    background_entrez <- unique(background_genes_match[["found"]]$entrezgene)
  }

  if (qgenes_match[["success"]] == -1) {
    return(-1)
  }

  query_entrezgene <- unique(qgenes_match[["found"]]$entrezgene)
  query_symbol <- unique(qgenes_match[["found"]]$symbol)

  onc_rep <- init_report(project_title = project_title,
                         project_owner = project_owner,
                         project_description = project_description,
                         ppi_add_nodes = ppi_add_nodes,
                         show_ppi = show_ppi,
                         show_drugs_in_ppi = show_drugs_in_ppi,
                         show_gene_summary = show_gene_summary,
                         show_disease = show_disease,
                         show_drug = show_drug,
                         show_enrichment = show_enrichment,
                         show_tcga_aberration = show_tcga_aberration,
                         show_tcga_coexpression = show_tcga_coexpression,
                         show_crispr_lof = show_crispr_lof,
                         min_geneset_size = min_geneset_size,
                         max_geneset_size = max_geneset_size,
                         simplify_go = simplify_go,
                         #show_gtex_coexp = show_gtex_coexp,
                         show_complex = show_complex,
                         show_subcell_comp = show_subcell_comp,
                         background_enrichment_description = background_enrichment_description,
                         p_value_cutoff_enrichment = p_value_cutoff_enrichment,
                         p_value_adjustment_method = p_value_adjustment_method,
                         q_value_cutoff_enrichment = q_value_cutoff_enrichment)

  if (length(query_symbol) > 20) {
    onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] <-
      onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] + as.integer((length(query_symbol) - 20)/ 8.5)
    onc_rep[["config"]][["co_expression_gtex"]][["plot_height"]] <-
      onc_rep[["config"]][["co_expression_gtex"]][["plot_height"]] + as.integer((length(query_symbol) - 20)/ 8.5)
  }

  if (show_disease == T) {
    onc_rep[["data"]][["disease"]][["target"]] <-
      oncoEnrichR::target_disease_associations(query_symbol,
                                               genedb = oncoEnrichR::genedb)
  }

  if (show_drug == T) {
    onc_rep[["data"]][["drug"]][["target"]] <-
      oncoEnrichR::target_drug_associations(query_symbol,
                                               genedb = oncoEnrichR::genedb)
  }

  for (c in names(oncoEnrichR::msigdb[["COLLECTION"]])) {
    for (subcat in names(oncoEnrichR::msigdb[["COLLECTION"]][[c]])) {
      if (c == "C5") {
        enr <- oncoEnrichR::get_go_enrichment(query_entrezgene,
                                  background_entrez = background_entrez,
                                  min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
                                  max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
                                  q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
                                  p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
                                  p_value_adjustment_method = onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
                                  simplify = onc_rep[["config"]][["enrichment"]][["simplify_go"]],
                                  ontology = subcat,
                                  genedb = oncoEnrichR::genedb)
        if (!is.null(enr)) {
          onc_rep[["data"]][["enrichment"]][["go"]] <-
            dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["go"]], enr) %>%
            dplyr::distinct()
        }
      }else{
        db = paste0("MSIGdb/", c, "/",subcat)
        enr <- oncoEnrichR::get_universal_enrichment(query_entrezgene,
                                         genedb = oncoEnrichR::genedb,
                                         background_entrez = background_entrez,
                                         min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
                                         max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
                                         q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
                                         p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
                                         p_value_adjustment_method = onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
                                         TERM2GENE = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2GENE,
                                         TERM2NAME = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2NAME,
                                         TERM2SOURCE = oncoEnrichR::msigdb$TERM2SOURCE,
                                         dbsource = db)
        if (!is.null(enr)) {
          onc_rep[["data"]][["enrichment"]][["msigdb"]] <-
            dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["msigdb"]], enr) %>%
            dplyr::distinct()
        }
      }
    }
  }


  db <- "WikiPathways"
  onc_rep[["data"]][["enrichment"]][["wikipathwaydb"]] <-
    oncoEnrichR::get_universal_enrichment(query_entrezgene,
                              genedb = oncoEnrichR::genedb,
                              background_entrez = background_entrez,
                              min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
                              max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
                              q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
                              p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
                              p_value_adjustment_method = onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
                              TERM2GENE = oncoEnrichR::wikipathwaydb$TERM2GENE,
                              TERM2NAME = oncoEnrichR::wikipathwaydb$TERM2NAME,
                              TERM2SOURCE = oncoEnrichR::wikipathwaydb$TERM2SOURCE,
                              dbsource = db)

  db <- "KEGG"
  onc_rep[["data"]][["enrichment"]][["keggdb"]] <-
    oncoEnrichR::get_universal_enrichment(query_entrezgene,
                                          genedb = oncoEnrichR::genedb,
                                          background_entrez = background_entrez,
                                          min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
                                          max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
                                          q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
                                          p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
                                          p_value_adjustment_method = onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
                                          TERM2GENE = oncoEnrichR::keggdb$TERM2GENE,
                                          TERM2NAME = oncoEnrichR::keggdb$TERM2NAME,
                                          TERM2SOURCE = oncoEnrichR::keggdb$TERM2SOURCE,
                                          dbsource = db)


  # if (show_gtex_coexp == T) {
  #   gtex_results <-
  #     oncoEnrichR::gtex_co_expression(query_symbol, genedb = oncoEnrichR::genedb, gids = gtex_atlasassay_groups)
  #   onc_rep[["data"]][["co_expression_gtex"]][["plots"]] <- gtex_results[["plots"]]
  #   onc_rep[["data"]][["co_expression_gtex"]][["df"]] <- gtex_results[["df"]]
  # }

  if (show_ppi == T) {
    onc_rep[["data"]][["ppi"]] <-
      oncoEnrichR::get_ppi_network(query_entrezgene,
                                   ppi_source = "STRING",
                                   genedb = oncoEnrichR::genedb,
                                   cancerdrugdb = oncoEnrichR::cancerdrugdb,
                                   settings = onc_rep[["config"]][["ppi"]][["stringdb"]])
  }

  if (show_complex == T) {
    onc_rep[["data"]][["protein_complex"]][["complex"]] <-
      oncoEnrichR::annotate_protein_complex(query_symbol,
                                            genedb = oncoEnrichR::genedb,
                                            corum_db = oncoEnrichR::corumdb,
                                            uniprot_acc = oncoEnrichR::uniprot_xref)
  }

  if (show_subcell_comp == T) {
     subcellcomp_annotations <-
      oncoEnrichR::annotate_subcellular_compartments(query_symbol,
                                            genedb = oncoEnrichR::genedb,
                                            comppidb = oncoEnrichR::comppidb)

     onc_rep[["data"]][["subcellcomp"]][["all"]] <- subcellcomp_annotations[["all"]]
     onc_rep[["data"]][["subcellcomp"]][["grouped"]] <- subcellcomp_annotations[["grouped"]]
     onc_rep[["data"]][["subcellcomp"]][["anatogram"]] <- subcellcomp_annotations[["anatogram"]]

  }

  if (show_crispr_lof == T) {
    onc_rep[["data"]][["loss_of_fitness"]] <-
      oncoEnrichR::get_crispr_lof_scores(query_symbol, projectscoredb = oncoEnrichR::projectscoredb)

    if (onc_rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] >= 20) {
      onc_rep[["config"]][["loss_of_fitness"]][["plot_height"]]  <-
        onc_rep[["config"]][["loss_of_fitness"]][["plot_height"]] + as.integer((onc_rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] - 20)/8.5)
    }
  }


  if (show_tcga_aberration == T) {
    for (v in c("cna_homdel","cna_ampl")) {
      onc_rep[["data"]][["tcga"]][["aberration"]][["plot"]][[v]] <-
        oncoEnrichR::tcga_aberration_plot(query_entrezgene, qsource = "entrezgene",
                                          genedb = oncoEnrichR::genedb, vtype = v)
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <-
        oncoEnrichR::tcga_aberration_table(query_entrezgene, qsource = "entrezgene",
                                           genedb = oncoEnrichR::genedb, vtype = v)
    }

    for(psite in names(onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]])){
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]][[psite]][['top_mutated_genes']] <-
        oncoEnrichR::tcga_oncoplot_genes(query_symbol, qsource = "symbol", genedb = oncoEnrichR::genedb,
                                         site = psite)
    }
  }

  if (show_tcga_coexpression == T) {
    onc_rep[["data"]][["tcga"]][["co_expression"]] <-
      oncoEnrichR::tcga_co_expression(query_symbol, genedb = oncoEnrichR::genedb)
  }

  return(onc_rep)

}

#' Function that writes the contents in the oncoEnrichR report object to an interactive HTML report
#'
#' @param report object with oncoEnrichR report data (returned by oncoEnrichR::onco_enrich)
#' @param project_directory working directory
#' @param report_name prefix filename for report output
#' @param format file format of output (html/json)
#' @export

write <- function(report, project_directory, report_name, format = "html") {

  outfname <- list()
  outfname[["html"]] <- paste(report_name, "html", sep=".")
  disclaimer <- system.file("templates", "disclaimer.md", package = "oncoEnrichR")
  report_theme <- "default"

  ## TODO: check that report parameter is a valid oncoEnrichR result object

  assign("onc_enrich_report", report, envir = .GlobalEnv)

  if (format == "html") {
    rlogging::message("------")
    rlogging::message("Writing HTML file with report contents")
    markdown_input <- system.file("templates", "onco_enrich_report.Rmd", package = "oncoEnrichR")
    rmarkdown::render(markdown_input,
                      output_format = rmarkdown::html_document(theme = report_theme, toc = T, toc_depth = 3,
                                                               toc_float = T, number_sections = F,
                                                               includes = rmarkdown::includes(after_body = disclaimer)),
                      output_file = outfname[["html"]],
                      output_dir = project_directory,
                      clean = T,
                      intermediates_dir = project_directory,
                      quiet = T)
    rlogging::message(paste0("Output file: ", file.path(project_directory,outfname[["html"]])))
    rlogging::message("------")
  }else{
    rlogging::message("JSON output not yet implemented")
  }
}
