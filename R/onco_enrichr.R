
init_report <- function(project_title = "Project title",
                        project_owner = "Project owner",
                        project_description = "Project description",
                        ppi_min_string_score = 900,
                        ppi_add_nodes = 50,
                        bgset_description =
                          "All protein-coding genes",
                        p_value_cutoff_enrichment = 0.05,
                        p_value_adjustment_method = "BH",
                        q_value_cutoff_enrichment = 0.2,
                        min_geneset_size = 10,
                        max_geneset_size = 500,
                        min_subcellcomp_confidence = 1,
                        simplify_go = F,
                        show_ppi = T,
                        show_drugs_in_ppi = T,
                        show_disease = T,
                        show_drug = T,
                        show_enrichment = T,
                        show_tcga_aberration = T,
                        show_tcga_coexpression = T,
                        show_subcell_comp = T,
                        show_crispr_lof = T,
                        show_cell_tissue = T,
                        show_unknown_function = T,
                        show_prognostic_cancer_assoc = T,
                        show_complex = T) {

  ## report object
  rep <- list()

  ## two main elements of report object
  # 1. data - contains all annotations and enrichment results
  # 2. config - contains all settings and underlying data sources
  for (e in c("data","config")) {
    rep[[e]] <- list()
  }

  ## config/resources - release notes - software and database versions
  rep[["config"]][["resources"]] <-  oncoEnrichR::release_notes

  ## config/show - logicals indicating which sections/analyses of the report to include
  rep[["config"]][["show"]] <- list()
  rep[["config"]][["show"]][["ppi"]] <- show_ppi
  rep[["config"]][["show"]][["disease"]] <- show_disease
  rep[["config"]][["show"]][["drug"]] <- show_drug
  rep[["config"]][["show"]][["enrichment"]] <- show_enrichment
  rep[["config"]][["show"]][["protein_complex"]] <- show_complex
  rep[["config"]][["show"]][["tcga_aberration"]] <- show_tcga_aberration
  rep[["config"]][["show"]][["tcga_coexpression"]] <- show_tcga_coexpression
  rep[["config"]][["show"]][["subcellcomp"]] <- show_subcell_comp
  rep[["config"]][["show"]][["loss_of_fitness"]] <- show_crispr_lof
  rep[["config"]][["show"]][["cell_tissue"]] <- show_cell_tissue
  rep[["config"]][["show"]][["unknown_function"]] <- show_unknown_function
  rep[["config"]][["show"]][["cancer_prognosis"]] <-
    show_prognostic_cancer_assoc

  ## config - report metadata (project owner, background, project_title)
  rep[["config"]][["project_title"]] <- project_title
  rep[["config"]][["project_description"]] <- project_description
  rep[["config"]][["project_owner"]] <- project_owner

  ## config/disease - color codes and
  ## thresholds for quantitative target-disease associations
  rep[["config"]][["disease"]] <- list()
  rep[["config"]][["disease"]][["breaks"]] <- c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  rep[["config"]][["disease"]][["colors"]] <- c("#b8b8ba","#EFF3FF","#C6DBEF",
                                                "#9ECAE1","#6BAED6","#4292C6",
                                                "#2171B5","#084594")

  ## config/loss_of_fitness - plot height for hits in CRISPR screens
  rep[["config"]][["loss_of_fitness"]] <- list()
  rep[["config"]][["loss_of_fitness"]][["plot_height"]] <- 10

  ## config/ppi - protein-protein interaction settings
  rep[["config"]][["ppi"]] <- list()
  rep[["config"]][["ppi"]][["stringdb"]] <- list()
  rep[["config"]][["ppi"]][["stringdb"]][["minimum_score"]] <- ppi_min_string_score
  rep[["config"]][["ppi"]][["stringdb"]][["visnetwork_shape"]] <- "dot"
  rep[["config"]][["ppi"]][["stringdb"]][["visnetwork_shadow"]] <- T
  rep[["config"]][["ppi"]][["stringdb"]][["show_drugs"]] <- show_drugs_in_ppi
  rep[["config"]][["ppi"]][["stringdb"]][["add_nodes"]] <- ppi_add_nodes
  rep[["config"]][["ppi"]][["stringdb"]][["query_type"]] <- "network"

  ## config/enrichment - general functional enrichment settings
  rep[["config"]][["enrichment"]] <- list()
  rep[["config"]][["enrichment"]][["p_value_cutoff"]] <-
    p_value_cutoff_enrichment
  rep[["config"]][["enrichment"]][["q_value_cutoff"]] <-
    q_value_cutoff_enrichment
  rep[["config"]][["enrichment"]][["p_adjust_method"]] <-
    p_value_adjustment_method
  rep[["config"]][["enrichment"]][["min_gs_size"]] <-
    min_geneset_size
  rep[["config"]][["enrichment"]][["max_gs_size"]] <-
    max_geneset_size
  rep[["config"]][["enrichment"]][["simplify_go"]] <-
    simplify_go
  rep[["config"]][["enrichment"]][["bgset_description"]] <-
    bgset_description

  ## specifiy plot height (tile/heatmap) - genes along y-axis
  ## tumor types / tissues / celltypes along x-axis
  rep[["config"]][["tcga_aberration"]] <- list()
  rep[["config"]][["tcga_aberration"]][["plot_height"]] <- 14

  ## config/prognosis - settings for prognostic associations
  rep[["config"]][["prognosis"]] <- list()
  rep[["config"]][["prognosis"]][["breaks"]] <-
    c(3, 5.2, 7.4, 9.5, 11.7, 13.9)
  rep[["config"]][["prognosis"]][["colors_unfavorable"]] <-
    c("#fee5d9","#fcbba1","#fc9272",
      "#fb6a4a","#ef3b2c","#cb181d","#99000d")
  rep[["config"]][["prognosis"]][["colors_favorable"]] <-
    c("#edf8e9","#c7e9c0","#a1d99b",
      "#74c476","#41ab5d","#238b45","#005a32")

  rep[["config"]][["cell_tissue"]] <- list()
  rep[["config"]][["cell_tissue"]][["tissue_enrichment_levels"]] <-
    c('Tissue enriched',
      'Group enriched',
      'Tissue enhanced',
      'Mixed',
      'Low tissue specificity',
      'Not detected')
  rep[["config"]][["cell_tissue"]][["ctype_enrichment_levels"]] <-
    c('Cell type enriched',
      'Group enriched',
      'Cell type enhanced',
      'Mixed',
      'Low cell type specificity',
      'Not detected')
  rep[["config"]][["cell_tissue"]][["enrichment_colors"]] <-
    c("#084594","#2171B5","#4292C6",
      "#6BAED6","#9ECAE1","#b8b8ba")

  rep[["config"]][["unknown_function"]] <- list()
  rep[["config"]][["unknown_function"]][["rank"]] <-
    c(1, 2, 3, 4, 5, 6)
  rep[["config"]][["unknown_function"]][["colors"]] <-
    c("#99000d","#cb181d","#ef3b2c","#fb6a4a","#fc9272","#fcbba1")

  rep[['config']][['subcellcomp']] <- list()
  rep[['config']][['subcellcomp']][['minimum_confidence']] <-
    min_subcellcomp_confidence

  ## initialize all data elements
  for (analysis in c("tcga",
                     "disease",
                     "ppi",
                     "tcga",
                     "enrichment",
                     "drug",
                     "protein_complex",
                     "subcellcomp",
                     "loss_of_fitness",
                     "cell_tissue",
                     "cancer_prognosis",
                     "unknown_function")) {
    rep[["data"]][[analysis]] <- list()
  }

  ## prognosis/survival - gene expression (HPA)
  rep[["data"]][["cancer_prognosis"]][['assocs']] <- data.frame()

  ## tissue and cell type enrichment
  rep[['data']][['cell_tissue']] <- list()
  rep[['data']][['cell_tissue']][['tissue_overview']] <- list()
  rep[['data']][['cell_tissue']][['tissue_enrichment']] <- list()
  rep[['data']][['cell_tissue']][['tissue_enrichment']][['per_gene']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['tissue_enrichment']][['per_type']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['tissue_overview']][['category_df']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['tissue_overview']][['exp_dist_df']] <-
    data.frame()

  rep[['data']][['cell_tissue']][['scRNA_overview']] <- list()
  rep[['data']][['cell_tissue']][['scRNA_enrichment']] <- list()
  rep[['data']][['cell_tissue']][['scRNA_enrichment']][['per_gene']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['scRNA_enrichment']][['per_type']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['scRNA_overview']][['category_df']] <-
    data.frame()
  rep[['data']][['cell_tissue']][['scRNA_overview']][['exp_dist_df']] <-
    data.frame()

  ## drug targets
  rep[["data"]][["drug"]][["target"]] <- data.frame()

  ## disease associations
  rep[["data"]][["disease"]][["target"]] <- data.frame()
  rep[["data"]][["disease"]][["target_assoc"]] <- data.frame()
  rep[["data"]][["disease"]][["assoc_pr_gene"]] <- list()

  ## protein-protein interactions
  rep[["data"]][["ppi"]][["complete_network"]] <- NULL
  rep[["data"]][["ppi"]][["hubscores"]] <- data.frame()
  rep[["data"]][["ppi"]][["community_network"]] <- NULL

  ## functional enrichment
  for (c in c("go","msigdb","wikipathway","kegg","netpath")) {
    rep[["data"]][["enrichment"]][[c]] <- data.frame()
  }

  ## unknown function
  rep[["data"]][["unknown_function"]][["hits_df"]] <- data.frame()

  ## CRISPR/Cas9 LOF hits
  rep[["data"]][["loss_of_fitness"]][["hits_df"]] <- data.frame()
  rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] <- 0

  ## TCGA co-expression
  rep[["data"]][["tcga"]][["co_expression"]] <- data.frame()

  ## Protein complexes
  rep[["data"]][["protein_complex"]][["complex"]] <- data.frame()

  ## Subcellular localizations
  rep[["data"]][["subcellcomp"]][["all"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["grouped"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["anatogram"]] <- data.frame()

  ## TCGA aberrations
  rep[["data"]][["tcga"]][["aberration"]] <- list()
  rep[["data"]][["tcga"]][["aberration"]][["table"]] <- list()
  rep[["data"]][["tcga"]][["aberration"]][["matrix"]] <- list()
  rep[["data"]][["tcga"]][["aberration"]][["plot"]] <- list()

  for (v in c("cna_homdel","cna_ampl","snv_indel")) {
    if(v != "snv_indel"){
      rep[["data"]][["tcga"]][["aberration"]][["matrix"]][[v]] <- NULL
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- data.frame()
    }
    else{
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- list()
      i <- 1
      while(i <= nrow(oncoEnrichR::maf_codes)){
        site <- oncoEnrichR::maf_codes[i,]$primary_site
        code <- oncoEnrichR::maf_codes[i,]$code
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]] <-
          list()
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]][['code']] <-
          code
        rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]][[site]][['top_mutated_genes']] <-
          data.frame()
        i <- i + 1
      }
    }
  }

  return(rep)
}

#' Function that interrogates and analyzes a list of human protein-coding genes for cancer relevance
#'
#' @param query character vector with gene/query identifiers
#' @param query_id_type character indicating source of query (one of "uniprot_acc", "symbol", "entrezgene", or "ensembl_gene_id")
#' @param ignore_id_err logical indicating if analysis should continue when uknown query identifiers are encountered
#' @param project_title project title (title of report)
#' @param project_owner name of project owner
#' @param project_description project background information
#' @param bgset character vector with gene identifiers, used as reference/background for enrichment/over-representation analysis
#' @param bgset_id_type character indicating source of background ("uniprot_acc","symbol","entrezgene","ensembl_gene_id")
#' @param bgset_description character indicating type of background (e.g. "All lipid-binding proteins (n = 200)")
#' @param p_value_cutoff_enrichment cutoff p-value for enrichment/over-representation analysis
#' @param p_value_adjustment_method one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param q_value_cutoff_enrichment cutoff q-value for enrichment analysis
#' @param min_geneset_size minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param max_geneset_size maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param min_subcellcomp_confidence minimum confidence level of subcellular compartment annotations (range from 1 to 6, 6 is strongest)
#' @param simplify_go remove highly similar GO terms in results from GO enrichment/over-representation analysis
#' @param ppi_add_nodes number of nodes to add to target set when computing the protein-protein interaction network (STRING)
#' @param ppi_score_threshold minimum score (0-1000) for retrieval of protein-protein interactions (STRING)
#' @param show_ppi logical indicating if report should contain protein-protein interaction data (STRING)
#' @param show_drugs_in_ppi logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform)
#' @param show_enrichment logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME etc.)
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of queryset with oncogenes/tumor suppressor genes
#' @param show_tissue_cell logical indicating if report should contain tissue-specificity and single cell-type specificity assessments
#' of target genes, using data from the Human Protein Atlas
#' @param show_unknown_function logical indicating if report should highlight target genes with unknown or poorly defined functions
#' @param show_prognostic_cancer_assoc  logical indicating if mRNA-based (single-gene) prognostic associations to cancer types should be listed
#' @param show_subcell_comp logical indicating if report should list subcellular compartment annotations (ComPPI)
#' @param show_crispr_lof logical indicating if report should list results from CRISPR/Cas9 loss-of-fitness screens (Project Score)
#' @param show_complex logical indicating if report should list proteins in known protein complexes (CORUM)
#' @export
#'
onco_enrich <- function(query,
                   query_id_type = "symbol",
                   ignore_id_err = FALSE,
                   project_title = "Project title",
                   project_owner = "Project owner",
                   project_description = "Project description",
                   bgset = NULL,
                   bgset_id_type = "symbol",
                   bgset_description = "All protein-coding genes",
                   p_value_cutoff_enrichment = 0.05,
                   p_value_adjustment_method = "BH",
                   q_value_cutoff_enrichment = 0.2,
                   min_geneset_size = 10,
                   max_geneset_size = 500,
                   min_subcellcomp_confidence = 1,
                   simplify_go = F,
                   ppi_add_nodes = 50,
                   ppi_score_threshold = 900,
                   show_ppi = T,
                   show_drugs_in_ppi = F,
                   show_disease = T,
                   show_drug = T,
                   show_enrichment = T,
                   show_tcga_aberration = T,
                   show_tcga_coexpression = T,
                   show_cell_tissue = T,
                   show_unknown_function = T,
                   show_prognostic_cancer_assoc = T,
                   show_subcell_comp = T,
                   show_crispr_lof = T,
                   show_complex = T) {
  stopifnot(is.character(query))
  stopifnot(p_value_adjustment_method %in% c("holm", "hochberg",
                                             "hommel", "bonferroni",
                                             "BH", "BY",
                                             "fdr", "none"))
  if (length(query) > 800 | length(query) < 20) {
    rlogging::message(paste0("ERROR: oncoEnrichR needs minimum 20 query ",
    "identifiers, and accepts a maximum of 800. Query contained n = ",
    length(query), " identifiers"))
    return(NULL)
  }
  stopifnot(query_id_type == "symbol" | query_id_type == "entrezgene" |
              query_id_type == "uniprot_acc" | query_id_type == "ensembl_gene_id")
  stopifnot(ppi_score_threshold > 0 & ppi_score_threshold <= 1000)
  stopifnot(p_value_cutoff_enrichment > 0 & p_value_cutoff_enrichment < 1)
  stopifnot(q_value_cutoff_enrichment > 0 & q_value_cutoff_enrichment < 1)
  stopifnot(min_subcellcomp_confidence >= 1 & min_subcellcomp_confidence <= 6)
  stopifnot(ppi_add_nodes <= 50)

  qgenes_match <-
    oncoEnrichR:::verify_query_genes(query,
                                     q_id_type = query_id_type,
                                    ignore_id_err = ignore_id_err,
                                    genedb = oncoEnrichR::genedb,
                                    uniprot_acc = oncoEnrichR::uniprot_xref)

  if (qgenes_match[["success"]] == -1) {
    return(NULL)
  }


  background_entrez <- NULL
  background_genes_match <- NULL
  if (!is.null(bgset)) {
    background_genes_match <-
      oncoEnrichR:::verify_query_genes(bgset,
                                       q_id_type = bgset_id_type,
                                      genedb = oncoEnrichR::genedb,
                                      qtype = "background",
                                      uniprot_acc = oncoEnrichR::uniprot_xref)
    if (background_genes_match[["success"]] == -1) {
      return(NULL)
    }
    background_entrez <- unique(background_genes_match[["found"]]$entrezgene)
  }


  query_entrezgene <- unique(qgenes_match[["found"]]$entrezgene)
  query_symbol <- unique(qgenes_match[["found"]]$symbol)

  onc_rep <- init_report(project_title = project_title,
                         project_owner = project_owner,
                         project_description = project_description,
                         ppi_add_nodes = ppi_add_nodes,
                         show_ppi = show_ppi,
                         show_drugs_in_ppi = show_drugs_in_ppi,
                         show_disease = show_disease,
                         show_drug = show_drug,
                         show_enrichment = show_enrichment,
                         show_tcga_aberration = show_tcga_aberration,
                         show_tcga_coexpression = show_tcga_coexpression,
                         show_cell_tissue = show_cell_tissue,
                         show_unknown_function = show_unknown_function,
                         show_prognostic_cancer_assoc =
                           show_prognostic_cancer_assoc,
                         show_crispr_lof = show_crispr_lof,
                         min_geneset_size = min_geneset_size,
                         max_geneset_size = max_geneset_size,
                         min_subcellcomp_confidence = min_subcellcomp_confidence,
                         simplify_go = simplify_go,
                         show_complex = show_complex,
                         show_subcell_comp = show_subcell_comp,
                         bgset_description =
                           bgset_description,
                         p_value_cutoff_enrichment = p_value_cutoff_enrichment,
                         p_value_adjustment_method = p_value_adjustment_method,
                         q_value_cutoff_enrichment = q_value_cutoff_enrichment)

  if (length(query_symbol) > 20) {
    onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] <-
      onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] +
      as.integer((length(query_symbol) - 20)/ 7.5)
  }

  if (show_disease == T) {
    onc_rep[["data"]][["disease"]][["target"]] <-
      oncoEnrichR:::target_disease_associations(
        query_symbol,
        genedb = oncoEnrichR::genedb)
  }

  if (show_drug == T) {
    onc_rep[["data"]][["drug"]][["target"]] <-
      oncoEnrichR:::target_drug_associations(
        query_symbol,
        genedb = oncoEnrichR::genedb)
  }

  for (c in names(oncoEnrichR::msigdb[["COLLECTION"]])) {
    for (subcat in names(oncoEnrichR::msigdb[["COLLECTION"]][[c]])) {
      if (c == "C5" & subcat != "HPO") {
        enr <- oncoEnrichR:::get_go_enrichment(
          query_entrezgene,
          background_entrez = background_entrez,
          min_geneset_size =
            onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
          max_geneset_size =
            onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
          q_value_cutoff =
            onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
          p_value_cutoff =
            onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
          p_value_adjustment_method =
            onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
          simplify = onc_rep[["config"]][["enrichment"]][["simplify_go"]],
          ontology = subcat,
          genedb = oncoEnrichR::genedb)
        if (!is.null(enr)) {
          onc_rep[["data"]][["enrichment"]][["go"]] <-
            dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["go"]], enr) %>%
            dplyr::distinct()
        }
      }else{
        if(c != "C5"){
          db = paste0("MSIGdb/", c, "/",subcat)
          enr <- oncoEnrichR:::get_universal_enrichment(
            query_entrezgene,
            genedb = oncoEnrichR::genedb,
            background_entrez = background_entrez,
            min_geneset_size =
              onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
            max_geneset_size =
              onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
            q_value_cutoff =
              onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
            p_value_cutoff =
              onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
            p_value_adjustment_method =
              onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
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
  }


  db <- "WikiPathways"
  onc_rep[["data"]][["enrichment"]][["wikipathway"]] <-
    oncoEnrichR:::get_universal_enrichment(
      query_entrezgene,
      genedb = oncoEnrichR::genedb,
      background_entrez = background_entrez,
      min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
      max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
      q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
      p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
      p_value_adjustment_method =
        onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
      TERM2GENE = oncoEnrichR::wikipathwaydb$TERM2GENE,
      TERM2NAME = oncoEnrichR::wikipathwaydb$TERM2NAME,
      TERM2SOURCE = oncoEnrichR::wikipathwaydb$TERM2SOURCE,
      dbsource = db)

  db <- "NetPath"
  onc_rep[["data"]][["enrichment"]][["netpath"]] <-
    oncoEnrichR:::get_universal_enrichment(
      query_entrezgene,
      genedb = oncoEnrichR::genedb,
      background_entrez = background_entrez,
      min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
      max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
      q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
      p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
      p_value_adjustment_method =
        onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
      TERM2GENE = oncoEnrichR::netpathdb$TERM2GENE,
      TERM2NAME = oncoEnrichR::netpathdb$TERM2NAME,
      TERM2SOURCE = oncoEnrichR::netpathdb$TERM2SOURCE,
      dbsource = db)

  db <- "KEGG"
  onc_rep[["data"]][["enrichment"]][["kegg"]] <-
    oncoEnrichR:::get_universal_enrichment(
      query_entrezgene,
      genedb = oncoEnrichR::genedb,
      background_entrez = background_entrez,
      min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
      max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
      q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
      p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
      p_value_adjustment_method =
        onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
      TERM2GENE = oncoEnrichR::keggdb$TERM2GENE,
      TERM2NAME = oncoEnrichR::keggdb$TERM2NAME,
      TERM2SOURCE = oncoEnrichR::keggdb$TERM2SOURCE,
      dbsource = db)

  if (show_ppi == T) {
    ## TEST TO CHECK OMNIPATHDB IS LIVE IS NOT WORKING (NOT SURE WHY)
    # service_is_down <- unique(is.na(pingr::ping("string-db.org")))
    # if(service_is_down){
    #   message("EXCEPTION: https://string-db.org is NOT responding - ",
    #           "skipping retrievel of protein-protein interactions")
    #
    #   onc_rep[["config"]][["show"]][["ppi"]] <- FALSE
    # }else{
    onc_rep[["data"]][["ppi"]] <-
      oncoEnrichR:::get_ppi_network(
        query_entrezgene,
        ppi_source = "STRING",
        genedb = oncoEnrichR::genedb,
        cancerdrugdb = oncoEnrichR::cancerdrugdb,
        settings = onc_rep[["config"]][["ppi"]][["stringdb"]])
    #}
  }

  if (show_complex == T) {
    ## TEST TO CHECK OMNIPATHDB IS LIVE IS NOT WORKING (NOT SURE WHY)
    # service_is_down <- unique(is.na(pingr::ping("omnipathdb.org")))
    # if(service_is_down){
    #   rlogging::message("EXCEPTION: https://omnipathdb.org is NOT responding - ",
    #                     "skipping retrievel of protein complexes")
    #
    #   onc_rep[["config"]][["show"]][["protein_complex"]] <- FALSE
    # }else{
    onc_rep[["data"]][["protein_complex"]][["complex"]] <-
      oncoEnrichR:::annotate_protein_complex(
        query_symbol,
        genedb = oncoEnrichR::genedb,
        corum_db = oncoEnrichR::corumdb,
        uniprot_acc = oncoEnrichR::uniprot_xref)
    #}
  }

  if (show_unknown_function == T) {
    onc_rep[["data"]][["unknown_function"]][["hits_df"]] <-
      oncoEnrichR:::get_genes_unknown_function(
        query_symbol,
        genedb = oncoEnrichR::genedb,
        poorly_defined_genes = oncoEnrichR::poorly_defined_genes
      )
  }

  if (show_subcell_comp == T) {
     subcellcomp_annotations <-
      oncoEnrichR:::annotate_subcellular_compartments(
        query_symbol,
        minimum_confidence = min_subcellcomp_confidence,
        genedb = oncoEnrichR::genedb,
        comppidb = oncoEnrichR::subcellcomp$comppidb)

     onc_rep[["data"]][["subcellcomp"]][["all"]] <-
       subcellcomp_annotations[["all"]]
     onc_rep[["data"]][["subcellcomp"]][["grouped"]] <-
       subcellcomp_annotations[["grouped"]]
     onc_rep[["data"]][["subcellcomp"]][["anatogram"]] <-
       subcellcomp_annotations[["anatogram"]]

  }

  if (show_crispr_lof == T) {
    onc_rep[["data"]][["loss_of_fitness"]] <-
      oncoEnrichR:::get_crispr_lof_scores(
        query_symbol,
        projectscoredb = oncoEnrichR::projectscoredb)

    if (onc_rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] >= 20) {
      onc_rep[["config"]][["loss_of_fitness"]][["plot_height"]]  <-
        onc_rep[["config"]][["loss_of_fitness"]][["plot_height"]] +
        as.integer((onc_rep[["data"]][["loss_of_fitness"]][["n_genes_with_hits"]] - 20)/8.5)
    }
  }


  if (show_tcga_aberration == T) {
    for (v in c("cna_homdel","cna_ampl")) {
      onc_rep[["data"]][["tcga"]][["aberration"]][["matrix"]][[v]] <-
        oncoEnrichR:::tcga_aberration_matrix(
          query_entrezgene,
          qsource = "entrezgene",
          genedb = oncoEnrichR::genedb,
          vtype = v)
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <-
        oncoEnrichR:::tcga_aberration_table(
          query_entrezgene,
          qsource = "entrezgene",
          genedb = oncoEnrichR::genedb,
          vtype = v)
    }

    for(psite in names(onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]])){
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]][[psite]][['top_mutated_genes']] <-
        oncoEnrichR:::tcga_oncoplot_genes(
          query_symbol,
          qsource = "symbol",
          genedb = oncoEnrichR::genedb,
          site = psite)
    }
  }

  if (show_tcga_coexpression == T) {
    onc_rep[["data"]][["tcga"]][["co_expression"]] <-
      oncoEnrichR:::tcga_co_expression(
        query_symbol,
        genedb = oncoEnrichR::genedb)
  }

  if(show_prognostic_cancer_assoc == T){
    onc_rep[["data"]][["cancer_prognosis"]][['assocs']] <-
      oncoEnrichR:::hpa_prognostic_genes(
        query_symbol,
        genedb = oncoEnrichR::genedb)

  }

  if(show_cell_tissue == T){
    onc_rep[["data"]][["cell_tissue"]][['tissue_overview']] <-
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        query_symbol,
        genedb = oncoEnrichR::genedb)

    onc_rep[["data"]][["cell_tissue"]][['tissue_enrichment']] <-
      oncoEnrichR:::gene_tissue_cell_enrichment(
        query_entrezgene,
        resolution = "tissue",
        background_entrez = background_entrez,
        genedb = oncoEnrichR::genedb)

    onc_rep[["data"]][["cell_tissue"]][['scRNA_overview']] <-
      oncoEnrichR:::gene_tissue_cell_spec_cat(
        query_symbol,
        resolution = "single_cell",
        genedb = oncoEnrichR::genedb)

    onc_rep[["data"]][["cell_tissue"]][['scRNA_enrichment']] <-
      oncoEnrichR:::gene_tissue_cell_enrichment(
        query_entrezgene,
        background_entrez = background_entrez,
        resolution = "single_cell",
        genedb = oncoEnrichR::genedb)

  }

  return(onc_rep)

}

#' Function that writes the contents of the oncoEnrichR report object to either
#' A) interactive HTML report, or B) Excel workbook
#'
#' @param report object with oncoEnrichR report data (returned by oncoEnrichR::onco_enrich)
#' @param file filename for report output
#' @param ignore_file_extension logical to accept any type of filaname extensions (for Galaxy integration)
#' @param overwrite logical indicating if file contents may be overwritten
#' @param format file format of output (html/excel)
#' @export

write <- function(report,
                  #project_directory,
                  file = "testReport.html",
                  ignore_file_extension = T,
                  overwrite = T,
                  format = "html") {


  val <- assertthat::validate_that(
    format %in% c("html","excel")
  )
  if(!is.logical(val)){
    message(val)
  }

  val <- assertthat::validate_that(
      is.character(file)
  )
  if(!is.logical(val)){
    message(val)
  }

  val <- assertthat::validate_that(
    is.logical(overwrite)
  )
  if(!is.logical(val)){
    message(val)
  }

  output_directory <- dirname(file)
  file_basename <- basename(file)
  if(output_directory != "."){
    val <- assertthat::validate_that(
      dir.exists(output_directory)
    )
    if(!is.logical(val)){
      message(val)
    }
  }

  if(overwrite == F){
    val <- assertthat::validate_that(
      file.exists(file)
    )
    if(!is.logical(val)){
      message(val)
    }
  }

  if(ignore_file_extension == F){
    if(format == "html" & tools::file_ext(file) != "html"){
      message("oncoEnrichR HTML output: File name must end with .html, not ",file)
    }

    if(format == "excel" & tools::file_ext(file) != "xlsx"){
      message("oncoEnrichR Excel output: File name must end with .xlsx, not ",file)
    }
  }

  ## TODO: check that report parameter is a valid oncoEnrichR result object

  if(!is.null(report)){
    assign("onc_enrich_report",
           report, envir = .GlobalEnv)
  }else{
    message("ERROR: report object is NULL - cannot write report contents")
  }

  if (format == "html") {

    disclaimer <- system.file("templates",
                              "disclaimer.md",
                              package = "oncoEnrichR")
    report_theme <- "default"

    rlogging::message("------")
    rlogging::message("Writing HTML file with report contents")
    markdown_input <- system.file("templates", "onco_enrich_report.Rmd",
                                  package = "oncoEnrichR")
    rmarkdown::render(
      markdown_input,
      output_format = rmarkdown::html_document(
        theme = report_theme, toc = T, toc_depth = 3,
        toc_float = T, number_sections = F,
        includes = rmarkdown::includes(after_body = disclaimer)),
      output_file = file_basename,
      output_dir = output_directory,
      clean = T,
      intermediates_dir = output_directory,
      quiet = T)
    rlogging::message(paste0("Output file: ",
                             file))
    rlogging::message("------")
  }
  if (format == "excel") {

    wb <- openxlsx::createWorkbook()
    rlogging::message("------")
    rlogging::message("Writing Excel workbook with report contents")

    table_style_index <- 15
    for(elem in c("disease",
                  "tcga_aberration",
                  "tcga_coexpression",
                  "cancer_prognosis",
                  "enrichment",
                  "drug",
                  "protein_complex",
                  "subcellcomp",
                  "cell_tissue",
                  "unknown_function")){

      show_elem <- elem

      if(report[['config']][['show']][[show_elem]] == FALSE){
        next
      }

      wb <- oncoEnrichR:::add_excel_sheet(
        report = report,
        workbook = wb,
        analysis_output = elem,
        tableStyle =
          paste0("TableStyleMedium",
                table_style_index)
      )

      ## Excel Table style: TableStyleMedium - 15 to 21                                                                                               table_style_index))
      if(table_style_index == 21){
        table_style_index <- 14
      }
      table_style_index <- table_style_index + 1

    }

    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    rlogging::message(paste0("Output file: ",file))
    rlogging::message("------")
  }
  else{
    if(format == "json"){
      rlogging::message("JSON output not yet implemented")
    }
  }
}
