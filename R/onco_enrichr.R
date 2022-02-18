#' Load oncoEnrichR annotation database
#'
#' @param github If TRUE (default), load database from GitHub repository
#' @param cache_dir path to local cache for more efficient subsequent loads
#'
#' @export
#'
load_db <- function(github = T,
                    cache_dir = NULL){

  logger <- log4r::logger(threshold = "INFO",
                          appenders = log4r::console_appender(log4r_layout))

  log4r_info(logger, paste0("Loading oncoEnrichR annotation databases"))

  if(!is.null(cache_dir)){
    val <- assertthat::validate_that(
      dir.exists(cache_dir)
    )
    if(!is.logical(val)){
      log4r_info(paste0("ERROR: ",val))
      return()
    }

  }

  oe_version <- paste0("v", utils::packageVersion("oncoEnrichR"))
  db_url <- paste0("https://github.com/sigven/oncoEnrichR/raw/",
                     oe_version, "/db/")
  if(!RCurl::url.exists(db_url)){
    if(is.null(cache_dir)){
      log4r_warn(logger, paste0("Database for version '", oe_version,
                              "' is not tagged on GitHub - using master branch"))
    }

    db_url <- paste0("https://github.com/sigven/oncoEnrichR/raw/",
                     "master/db/")
  }

  oedb <- list()

  cache_in_use_log_printed <- 0

  for(db in c("cancerdrugdb",
              "tcgadb",
              "projectsurvivaldb",
              "pathwaydb",
              "genedb",
              "subcelldb",
              "tissuecelldb",
              "projectscoredb",
              "tftargetdb",
              "cancerdrugdb",
              "otdb",
              "hpa",
              "ligandreceptordb",
              "release_notes")){


    cache_found <- 0

    if(!is.null(cache_dir)){
      db_cache_path <- file.path(
        cache_dir,
        oe_version,
        paste0(db,".rds"))

      if(file.exists(db_cache_path)){
        oedb[[db]] <- readRDS(file = db_cache_path)


        if(cache_in_use_log_printed == 0){
          log4r_info(logger, paste0("Reading from local cache - ",
                                    file.path(cache_dir, oe_version)))
          cache_in_use_log_printed <- 1
        }
        cache_found <- 1

      }
    }

    if(cache_found == 0){

      subdb_url <- paste0(db_url,
                    db,".rds")
      if(RCurl::url.exists(subdb_url)){
        oedb[[db]] <- readRDS(url(subdb_url,"rb"))
      }else{
        log4r_err(logger, paste0("Could not retrieve data from ", subdb_url))
      }
    }


    checksum_db <- R.cache::getChecksum(oedb[[db]])

    ##subcelldb's checksum does not recover correctly with all
    ##list entries involved, choose comppi only
    if(db == 'subcelldb'){
      if('comppidb' %in% names(oedb[[db]])){
        checksum_db <- R.cache::getChecksum(oedb[[db]][['comppidb']])
      }
    }


    if(checksum_db ==
       oncoEnrichR::db_props[oncoEnrichR::db_props$name == db,"checksum"]){

      log4r_info(logger, paste0("'",
                                db, "' - ",
                                oe_version, " - ",
                                checksum_db,
                                " - verifies correctly"))
    }else{
      log4r_err(logger, paste0("'",
                                db, "' - ",
                                oe_version, " - ",
                                checksum_db,
                                " - does not verify correctly"))
    }


    if (cache_found == 0 & !is.null(cache_dir)){

      if(!dir.exists(file.path(cache_dir, oe_version))){
        system(paste0('mkdir ', file.path(cache_dir, oe_version)),
               intern = F)
      }

      saveRDS(oedb[[db]], file =
                file.path(cache_dir, oe_version,
                          paste0(db, ".rds")))

    }
  }


  return(oedb)

}


#' Function that initiates oncoEnrichR report structure
#'
#' @param oeDB oncoEnrichR annotation database object
#' @param project_title project title (title of report)
#' @param project_owner name of project owner
#' @param html_floating_toc logical - float the table of contents to the left of the main document content. The floating table of contents will always be visible even when the document is scrolled
#' @param html_report_theme Bootswatch theme for HTML report (any of "bootstrap","cerulean","cosmo","default","flatly","journal","lumen","paper","sandstone","simplex","spacelab","united","yeti")
#' @param query_id_type character indicating source of query (one of "uniprot_acc", "symbol",
#' "entrezgene", or "ensembl_gene","ensembl_mrna","refseq_mrna","ensembl_protein","refseq_protein")
#' @param ignore_id_err logical indicating if analysis should continue when uknown query identifiers are encountered
#' @param project_description project background information
#' @param ppi_min_string_score minimum score (0-1000) for retrieval of protein-protein interactions (STRING)
#' @param ppi_add_nodes number of nodes to add to target set when computing the protein-protein interaction network (STRING)
#' @param bgset_description character indicating type of background (e.g. "All lipid-binding proteins (n = 200)")
#' @param bgset_id_type character indicating source of query (one of "uniprot_acc", "symbol",
#' "entrezgene", or "ensembl_gene","ensembl_mrna","refseq_mrna","ensembl_protein","refseq_protein")
#' @param p_value_cutoff_enrichment cutoff p-value for enrichment/over-representation analysis
#' @param p_value_adjustment_method one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param q_value_cutoff_enrichment cutoff q-value for enrichment analysis
#' @param min_geneset_size minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param max_geneset_size maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param num_terms_enrichment_plot number of top enriched Gene Ontology terms (max) to show in enrichment barplot
#' @param min_subcellcomp_confidence minimum confidence level of subcellular compartment annotations (range from 1 to 6, 6 is strongest)
#' @param subcellcomp_show_cytosol logical indicating if subcellular heatmap should highlight cytosol as a subcellular protein location or not
#' @param min_confidence_reg_interaction minimum confidence level for regulatory interactions (TF-target) retrieved from DoRothEA ('A','B','C', or 'D')
#' @param simplify_go remove highly similar GO terms in results from GO enrichment/over-representation analysis
#' @param show_ppi logical indicating if report should contain protein-protein interaction data (STRING)
#' @param show_drugs_in_ppi logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform, association_score >= 0.05, support from at least two data types)
#' @param show_top_diseases_only logical indicating if report should contain top (n = 20) disease associations only pr. query gene (Open Targets Platform)
#' @param show_cancer_hallmarks logical indicating if report should contain annotations/evidence of cancer hallmarks per query gene (COSMIC/Open Targets Platform)
#' @param show_drug logical indicating if report should contain targeted cancer drug information
#' @param show_enrichment logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME, NetPath, WikiPathways)
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of query set with oncogenes/tumor suppressor genes
#' @param show_cell_tissue logical indicating if report should contain tissue-specificity and single cell-type specificity assessments (Human Protein Atlas)
#' of target genes, using data from the Human Protein Atlas
#' @param show_ligand_receptor logical indicating if report should contain ligand-receptor interactions (CellChatDB)
#' @param show_regulatory_interactions logical indicating if report should contain data on transcription factor (TF) - target interactions relevant for the query set (DoRothEA)
#' @param show_unknown_function logical indicating if report should highlight target genes with unknown or poorly defined functions (GO/Uniprot KB/NCBI)
#' @param show_prognostic_cancer_assoc  logical indicating if mRNA-based (single-gene) prognostic associations to cancer types should be listed (Human Protein Atlas/TCGA)
#' @param show_subcell_comp logical indicating if report should provide subcellular compartment annotations (ComPPI)
#' @param show_crispr_lof logical indicating if report should provide fitness scores and target priority scores from CRISPR/Cas9 loss-of-fitness screens (Project Score)
#' @param show_complex logical indicating if report should provide target memberships in known protein complexes (ComplexPortal/Compleat/PDB/CORUM)
#'
#' @keywords internal
#'

init_report <- function(oeDB,
                        project_title = "Project title",
                        project_owner = "Project owner",
                        html_floating_toc = T,
                        html_report_theme = "default",
                        query_id_type = "symbol",
                        ignore_id_err = TRUE,
                        project_description = "Project description",
                        ppi_min_string_score = 900,
                        ppi_add_nodes = 50,
                        bgset_description =
                          "All annotated protein-coding genes",
                        bgset_id_type = "symbol",
                        p_value_cutoff_enrichment = 0.05,
                        p_value_adjustment_method = "BH",
                        q_value_cutoff_enrichment = 0.2,
                        min_geneset_size = 10,
                        max_geneset_size = 500,
                        num_terms_enrichment_plot = 20,
                        min_subcellcomp_confidence = 1,
                        subcellcomp_show_cytosol = F,
                        min_confidence_reg_interaction = "D",
                        simplify_go = F,
                        show_ppi = T,
                        show_drugs_in_ppi = T,
                        show_disease = T,
                        show_top_diseases_only = T,
                        show_cancer_hallmarks = T,
                        show_drug = T,
                        show_enrichment = T,
                        show_tcga_aberration = T,
                        show_tcga_coexpression = T,
                        show_subcell_comp = T,
                        show_crispr_lof = T,
                        show_cell_tissue = F,
                        show_ligand_receptor = T,
                        show_regulatory_interactions = T,
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
  rep[["config"]][["resources"]] <-  oeDB$release_notes

  ## config/show - logicals indicating which sections/analyses of the report to include
  rep[["config"]][["show"]] <- list()
  rep[["config"]][["show"]][["query"]] <- TRUE
  rep[["config"]][["show"]][["ppi"]] <- show_ppi
  rep[["config"]][["show"]][["disease"]] <- show_disease
  rep[["config"]][["show"]][["drug"]] <- show_drug
  rep[["config"]][["show"]][["drug_tractability"]] <- show_drug
  rep[["config"]][["show"]][["enrichment"]] <- show_enrichment
  rep[["config"]][["show"]][["protein_complex"]] <- show_complex
  rep[["config"]][["show"]][["tcga_aberration"]] <- show_tcga_aberration
  rep[["config"]][["show"]][["tcga_coexpression"]] <- show_tcga_coexpression
  rep[["config"]][["show"]][["subcellcomp"]] <- show_subcell_comp
  rep[["config"]][["show"]][["crispr_ps"]] <- show_crispr_lof
  rep[["config"]][["show"]][["crispr_ps_fitness"]] <- show_crispr_lof
  rep[["config"]][["show"]][["crispr_ps_prioritized"]] <- show_crispr_lof
  rep[["config"]][["show"]][["cell_tissue"]] <- show_cell_tissue
  rep[["config"]][["show"]][["regulatory"]] <- show_regulatory_interactions
  rep[["config"]][["show"]][["ligand_receptor"]] <- show_ligand_receptor
  rep[["config"]][["show"]][["cancer_hallmark"]] <- show_cancer_hallmarks
  rep[["config"]][["show"]][["unknown_function"]] <- show_unknown_function
  rep[["config"]][["show"]][["cancer_prognosis"]] <-
    show_prognostic_cancer_assoc

  ## config - report metadata (project owner, background, project_title)
  rep[["config"]][["project_title"]] <- project_title
  rep[["config"]][["project_description"]] <- project_description
  rep[["config"]][["project_owner"]] <- project_owner
  rep[["config"]][["rmarkdown"]] <- list()
  rep[["config"]][["rmarkdown"]][["floating_toc"]] <- html_floating_toc
  rep[["config"]][["rmarkdown"]][["theme"]] <- html_report_theme

  ## config - query (id type and option to ignore errors)
  rep[["config"]][["query"]] <- list()
  rep[["config"]][["query"]][["id_type"]] <- query_id_type
  rep[["config"]][["query"]][["ignore_err"]] <- ignore_id_err

  rep[["config"]][["bgset"]] <- list()
  rep[["config"]][["bgset"]][["id_type"]] <- bgset_id_type

  rep[["config"]][["regulatory"]] <- list()
  rep[["config"]][["regulatory"]][["target_levels"]] <-
    c("TF_TARGET_A","TF_TARGET_B","TF_TARGET_C","TF_TARGET_D",
      "TARGET_A","TARGET_B","TARGET_C","TARGET_D")
  rep[["config"]][["regulatory"]][["target_colors"]] <-
    c("#08306b","#08519c","#2171b5", "#4292c6",
      "#08306b","#08519c","#2171b5", "#4292c6")

  rep[["config"]][["regulatory"]][["tf_levels"]] <-
    c("TF_TARGET_A","TF_TARGET_B","TF_TARGET_C","TF_TARGET_D",
      "TF_A","TF_B","TF_C","TF_D")
  rep[["config"]][["regulatory"]][["tf_colors"]] <-
    c("#08306b","#08519c","#2171b5", "#4292c6",
      "#08306b","#08519c","#2171b5", "#4292c6")

  rep[["config"]][["regulatory"]][["min_confidence"]] <-
    min_confidence_reg_interaction

  ## config - disease - color codes and
  ## thresholds for quantitative target-disease associations
  rep[["config"]][["disease"]] <- list()
  rep[["config"]][["disease"]][["breaks"]] <-
     c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  rep[["config"]][["disease"]][["colors"]] <-
    c("#b8b8ba",
      "#c6dbef","#9ecae1",
      "#6baed6","#4292c6",
      "#2171b5","#08519c",
      "#08306b")
  rep[['config']][['disease']][['show_top_diseases']] <-
    show_top_diseases_only

  rep[["config"]][["drug"]] <- list()
  rep[["config"]][["drug"]][["ab_levels_colors"]] <-
    c("#b8b8ba","#8c96c6","#8856a7","#810f7c")
  rep[["config"]][["drug"]][["ab_levels"]] <-
    c("Unknown", "Predicted_Tractable_Medium_to_low_confidence",
      "Predicted_Tractable_High_confidence","Clinical_Precedence")

  rep[["config"]][["drug"]][["sm_levels_colors"]] <-
    c("#b8b8ba","#8c96c6","#8856a7","#810f7c")
  rep[["config"]][["drug"]][["sm_levels"]] <-
    c("Unknown", "Predicted_Tractable",
      "Discovery_Precedence","Clinical_Precedence")

  ## config/loss_of_fitness - plot height for hits in CRISPR screens
  rep[["config"]][["crispr_ps"]] <- list()
  rep[["config"]][["crispr_ps"]][["plot_height_fitness"]] <- 10

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
  rep[["config"]][["enrichment"]][["num_terms_enrichment_plot"]] <-
    num_terms_enrichment_plot


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

  rep[["config"]][["unknown_function"]][['num_candidates']] <-
    oeDB[['genedb']]$all %>%
    dplyr::filter(!is.na(.data$unknown_function_rank)) %>%
    nrow()

  rep[['config']][['subcellcomp']] <- list()
  rep[['config']][['subcellcomp']][['minimum_confidence']] <-
    min_subcellcomp_confidence
  rep[['config']][['subcellcomp']][['show_cytosol']] <-
    subcellcomp_show_cytosol
  rep[['config']][['subcellcomp']][['gganatogram_legend']] <-
    oeDB$subcelldb$gganatogram_legend


  ## initialize all data elements
  for (analysis in c("query",
                     "tcga",
                     "disease",
                     "ppi",
                     "tcga",
                     "enrichment",
                     "drug",
                     "regulatory",
                     "cancer_hallmark",
                     "protein_complex",
                     "ligand_receptor",
                     "subcellcomp",
                     "crispr_ps",
                     "cell_tissue",
                     "cancer_prognosis",
                     "unknown_function")) {
    rep[["data"]][[analysis]] <- list()
  }

  ## query verification
  rep[["data"]][["query"]][["target"]] <- data.frame()
  rep[["data"]][["query"]][["validation_status"]] <- "perfect_go"

  ## prognosis/survival - gene expression (HPA)
  rep[["data"]][["cancer_prognosis"]][['hpa']] <- list()
  rep[["data"]][["cancer_prognosis"]][['hpa']][['assocs']] <- data.frame()

  ## regulatory interactions (DoRothEA)
  rep[["data"]][["regulatory"]][["interactions"]] <- list()
  rep[["data"]][["regulatory"]][["interactions"]][["global"]] <- data.frame()
  rep[["data"]][["regulatory"]][["interactions"]][["pancancer"]] <- data.frame()

  rep[["data"]][["regulatory"]][["network"]] <- list()
  rep[["data"]][["regulatory"]][["network"]][["edges"]] <- data.frame()
  rep[["data"]][["regulatory"]][["network"]][["nodes"]] <- data.frame()

  ## prognosis/survival - gene expression, mutation, CNA (CSHL)
  rep[["data"]][["cancer_prognosis"]][['km_cshl']] <- list()
  rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']] <- list()
  rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']][['cna']] <- data.frame()
  rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']][['mut']] <- data.frame()
  rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']][['exp']] <- data.frame()

  ## ligand-receptor interactions
  rep[["data"]][["ligand_receptor"]][["cell_cell_contact"]] <- data.frame()
  rep[["data"]][["ligand_receptor"]][["ecm_receptor"]] <- data.frame()
  rep[["data"]][["ligand_receptor"]][["secreted_signaling"]] <- data.frame()



  ## cancer hallmarks
  rep[["data"]][["cancer_hallmark"]][["target"]] <- data.frame()

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
  rep[["data"]][["drug"]][["target_drugs"]] <- data.frame()
  rep[["data"]][["drug"]][["tractability_sm"]] <- data.frame()
  rep[["data"]][["drug"]][["tractability_ab"]] <- data.frame()

  ## disease associations
  rep[["data"]][["disease"]][["target"]] <- data.frame()
  rep[["data"]][["disease"]][["target_assoc"]] <- data.frame()
  rep[["data"]][["disease"]][["assoc_pr_gene"]] <- list()
  rep[["data"]][["disease"]][["ttype_matrix"]] <- matrix()

  ## protein-protein interactions
  rep[["data"]][["ppi"]][["source"]] <- "STRING"
  rep[["data"]][["ppi"]][["complete_network"]] <- NULL
  rep[["data"]][["ppi"]][["hubscores"]] <- data.frame()
  rep[["data"]][["ppi"]][["community_network"]] <- NULL

  ## functional enrichment
  for (c in c("go","msigdb","wikipathways","kegg","netpath")) {
    rep[["data"]][["enrichment"]][[c]] <- data.frame()
  }

  ## unknown function
  rep[["data"]][["unknown_function"]][["hits_df"]] <- data.frame()

  ## Fitness scores and target priority scores from CRISPR/Cas9 screens
  ## (Project Score (ps))
  rep[["data"]][["crispr_ps"]][['fitness_scores']] <- list()
  rep[["data"]][["crispr_ps"]][['fitness_scores']][["targets"]] <- data.frame()
  rep[["data"]][["crispr_ps"]][['fitness_scores']][["n_targets"]] <- 0
  rep[["data"]][["crispr_ps"]][['target_priority_scores']] <- list()
  rep[["data"]][["crispr_ps"]][['target_priority_scores']][['targets']] <-
    data.frame()
  rep[["data"]][["crispr_ps"]][['target_priority_scores']][['n_pri_targets']] <-
    0


  ## TCGA co-expression
  rep[["data"]][["tcga"]][["co_expression"]] <- data.frame()

  ## Protein complexes
  rep[["data"]][["protein_complex"]][["humap2"]] <- data.frame()
  rep[["data"]][["protein_complex"]][["omnipath"]] <- data.frame()

  ## Subcellular localizations
  rep[["data"]][["subcellcomp"]][["all"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["grouped"]] <- data.frame()
  rep[["data"]][["subcellcomp"]][["anatogram"]] <- data.frame()

  ## TCGA aberrations
  rep[["data"]][["tcga"]][["recurrent_variants"]] <- data.frame()
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
      while(i <= nrow(oeDB$tcgadb$maf_codes)){
        site <- oeDB$tcgadb$maf_codes[i,]$primary_site
        code <- oeDB$tcgadb$maf_codes[i,]$code
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

#' Interrogate a gene list for cancer relevance
#'
#' @param query character vector with gene/query identifiers
#' @param oeDB oncoEnrichR annotation database - as returned from oncoEnrichR::load_db()
#' @param query_id_type character indicating source of query (one of "uniprot_acc", "symbol",
#' "entrezgene", or "ensembl_gene","ensembl_mrna","refseq_mrna","ensembl_protein","refseq_protein")
#' @param html_floating_toc logical - float the table of contents to the left of the main document content (HTML report). The floating table of contents will always be visible even when the document is scrolled
#' @param html_report_theme Bootswatch theme for HTML report (any of "bootstrap","cerulean","cosmo","default","flatly","journal","lumen","paper","sandstone","simplex","spacelab","united","yeti")
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
#' @param num_terms_enrichment_plot number of top enriched Gene Ontology terms (max) to show in enrichment barplot
#' @param min_subcellcomp_confidence minimum confidence level of subcellular compartment annotations (range from 1 to 6, 6 is strongest)
#' @param subcellcomp_show_cytosol logical indicating if subcellular heatmap should show highlight proteins located in the cytosol or not
#' @param min_confidence_reg_interaction minimum confidence level for regulatory interactions (TF-target) retrieved from DoRothEA ('A','B','C', or 'D')
#' @param simplify_go remove highly similar GO terms in results from GO enrichment/over-representation analysis
#' @param ppi_add_nodes number of nodes to add to target set when computing the protein-protein interaction network (STRING)
#' @param ppi_score_threshold minimum score (0-1000) for retrieval of protein-protein interactions (STRING)
#' @param show_ppi logical indicating if report should contain protein-protein interaction data (STRING)
#' @param show_drugs_in_ppi logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform, association_score >= 0.05, support from at least two data types)
#' @param show_top_diseases_only logical indicating if report should contain top (n = 20) disease associations only pr. query gene (Open Targets Platform)
#' @param show_cancer_hallmarks logical indicating if report should contain annotations/evidence of cancer hallmarks per query gene (COSMIC/Open Targets Platform)
#' @param show_drug logical indicating if report should contain targeted cancer drug information
#' @param show_enrichment logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME, NetPath, WikiPathways)
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of query set with oncogenes/tumor suppressor genes
#' @param show_cell_tissue logical indicating if report should contain tissue-specificity and single cell-type specificity assessments (Human Protein Atlas)
#' of target genes, using data from the Human Protein Atlas
#' @param show_ligand_receptor logical indicating if report should contain ligand-receptor interactions (CellChatDB)
#' @param show_regulatory_interactions logical indicating if report should contain data on transcription factor (TF) - target interactions relevant for the query set (DoRothEA)
#' @param show_unknown_function logical indicating if report should highlight target genes with unknown or poorly defined functions (GO/Uniprot KB/NCBI)
#' @param show_prognostic_cancer_assoc  logical indicating if mRNA-based (single-gene) prognostic associations to cancer types should be listed (Human Protein Atlas/TCGA)
#' @param show_subcell_comp logical indicating if report should provide subcellular compartment annotations (ComPPI)
#' @param show_crispr_lof logical indicating if report should provide fitness scores and target priority scores from CRISPR/Cas9 loss-of-fitness screens (Project Score)
#' @param show_complex logical indicating if report should provide target memberships in known protein complexes (ComplexPortal/Compleat/PDB/CORUM)
#' @param ... arguments for Galaxy/web-based processing
#'
#' @export
#'
onco_enrich <- function(query = NULL,
                        oeDB = NULL,
                        query_id_type = "symbol",
                        ignore_id_err = TRUE,
                        html_floating_toc = T,
                        html_report_theme = "default",
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
                        num_terms_enrichment_plot = 20,
                        min_subcellcomp_confidence = 1,
                        subcellcomp_show_cytosol = FALSE,
                        min_confidence_reg_interaction = "D",
                        simplify_go = TRUE,
                        ppi_add_nodes = 50,
                        ppi_score_threshold = 900,
                        show_ppi = TRUE,
                        show_drugs_in_ppi = TRUE,
                        show_disease = TRUE,
                        show_top_diseases_only = TRUE,
                        show_cancer_hallmarks = TRUE,
                        show_drug = TRUE,
                        show_enrichment = TRUE,
                        show_tcga_aberration = TRUE,
                        show_tcga_coexpression = TRUE,
                        show_cell_tissue = FALSE,
                        show_ligand_receptor = TRUE,
                        show_regulatory_interactions = TRUE,
                        show_unknown_function = TRUE,
                        show_prognostic_cancer_assoc = TRUE,
                        show_subcell_comp = TRUE,
                        show_crispr_lof = TRUE,
                        show_complex = TRUE,
                        ...) {


  logger <- log4r::logger(
    threshold = "INFO",
    appenders = log4r::console_appender(log4r_layout))

  dot_args <- list(...)

  if(is.null(oeDB)){
    log4r_info(logger, paste0(
      "ERROR: mandatory argument 'oeDB' cannot be NULL"))
    return()
  }
  oedb_val <- validate_db(oeDB)
  if(oedb_val != 0){
    return()
  }

  if(is.null(query)){
    log4r_info(logger, paste0(
      "ERROR: mandatory argument 'query' cannot be NULL"))
    return()
  }
  if(!is.character(query)){
    log4r_info(logger, paste0(
      "ERROR: mandatory argument 'query' is of wrong type (not character)"))
    return()
  }

  if(!is.null(bgset)){
    val <- assertthat::validate_that(
      is.character(bgset))
    if(!is.logical(val)){
      log4r_info(logger, paste0(
        "ERROR: ", val))
      return()
    }
  }

  val <- assertthat::validate_that(length(query) >= 2)
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR: query set must contain at least two entries - length of query is: ", length(query)))
    return()
  }

  val <- assertthat::validate_that(
    p_value_adjustment_method %in% c("holm", "hochberg",
                                     "hommel", "bonferroni",
                                     "BH", "BY",
                                     "fdr", "none")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'p_value_adjustment_method' must take on any of the following values: ",
      "'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'",
      " (value provided was '", p_value_adjustment_method,"')"))
    return()
  }

  val <- assertthat::validate_that(
    min_confidence_reg_interaction %in% c("A", "B",
                                          "C", "D")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'min_confidence_reg_interaction' must take on any of the following values: 'A', 'B', 'C', 'D'",
      " (value provided was '", min_confidence_reg_interaction,"')"))
    return()
  }

  val <- assertthat::validate_that(
    html_report_theme %in% c("bootstrap","cerulean","cosmo","default",
                        "flatly","journal","lumen","paper","sandstone",
                        "simplex","spacelab","united","yeti")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'html_report_theme' must take on any of the following values: ",
      "'bootstrap', 'cerulean', 'cosmo', 'default', 'flatly', 'journal', 'lumen',",
      "'paper', 'sandstone', 'simplex', 'spacelab', 'united', 'yeti'",
      " (value provided was '", p_value_adjustment_method,"')"))
    return()
  }

  ## Number of allowed query genes
  oncoenrichr_query_limit <- 600

  if(length(names(dot_args)) > 0){
    if("galaxy" %in% names(dot_args))
      if(is.logical(dot_args$galaxy)){
        if(dot_args$galaxy == T){
          oncoenrichr_query_limit <- 500
        }
      }
  }

  if (length(query) > oncoenrichr_query_limit) {
    log4r_info(logger,
      paste0("WARNING: oncoEnrichR is limited to the analysis of n = ",
             oncoenrichr_query_limit," entries. Query contained n = ",
             length(query), " identifiers, limiting to top ",
             oncoenrichr_query_limit," genes"))
    query <- utils::head(unique(query),
                         oncoenrichr_query_limit)
  }

  val <- assertthat::validate_that(
    query_id_type %in% c("symbol", "entrezgene", "refseq_mrna", "ensembl_mrna",
                         "refseq_protein", "ensembl_protein", "uniprot_acc",
                         "ensembl_gene")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'query_id_type' must take on of the following values: ",
      "'symbol', 'entrezgene', 'refseq_mrna', 'ensembl_mrna', 'flatly'",
      "'refseq_protein', 'ensembl_protein', 'uniprot_acc', 'ensembl_gene'",
      " (value provided was '", query_id_type,"')"))
    return()
  }

  val <- assertthat::validate_that(
    bgset_id_type %in% c("symbol", "entrezgene", "refseq_mrna", "ensembl_mrna",
                         "refseq_protein", "ensembl_protein", "uniprot_acc",
                         "ensembl_gene")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'bgset_id_type' must take on any of the following values: ",
      "'symbol', 'entrezgene', 'refseq_mrna', 'ensembl_mrna', 'flatly'",
      "'refseq_protein', 'ensembl_protein', 'uniprot_acc', 'ensembl_gene'",
      " (value provided was '", bgset_id_type,"')"))
    return()
  }

  val <-
    (ppi_score_threshold %% 1 == 0) & ## check that number is whole integer
      (ppi_score_threshold > 0) &
      (ppi_score_threshold <= 1000)

  if(val == F){
    log4r_info(logger, paste0(
      "ERROR - 'ppi_score_threshold' must be an integer/whole number and take a value from 1 to 1000 ",
      "(current type and value: '",typeof(ppi_score_threshold),"' - ",
      ppi_score_threshold,")")
    )
    return()
  }
  val <-
    num_terms_enrichment_plot %% 1 == 0 &
      num_terms_enrichment_plot >= 10 &
      num_terms_enrichment_plot <= 30

  if(val == F){
    log4r_info(logger, paste0(
      "ERROR - 'num_terms_enrichment_plot' must be an integer/whole number and take a value from 10 to 30 ",
      "(current type and value: '",typeof(num_terms_enrichment_plot),"' - ",
      num_terms_enrichment_plot,")")
    )
    return()
  }

  val <-
    min_subcellcomp_confidence %% 1 == 0 &
      min_subcellcomp_confidence >= 1 &
      min_subcellcomp_confidence <= 6

  if(val == F){
    log4r_info(logger, paste0(
      "ERROR - 'min_subcellcomp_confidence' must be an integer/whole number and take a value from 1 to 6 ",
      "(current type and value: '",typeof(min_subcellcomp_confidence),"' - ",
      min_subcellcomp_confidence,")")
    )
    return()
  }

  val <- assertthat::validate_that(
    is.numeric(p_value_cutoff_enrichment) &
      p_value_cutoff_enrichment > 0 & p_value_cutoff_enrichment < 1)

  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'p_value_cutoff_enrichment' must be of type numeric and be greater than 0 and less than 1 ",
      "(current type and value: '",typeof(p_value_cutoff_enrichment),"' - ",
      p_value_cutoff_enrichment,")")
    )
  }

  val <- assertthat::validate_that(
    is.numeric(q_value_cutoff_enrichment) &
      q_value_cutoff_enrichment > 0 & q_value_cutoff_enrichment < 1)

  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR - 'q_value_cutoff_enrichment' must be of type numeric and be greater than 0 and less than 1 ",
      "(current type and value: '",typeof(q_value_cutoff_enrichment),"' - ",
      q_value_cutoff_enrichment,")")
    )
  }


  stopifnot(ppi_add_nodes <= 50)



  ## Initialize the oncoEnrichR report structure
  onc_rep <- init_report(
    oeDB = oeDB,
    project_title = project_title,
    project_owner = project_owner,
    html_report_theme = html_report_theme,
    html_floating_toc = html_floating_toc,
    query_id_type = query_id_type,
    ignore_id_err = ignore_id_err,
    project_description = project_description,
    ppi_min_string_score = ppi_score_threshold,
    ppi_add_nodes = ppi_add_nodes,
    bgset_description =
      bgset_description,
    bgset_id_type = bgset_id_type,
    p_value_cutoff_enrichment = p_value_cutoff_enrichment,
    p_value_adjustment_method = p_value_adjustment_method,
    q_value_cutoff_enrichment = q_value_cutoff_enrichment,
    min_geneset_size = min_geneset_size,
    max_geneset_size = max_geneset_size,
    num_terms_enrichment_plot = num_terms_enrichment_plot,
    min_subcellcomp_confidence = min_subcellcomp_confidence,
    subcellcomp_show_cytosol = subcellcomp_show_cytosol,
    min_confidence_reg_interaction = min_confidence_reg_interaction,
    simplify_go = simplify_go,
    show_ppi = show_ppi,
    show_drugs_in_ppi = show_drugs_in_ppi,
    show_disease = show_disease,
    show_top_diseases_only = show_top_diseases_only,
    show_cancer_hallmarks = show_cancer_hallmarks,
    show_drug = show_drug,
    show_enrichment = show_enrichment,
    show_tcga_aberration = show_tcga_aberration,
    show_tcga_coexpression = show_tcga_coexpression,
    show_subcell_comp = show_subcell_comp,
    show_crispr_lof = show_crispr_lof,
    show_cell_tissue = show_cell_tissue,
    show_ligand_receptor = show_ligand_receptor,
    show_regulatory_interactions = show_regulatory_interactions,
    show_unknown_function = show_unknown_function,
    show_prognostic_cancer_assoc =
      show_prognostic_cancer_assoc,
    show_complex = show_complex)

  ## validate query gene set
  qgenes_match <-
    validate_query_genes(
      qgenes = query,
      q_id_type = query_id_type,
      ignore_id_err = ignore_id_err,
      genedb = oeDB[['genedb']][['all']],
      transcript_xref_db = oeDB[['genedb']][['transcript_xref']],
      logger = logger)

  ## assign validation result to report object
  onc_rep[['data']][['query']][['target']] <-
    qgenes_match[['all']]
  onc_rep[["data"]][["query"]][["validation_status"]] <-
    qgenes_match[["match_status"]]

  ## if query list contains non-validated entries and 'ignore_id_err' is
  ## set to FALSE, oncoEnrichR analysis will halt (no analysis modules
  ## will be included in report)
  if (qgenes_match[["match_status"]] == "imperfect_stop") {
    for(e in c("disease",
               "drug",
               "enrichment",
               "protein_complex",
               "ppi",
               "tcga_aberration",
               "tcga_coexpression",
               "cell_tissue",
               "cancer_hallmark",
               "crispr_ps",
               "regulatory_interactions",
               "ligand_receptor",
               "subcellcomp",
               "unknown_function",
               "cancer_prognosis")){
      onc_rep[['config']][['show']][[e]] <- F
    }
    return(onc_rep)
  }

  ## If list of genes are less than 5, pathway and GO enrichment
  ## is not performed (too few genes)
  if(NROW(qgenes_match[["found"]]) < 5){
    onc_rep[['config']][['show']][['enrichment']] <- F
    log4r_info(logger,
      paste0("WARNING: Function (GO) and pathway enrichment is NOT performed for gene sets of size < 5. Query contained n = ",
             NROW(qgenes_match[["found"]])," valid entries"))
  }


  ## validate background gene set (if provided)
  background_entrez <- NULL
  background_genes_match <- NULL
  if (!is.null(bgset)) {
    background_genes_match <-
      validate_query_genes(
        bgset,
        q_id_type = bgset_id_type,
        genedb = oeDB[['genedb']][['all']],
        transcript_xref_db = oeDB[['genedb']][['transcript_xref_db']],
        qtype = "background",
        logger = logger)
    if (background_genes_match[["match_status"]] == "imperfect_stop") {
      log4r_info(logger, paste0("WARNING: Background geneset not defined properly - ",
                        "using all protein-coding genes instead"))
      background_entrez <- unique(oeDB[['genedb']][['all']]$entrezgene)
    }else{
      background_entrez <- unique(background_genes_match[["found"]]$entrezgene)
    }
  }



  query_entrezgene <- unique(qgenes_match[["found"]]$entrezgene)
  query_symbol <- unique(qgenes_match[["found"]]$symbol)

  if (length(query_symbol) > 20) {
    onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] <-
      onc_rep[["config"]][["tcga_aberration"]][["plot_height"]] +
      as.integer((length(query_symbol) - 20)/ 7.5)
  }

  ## Include gene-disease annotations in the report, and
  ## tumor-type specific rankings
  if (show_disease == T) {
    onc_rep[["data"]][["disease"]] <-
      target_disease_associations(
        qgenes = query_symbol,
        show_top_diseases_only = show_top_diseases_only,
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)
  }

  ## Include gene-drug annotations in the report (targeted cancer drugs)
  if (show_drug == T) {
    onc_rep[["data"]][["drug"]] <-
      target_drug_associations(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        cancerdrugdb = oeDB[['cancerdrugdb']],
        logger = logger)
  }

  ## Include ligand-receptor interactions in the report
  if (show_ligand_receptor == T) {
    onc_rep[["data"]][["ligand_receptor"]] <-
      annotate_ligand_receptor_interactions(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        ligandreceptordb = oeDB[['ligandreceptordb']],
        logger = logger)
  }

  ## Include enrichment analyses in the report (pathway, GO, MSigDB)
  if(onc_rep[['config']][['show']][['enrichment']] == T){
    for (c in names(oeDB[['pathwaydb']][["msigdb"]][["COLLECTION"]])) {
      for (subcat in names(oeDB[['pathwaydb']][["msigdb"]][["COLLECTION"]][[c]])) {
        if (c == "C5" & subcat != "HPO") {
          enr <- get_go_enrichment(
            query_entrez = as.character(query_entrezgene),
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
            genedb = oeDB[['genedb']][['all']],
            logger = logger)
          if (!is.null(enr)) {
            onc_rep[["data"]][["enrichment"]][["go"]] <-
              dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["go"]], enr) %>%
              dplyr::distinct()
          }
        }else{
          if(c != "C5"){
            db = paste0("MSIGdb/", c, "/",subcat)
            enr <- get_universal_enrichment(
              query_entrez = as.character(query_entrezgene),
              genedb = oeDB[['genedb']][['all']],
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
              TERM2GENE = oeDB[['pathwaydb']][['msigdb']]$COLLECTION[[c]][[subcat]]$TERM2GENE,
              TERM2NAME = oeDB[['pathwaydb']][['msigdb']]$COLLECTION[[c]][[subcat]]$TERM2NAME,
              TERM2SOURCE = oeDB[['pathwaydb']][['msigdb']]$TERM2SOURCE,
              dbsource = db,
              logger = logger)
            if (!is.null(enr)) {
              onc_rep[["data"]][["enrichment"]][["msigdb"]] <-
                dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["msigdb"]], enr) %>%
                dplyr::distinct() %>%
                dplyr::arrange(.data$qvalue)
            }
          }
        }
      }
    }

    for(pwaydb in c('wikipathways','netpath','kegg')){
      db <- "NetPath"
      if(pwaydb == "wikipathways"){
        db <- "WikiPathways"
      }
      if(pwaydb == "kegg"){
        db <- "KEGG"
      }
      onc_rep[["data"]][["enrichment"]][[pwaydb]] <-
        get_universal_enrichment(
          as.character(query_entrezgene),
          genedb = oeDB[['genedb']][['all']],
          background_entrez = background_entrez,
          min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
          max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
          q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
          p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
          p_value_adjustment_method =
            onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
          TERM2GENE = oeDB[['pathwaydb']][[pwaydb]]$TERM2GENE,
          TERM2NAME = oeDB[['pathwaydb']][[pwaydb]]$TERM2NAME,
          TERM2SOURCE = oeDB[['pathwaydb']][[pwaydb]]$TERM2SOURCE,
          dbsource = db,
          logger = logger)

    }
  }

  ## Include protein-protein interactions in the report
  if (show_ppi == T) {
    ## TEST TO CHECK OMNIPATHDB IS LIVE IS NOT WORKING
    # service_is_down <- unique(is.na(pingr::ping("string-db.org")))
    onc_rep[["data"]][["ppi"]] <-
      get_ppi_network(
        qgenes = query_entrezgene,
        ppi_source = "STRING",
        genedb = oeDB[['genedb']][['all']],
        cancerdrugdb = oeDB[['cancerdrugdb']],
        settings = onc_rep[["config"]][["ppi"]][["stringdb"]],
        logger = logger)
  }

  ## include protein complex annotations in the report
  if (show_complex == T) {
    onc_rep[["data"]][["protein_complex"]] <-
      annotate_protein_complex(
        query_entrez = as.integer(query_entrezgene),
        genedb = oeDB[['genedb']][['all']],
        complex_db = oeDB[['genedb']][['proteincomplexdb']],
        transcript_xref_db = oeDB[['genedb']][['transcript_xref']],
        logger = logger)
  }

  ## include annotations regarding genes of unknown/poorly defined function
  if (show_unknown_function == T) {
    onc_rep[["data"]][["unknown_function"]][["hits_df"]] <-
      get_genes_unknown_function(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        logger = logger
      )
  }

  ## include subcellular enrichment and annotations
  if (show_subcell_comp == T) {
     subcellcomp_annotations <-
      annotate_subcellular_compartments(
        query_entrez = as.integer(query_entrezgene),
        minimum_confidence = min_subcellcomp_confidence,
        show_cytosol = subcellcomp_show_cytosol,
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

     onc_rep[["data"]][["subcellcomp"]][["all"]] <-
       subcellcomp_annotations[["all"]]
     onc_rep[["data"]][["subcellcomp"]][["grouped"]] <-
       subcellcomp_annotations[["grouped"]]
     onc_rep[["data"]][["subcellcomp"]][["anatogram"]] <-
       subcellcomp_annotations[["anatogram"]]

  }

  ##
  if (show_crispr_lof == T) {
    onc_rep[["data"]][["crispr_ps"]][["fitness_scores"]] <-
      get_crispr_lof_scores(
        qgenes = query_symbol,
        projectscoredb = oeDB[['projectscoredb']],
        logger = logger)

    if (onc_rep[["data"]][["crispr_ps"]][["fitness_scores"]][["n_targets"]] <= 10){
      onc_rep[["config"]][["crispr_ps"]][["plot_height_fitness"]] <- 5
    }

    if (onc_rep[["data"]][["crispr_ps"]][["fitness_scores"]][["n_targets"]] >= 20) {
      onc_rep[["config"]][["crispr_ps"]][["plot_height_fitness"]]  <-
        onc_rep[["config"]][["crispr_ps"]][["plot_height_fitness"]] +
        as.integer((onc_rep[["data"]][["crispr_ps"]][["fitness_scores"]][["n_targets"]] - 20)/8.5)
    }

    onc_rep[["data"]][["crispr_ps"]][["target_priority_scores"]] <-
      get_target_priority_scores(
        qgenes = query_symbol,
        projectscoredb = oeDB[['projectscoredb']],
        logger = logger)
  }


  if (show_tcga_aberration == T) {
    for (v in c("cna_homdel","cna_ampl")) {
      onc_rep[["data"]][["tcga"]][["aberration"]][["matrix"]][[v]] <-
        tcga_aberration_matrix(
          qgenes = as.integer(query_entrezgene),
          qsource = "entrezgene",
          genedb = oeDB[['genedb']][['all']],
          oeDB = oeDB,
          vtype = v,
          logger = logger)
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <-
        tcga_aberration_table(
          qgenes = as.integer(query_entrezgene),
          qsource = "entrezgene",
          genedb = oeDB[['genedb']][['all']],
          oeDB = oeDB,
          vtype = v,
          logger = logger)
    }

    onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
      oeDB$tcgadb[["recurrent_variants"]] %>%
      dplyr::inner_join(
        dplyr::select(qgenes_match$found, .data$symbol),
        by = c("SYMBOL" = "symbol")) %>%
      dplyr::distinct()

    if(nrow(onc_rep[["data"]][["tcga"]][["recurrent_variants"]]) > 0){
      cosmic_variants <-
        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] %>%
        dplyr::select(.data$VAR_ID, .data$COSMIC_MUTATION_ID) %>%
        dplyr::filter(!is.na(.data$COSMIC_MUTATION_ID)) %>%
        dplyr::distinct()

      if(nrow(cosmic_variants) > 0){

        cosmic_variants <- as.data.frame(
          cosmic_variants %>%
          tidyr::separate_rows(.data$COSMIC_MUTATION_ID, sep="&") %>%
          dplyr::mutate(
            COSMIC_MUTATION_ID = paste0(
              "<a href=\"https://cancer.sanger.ac.uk/cosmic/search?q=",
              .data$COSMIC_MUTATION_ID,"\" target='_blank'>",
              .data$COSMIC_MUTATION_ID,"</a>"
            )) %>%
          dplyr::group_by(.data$VAR_ID) %>%
          dplyr::summarise(
            COSMIC_MUTATION_ID =
              paste(
                .data$COSMIC_MUTATION_ID, collapse=", "
              ),
            .groups = "drop")
        )

        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
          onc_rep[["data"]][["tcga"]][["recurrent_variants"]] %>%
          dplyr::select(-.data$COSMIC_MUTATION_ID) %>%
          dplyr::left_join(cosmic_variants,
                           by = c("VAR_ID"))
      }

      onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] %>%
        dplyr::left_join(oeDB[['tcgadb']][['pfam']], by = "PFAM_ID") %>%
        dplyr::mutate(
          PROTEIN_DOMAIN = dplyr::if_else(
            !is.na(.data$PFAM_ID),
            paste0(
            "<a href=\"http://pfam.xfam.org/family/", .data$PFAM_ID,
            "\" target='_blank'>",
            .data$PFAM_DOMAIN_NAME,
            "</a>"),
            as.character(NA)
          )
        ) %>%
        dplyr::select(-c(.data$PFAM_DOMAIN_NAME, .data$PFAM_ID)) %>%
        dplyr::left_join(dplyr::select(oeDB[['genedb']][['all']],
                                       .data$symbol, .data$ensembl_gene_id),
                         by = c("SYMBOL" = "symbol")) %>%
        dplyr::rename(ENSEMBL_GENE_ID = .data$ensembl_gene_id) %>%
        dplyr::mutate(ENSEMBL_TRANSCRIPT_ID =
                        paste0("<a href='https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",
                               .data$ENSEMBL_GENE_ID,
                               ";t=",
                               .data$ENSEMBL_TRANSCRIPT_ID,"' target='_blank'>",
                               .data$ENSEMBL_TRANSCRIPT_ID,"</a>")) %>%
        dplyr::select(-.data$VAR_ID) %>%
        dplyr::select(.data$SYMBOL,
                      .data$CONSEQUENCE,
                      .data$PROTEIN_CHANGE,
                      .data$MUTATION_HOTSPOT,
                      .data$PROTEIN_DOMAIN,
                      .data$LOSS_OF_FUNCTION,
                      .data$ENSEMBL_GENE_ID,
                      .data$ENSEMBL_TRANSCRIPT_ID,
                      dplyr::everything())
    }

    for(psite in names(onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]])){
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]][[psite]][['top_mutated_genes']] <-
        tcga_oncoplot_genes(
          qgenes = query_symbol,
          qsource = "symbol",
          genedb = oeDB[['genedb']][['all']],
          oeDB = oeDB,
          site = psite,
          logger = logger)
    }
  }

  if (show_cancer_hallmarks == T){

    log4r_info(logger, "Retrieving genes with evidence of cancer hallmark properties")
    onc_rep[["data"]][["cancer_hallmark"]][["target"]] <-
      oeDB[['genedb']][["cancer_hallmark"]][["short"]] %>%
      dplyr::inner_join(
        dplyr::select(qgenes_match$found, .data$ensembl_gene_id),
        by = "ensembl_gene_id") %>%
      dplyr::select(-.data$ensembl_gene_id) %>%
      dplyr::distinct()

    n_genes_cancerhallmarks <- 0

    if(nrow(onc_rep[["data"]][["cancer_hallmark"]][["target"]]) > 0){
      n_genes_cancerhallmarks <- length(unique(onc_rep[["data"]][["cancer_hallmark"]][["target"]]$symbol))
    }
    log4r_info(logger, paste0("Number of query genes attributed with cancer hallmark properties: ",
                             n_genes_cancerhallmarks))

  }

  if (show_tcga_coexpression == T) {
    onc_rep[["data"]][["tcga"]][["co_expression"]] <-
      tcga_co_expression(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)
  }

  if (show_regulatory_interactions == T) {

    for(collection in c('global','pancancer')){
      onc_rep[["data"]][["regulatory"]][["interactions"]][[collection]] <-
        annotate_tf_targets(
          query_symbol,
          genedb = oeDB[['genedb']][['all']],
          tf_target_interactions = oeDB[['tftargetdb']],
          collection = collection,
          min_confidence_reg_interaction = min_confidence_reg_interaction,
          logger = logger)
    }

    if(NROW(onc_rep[["data"]][["regulatory"]][["interactions"]][["pancancer"]]) > 0){
      onc_rep[["data"]][["regulatory"]][["network"]] <-
        retrieve_tf_target_network(
          tf_target_interactions =
            onc_rep[["data"]][["regulatory"]][["interactions"]][["pancancer"]]
        )
    }
  }

  if(show_prognostic_cancer_assoc == T){
    onc_rep[["data"]][["cancer_prognosis"]][['hpa']][['assocs']] <-
      hpa_prognostic_genes(
        query_symbol,
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

    for(feature in c('exp', 'mut', 'cna', 'meth')){
      onc_rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']][[feature]] <-
        km_cshl_survival_genes(
          query_symbol,
          projectsurvivaldb = oeDB[['projectsurvivaldb']][[feature]],
          logger = logger,
          genetic_feature = feature)
    }

  }

  if(show_cell_tissue == T){
    onc_rep[["data"]][["cell_tissue"]][['tissue_overview']] <-
      gene_tissue_cell_spec_cat(
        qgenes = query_symbol,
        q_id_type = "symbol",
        resolution = "tissue",
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

    onc_rep[["data"]][["cell_tissue"]][['tissue_enrichment']] <-
      gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(query_entrezgene),
        resolution = "tissue",
        background_entrez = background_entrez,
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

    onc_rep[["data"]][["cell_tissue"]][['scRNA_overview']] <-
      gene_tissue_cell_spec_cat(
        qgenes = as.integer(query_entrezgene),
        q_id_type = "entrezgene",
        resolution = "single_cell",
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

    onc_rep[["data"]][["cell_tissue"]][['scRNA_enrichment']] <-
      gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(query_entrezgene),
        background_entrez = background_entrez,
        resolution = "single_cell",
        genedb = oeDB[['genedb']][['all']],
        oeDB = oeDB,
        logger = logger)

  }

  return(onc_rep)

}

#' Function that writes an oncoEnrichR report object to file
#'
#' @param report object with oncoEnrichR report data (returned by oeDB$onco_enrich)
#' @param file full filename for report output (e.g. "oe_report.html", "oe_report.xlsx")
#' @param ignore_file_extension logical to accept any type of filaname extensions (for Galaxy integration)
#' @param overwrite logical indicating if existing output files may be overwritten
#' @param format file format of output (html/excel)
#' @param ... options for Galaxy/non self-contained HTML. Only applicable for use in Galaxy
#'
#' @export

write <- function(report,
                  file = "testReport.html",
                  ignore_file_extension = F,
                  overwrite = F,
                  format = "html",
                  ...) {

  selfcontained <- T
  galaxy_run <- T
  html_extern_path <- NA
  logger <- log4r::logger(threshold = "INFO",
                          appenders = log4r::console_appender(log4r_layout))
  dot_args <- list(...)
  if(length(names(dot_args)) > 0){

    for(arg in names(dot_args)){
      if(!(arg == "selfcontained_html" | arg == "galaxy" | arg == "extra_files_path")){
        log4r_info(logger, paste0(
          "ERROR: argument '",arg,"' does not exist for oeDB$write()"))
        return()
      }
    }

    if("selfcontained_html" %in% names(dot_args)){
      if(is.logical(dot_args$selfcontained_html)){
        selfcontained <- dot_args$selfcontained_html
      }
    }
    if("galaxy" %in% names(dot_args)){
      if(is.logical(dot_args$galaxy)){
        galaxy_run <- dot_args$galaxy
      }
    }
    if("extra_files_path" %in% names(dot_args)){
      if(is.character(dot_args$extra_files_path)){
        html_extern_path <- dot_args$extra_files_path

        if(!dir.exists(html_extern_path)){
          system(paste0('mkdir ', html_extern_path))
        }
        # val <- assertthat::validate_that(
        #   dir.exists(html_extern_path)
        # )
        # if(!is.logical(val)){
        #   log4r_info(paste0("ERROR: ",val))
        #   return()
        # }
      }
    }
  }

  val <- assertthat::validate_that(
    format %in% c("html","excel")
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR: ",val))
    return()
  }

  val <- assertthat::validate_that(
    is.character(file)
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR: ",val))
    return()
  }

  val <- assertthat::validate_that(
    is.logical(overwrite)
  )
  if(!is.logical(val)){
    log4r_info(logger, paste0(
      "ERROR: ",val))
    return()
  }

  output_directory <- dirname(file)
  file_basename <- basename(file)
  file_basename_prefix <- stringr::str_replace(
    file_basename,"\\.(html|xlsx)$","")

  if(output_directory != "."){
    val <- assertthat::validate_that(
      dir.exists(output_directory)
    )
    if(!is.logical(val)){
      log4r_info(paste0("ERROR: ",val))
      return()
    }
  }

  if(overwrite == F){
    val <- assertthat::validate_that(
      file.exists(file)
    )
    if(is.logical(val)){
      log4r_info(logger,
        paste0("ERROR: output file (",file, ") exists, (overwrite = F)"))
      return()
    }
  }

  if(ignore_file_extension == F){
    if(format == "html" & tools::file_ext(file) != "html"){
      log4r_info(logger, paste0(
        "ERROR: oncoEnrichR HTML output: File name must end with .html, not ",
        tools::file_ext(file))
      )
      return()
    }

    if(format == "excel" & tools::file_ext(file) != "xlsx"){
      log4r_info(logger, paste0(
        "ERROR: oncoEnrichR Excel output: File name must end with .xlsx, not ",
        tools::file_ext(file))
      )
      return()
    }
  }

  ## TODO: check that report parameter is a valid oncoEnrichR result object
  if(!is.null(report)){
    assign("onc_enrich_report",
           report, envir = .GlobalEnv)
  }else{
    log4r_info(logger,
      "ERROR: oncoEnrichR report object is NULL - cannot write report contents")
    return()
  }


  if (format == "html") {

    log4r_info(logger, "------")
    log4r_info(logger, "Writing HTML file with report contents")

    if(selfcontained == F){

      oe_rmarkdown_template_dir <-
        system.file("templates", package = "oncoEnrichR")
      tmpdir <-
        file.path(
          output_directory,
          paste0("tmp_",
                 stringi::stri_rand_strings(
                   1, 20, pattern = "[A-Za-z0-9]")
          )
        )
      system(paste0('mkdir ', tmpdir))
      system(paste0('cp ', oe_rmarkdown_template_dir, "/* ",
                    tmpdir))

      floating_toc <- "false"
      if(onc_enrich_report$config$rmarkdown$floating_toc == T){
        floating_toc <- "true"
      }

      sink(file = file.path(tmpdir,"_site.yml"))
      cat("output:",
          "  html_document:",
          "    number_sections: false",
          "    toc: true",
          "    fig_width: 5",
          "    fig_height: 4",
          "    toc_depth: 3",
          paste0("    toc_float: ", floating_toc),
          "    highlight: null",
          "    mathjax: null",
          paste0("    theme: ", onc_enrich_report$config$rmarkdown$theme),
          "    include:",
          "      after_body: _disclaimer.md", sep="\n")
      sink()


      rmdown_html <- file.path(tmpdir,"_site", "index.html")
      rmdown_supporting1 <- file.path(tmpdir,"_site","index_files")
      rmdown_supporting2 <- file.path(tmpdir,"_site","site_libs")

      if(galaxy_run == T){

        if(dir.exists(file.path(
          html_extern_path, "site_libs"))){
          if(overwrite == F){
            log4r_info(logger, paste0(
              "ERROR: Cannot create HTML since 'site_libs' exist in output_directory",
              " and 'overwrite' is FALSE")
            )
            return()
          }else{
            system(paste0('rm -rf ',file.path(
              html_extern_path, "site_libs"
            )))
          }
        }
        if(dir.exists(file.path(
          html_extern_path, "index_files"))){
          if(overwrite == F){
            log4r_info(logger, paste0(
              "ERROR: Cannot create HTML since 'index_files' exist in output_directory",
              " and 'overwrite' is FALSE")
            )
            return()
          }
          else{
            system(paste0('rm -rf ',file.path(
              html_extern_path, "index_files"
            )))
          }
        }

        rmarkdown::render_site(
          input = tmpdir,
          quiet = T
        )

        # target_html <- file.path(output_directory, paste0(
        #   file_basename_prefix, ".html")
        # )

        if(file.exists(rmdown_html) & dir.exists(rmdown_supporting1) &
           dir.exists(rmdown_supporting2)){
          system(paste0('mv ', rmdown_html, ' ',
                        file))
          system(paste0('mv ', rmdown_supporting1," ", html_extern_path))
          system(paste0('mv ', rmdown_supporting2," ", html_extern_path))
          system(paste0('rm -rf ',tmpdir))

          log4r_info(logger, paste0("Output file: ",
                                          file))
          log4r_info(logger, "------")
        }
      }

    }else{

      disclaimer <- system.file("templates",
                                "_disclaimer.md",
                                package = "oncoEnrichR")

      report_theme <- onc_enrich_report$config$rmarkdown$theme
      toc_float <- onc_enrich_report$config$rmarkdown$floating_toc

      markdown_input <- system.file("templates", "index.Rmd",
                                    package = "oncoEnrichR")

      rmarkdown::render(
        markdown_input,
        output_format = rmarkdown::html_document(
          theme = report_theme,
          toc = T,
          fig_width = 5,
          highlight = NULL,
          mathjax = NULL,
          fig_height = 4,
          toc_depth = 3,
          toc_float = toc_float,
          number_sections = F,
          includes = rmarkdown::includes(after_body = disclaimer)),
        output_file = file_basename,
        output_dir = output_directory,
        clean = T,
        intermediates_dir = output_directory,
        quiet = T)

      log4r_info(logger, paste0("Output file (self-contained HTML): ",
                                      file))
      log4r_info(logger, "------")
    }


  }
  if (format == "excel") {

    wb <- openxlsx::createWorkbook()
    log4r_info(logger, "------")
    log4r_info(logger, "Writing Excel workbook with report contents")

    table_style_index <- 15
    for(elem in c("query",
                  "disease_association",
                  "unknown_function",
                  "cancer_hallmark",
                  "drug_known",
                  "drug_tractability",
                  "protein_complex",
                  "enrichment",
                  "regulatory",
                  "subcellcomp",
                  "cell_tissue",
                  "ppi",
                  "ligand_receptor",
                  "aberration",
                  "coexpression",
                  "prognostic_association",
                  "survival_association",
                  "crispr_ps_fitness",
                  "crispr_ps_prioritized"
                  )){

      show_elem <- elem
      if(elem == "disease_association"){
        show_elem <- "disease"
      }
      if(elem == "prognostic_association"){
        show_elem <- "cancer_prognosis"
      }
      if(elem == "survival_association"){
        show_elem <- "cancer_prognosis"
      }
      if(elem == "drug_known"){
        show_elem <- "drug"
      }
      if(elem == "aberration"){
        show_elem <- "tcga_aberration"
      }
      if(elem == "coexpression"){
        show_elem <- "tcga_coexpression"
      }

      if(report[['config']][['show']][[show_elem]] == FALSE){
        next
      }

      wb <- add_excel_sheet(
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
    log4r_info(logger, paste0("Output file: ",file))
    log4r_info(logger, "------")
  }
  else{
    if(format == "json"){
      log4r_info(logger, "JSON output not yet implemented")
    }
  }
}
