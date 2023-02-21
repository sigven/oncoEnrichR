#' Load oncoEnrichR data repository
#'
#' Function that fetches version-tagged oncoEnrichR data
#' repository from Google Drive
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should force downloaded
#' (i.e. set to TRUE to re-download even if data exists in cache)
#'
#' @returns
#' A `list` object with oncoEnrichR datasets, to be used as
#' the \emph{oeDB} argument for `onco_enrich()` and `write()`
#'
#' @export
#'
#'
load_db <- function(cache_dir = NA,
                    force_download = F) {


  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  lgr::lgr$info( paste0("Loading oncoEnrichR annotation databases"))

  if (is.na(cache_dir)) {
    lgr::lgr$fatal(paste0("Argument cache_dir = '",
                          cache_dir, "' is not defined"))
    stop()
  }

  if (!dir.exists(cache_dir)) {
    lgr::lgr$fatal(paste0("Argument cache_dir = '",
                          cache_dir, "' does not exist"))
    stop()
  }

  oe_version <- unique(db_id_ref$pVersion)

  cache_dir_ext <-
    file.path(cache_dir, ".oncoenrichr")
  if (!dir.exists(cache_dir_ext)) {
    dir.create(cache_dir_ext)
  }
  cache_version_dir <-
    file.path(cache_dir,
              ".oncoenrichr",
              paste0("v", oe_version))

  if (!dir.exists(cache_version_dir)) {
    dir.create(cache_version_dir)
    lgr::lgr$info( paste0("Data will be cached in ", cache_version_dir))
  }

  oedb <- list()
  drug_datasets <- list()
  oe_datasets <- c("cancerdrugdb",
                   "genedb",
                   "biogrid",
                   "hpa",
                   "ligandreceptordb",
                   "otdb",
                   "pfamdb",
                   "pathwaydb",
                   "depmapdb",
                   "survivaldb",
                   "release_notes",
                   "subcelldb",
                   "slparalogdb",
                   "tcgadb",
                   "tissuecelldb",
                   "tftargetdb")

  for (elem in oe_datasets) {

    fname_local <- file.path(
      cache_version_dir,
      db_id_ref[db_id_ref$name == elem,]$filename
    )

    fname_gd <- googledrive::as_id(
      db_id_ref[db_id_ref$name == elem,]$gid)

    md5checksum_package <-
      db_id_ref[db_id_ref$name == elem,]$md5Checksum

    #dat <- NULL
    if (file.exists(fname_local) & force_download == F) {
      cat(fname_local, '\n')
      oedb[[elem]] <- readRDS(fname_local)
      if (!is.null(oedb[[elem]])) {
        lgr::lgr$info(paste0(
          "Reading from cache_dir = '", cache_version_dir, "', argument force_download = F"))
        lgr::lgr$info(paste0("Object '",elem,"' sucessfully loaded"))
      }

    } else {

      googledrive::drive_deauth()

      lgr::lgr$info("Downloading remote oncoEnrichR dataset from Google Drive to cache_dir")
      dl <- googledrive::with_drive_quiet(
        googledrive::drive_download(
          fname_gd,
          path = fname_local,
          overwrite = TRUE)
      )

      md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
      md5checksum_local <- tools::md5sum(fname_local)
      names(md5checksum_local) <- NULL

      if (md5checksum_remote == md5checksum_local) {
        oedb[[elem]] <- readRDS(fname_local)
        if (!is.null(oedb[[elem]])) {
          lgr::lgr$info(paste0(
            "Reading from cache_dir = ' (", cache_version_dir, "'), argument force_download = F"))
          lgr::lgr$info(paste0("Object '", elem, "' sucessfully loaded"))
          lgr::lgr$info(paste0("md5 checksum is valid: ", md5checksum_remote))
        }
      } else {
        lgr::lgr$error(paste0("md5 checksum of local file (", md5checksum_local,
                              ") is inconsistent with remote file (",
                              md5checksum_remote,")"))
        stop()
      }

    }
  }

  return(oedb)
}



#' Load oncoEnrichR annotation database
#'
#' @param remote If TRUE (default), load database from version-tagged Zenodo repository
#' @param cache_dir path to local cache for more efficient subsequent loads
#'
#' @keywords internal
#'
# load_db_zenodo <- function(remote = T,
#                            cache_dir = NA) {
#
#   lgr::lgr$appenders$console$set_layout(
#     lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
#
#   lgr::lgr$info( paste0("Loading oncoEnrichR annotation databases"))
#
#   ## check that either remote is TRUE or cache_dir is provided
#   val <- remote == T | !is.na(cache_dir)
#   if (val == F) {
#     lgr::lgr$info(
#       "ERROR: Pull database either remotely from zenodo.org ('remote' = T), or provide a cache directory ('cache_dir') with pre-loaded data")
#       return(-1)
#   }
#
#   if (!is.na(cache_dir)) {
#     val <- assertthat::validate_that(
#       dir.exists(cache_dir)
#     )
#     if (!is.logical(val)) {
#       lgr::lgr$info( paste0("ERROR: Cache directory '",cache_dir, "' does not exist"))
#       return(-1)
#     }
#
#   }
#
#
#   ## if remote is TRUE
#   ## load from Zenodo (tagged)
#
#   read_dest <- NULL
#   write_dest <- NULL
#   oe_version <- paste0("v", utils::packageVersion("oncoEnrichR"))
#
#   write_to_cache <- T
#   zenodo_record_files <- data.frame()
#
#   if (remote == T) {
#
#     lgr::lgr$info(
#       paste0("Loading oncoEnrichR annotation datasets for version '",
#              oe_version,"' from Zenodo (https://zenodo.org/api/files/)"))
#     lgr::lgr$info("This may take several minutes depending on your connection/bandwidth")
#
#     zenodo <- zen4R::ZenodoManager$new()
#     zenodo_doi <- unique(oncoEnrichR::db_props$zenodo_doi)
#     oedb_rec <- zenodo$getRecordByDOI(zenodo_doi)
#     zenodo_record_files <- oedb_rec$listFiles(pretty = TRUE)
#
#     if (!is.null(cache_dir)) {
#
#       cache_version_dir <- file.path(cache_dir, oe_version)
#       write_dest <- cache_version_dir
#
#       if (!dir.exists(cache_version_dir)) {
#         system(paste0('mkdir ', cache_version_dir),
#                intern = F)
#         lgr::lgr$info( paste0("Data will be cached in ", cache_dir))
#
#       } else {
#         lgr::lgr$warn(paste0("An existing cache for '",oe_version, "' was found in ",
#                    cache_dir, " - will be overwritten (set 'remote' = FALSE to read from cache only)"))
#
#       }
#     } else {
#       write_to_cache <- F
#     }
#
#   } else {
#
#     write_to_cache <- F
#
#     if (!is.null(cache_dir)) {
#
#       cache_version_dir <- file.path(cache_dir, oe_version)
#
#       if (!dir.exists(cache_version_dir)) {
#         lgr::lgr$info( paste0("ERROR: No cache for 'v", oe_version, "' was found in ",cache_dir))
#         lgr::lgr$info( paste0("Set 'remote' = TRUE to reload data from zenodo.org and write to ", cache_dir))
#         return(-1)
#       } else {
#         read_dest <- cache_version_dir
#         lgr::lgr$info( paste0("Loading oncoEnrichR annotation datasets locally (cached) from: ",read_dest))
#       }
#
#     }
#   }
#
#   oedb <- list()
#   cache_in_use_log_printed <- 0
#
#   for (db in c("cancerdrugdb",
#               "genedb",
#               "hpa",
#               "ligandreceptordb",
#               "otdb",
#               "pfamdb",
#               "pathwaydb",
#               "depmapdb",
#               "survivaldb",
#               "release_notes",
#               "subcelldb",
#               "slparalogdb",
#               "tcgadb",
#               "tissuecelldb",
#               "tftargetdb")) {
#
#
#
#     if (remote == T & NROW(zenodo_record_files) > 0) {
#       zenodo_download_entry <-
#         zenodo_record_files[zenodo_record_files$filename == paste0(db,".rds"),]$download
#       read_dest <- stringr::str_replace(
#         zenodo_download_entry,
#         paste0(.Platform$file.sep , db, ".rds"),"")
#     }
#
#     db_dest <- file.path(read_dest, paste0(db,".rds"))
#     options(timeout=9999999)
#     if (remote == T) {
#       if (RCurl::url.exists(db_dest)) {
#         oedb[[db]] <- readRDS(url(db_dest,"rb"))
#       } else {
#         lgr::lgr$error(paste0("Could not retrieve data from ", db_dest))
#         return(-1)
#       }
#       checksum_db <- R.cache::getChecksum(oedb[[db]])
#       if (db == 'subcelldb') {
#         if ('COMPARTMENTSdb' %in% names(oedb[[db]])) {
#           checksum_db <- R.cache::getChecksum(oedb[[db]][['COMPARTMENTSdb']])
#         }
#       }
#       if (checksum_db ==
#          oncoEnrichR::db_props[oncoEnrichR::db_props$name == db,"checksum"]) {
#         lgr::lgr$info( paste0("'",
#                               db, "' - ",
#                               oe_version, " - ",
#                               checksum_db,
#                               " - verifies correctly"))
#       } else {
#         lgr::lgr$error(paste0("'",
#                               db, "' - ",
#                               oe_version, " - ",
#                               checksum_db,
#                               " - does not verify correctly"))
#       }
#       if (write_to_cache == T) {
#         cache_db_dest = file.path(write_dest, paste0(db,".rds"))
#         saveRDS(oedb[[db]], file = cache_db_dest)
#       }
#     } else {
#       if (file.exists(db_dest)) {
#         oedb[[db]] <- readRDS(file = db_dest)
#       } else {
#         lgr::lgr$error(paste0("Could not retrieve data from ", db_dest))
#         return()
#       }
#
#
#       checksum_db <- R.cache::getChecksum(oedb[[db]])
#
#       ##subcelldb's checksum does not recover correctly with all
#       ##list entries involved, choose COMPARTMENTS only
#       if (db == 'subcelldb') {
#         if ('COMPARTMENTSdb' %in% names(oedb[[db]])) {
#           checksum_db <- R.cache::getChecksum(oedb[[db]][['COMPARTMENTSdb']])
#         }
#       }
#
#
#       if (checksum_db ==
#          oncoEnrichR::db_props[oncoEnrichR::db_props$name == db,"checksum"]) {
#         lgr::lgr$info( paste0("'",
#                               db, "' - ",
#                               oe_version, " - ",
#                               checksum_db,
#                               " - verifies correctly"))
#       } else {
#         lgr::lgr$error(paste0("'",
#                               db, "' - ",
#                               oe_version, " - ",
#                               checksum_db,
#                               " - does not verify correctly"))
#       }
#
#     }
#   }
#
#
#   return(oedb)
#
# }


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
#' @param ppi_string_min_score minimum score (between 0 and 1) for confidence of retrieved protein-protein interactions (STRING)
#' @param ppi_string_network_type STRING network type ('physical' or 'functional')
#' @param ppi_biogrid_min_evidence minimum number of evidence items for protein-protein interactions shown (BIOGRID)
#' @param ppi_add_nodes number of nodes to add to target set when computing the protein-protein interaction network (STRING/BIOGRID)
#' @param ppi_node_shadow show shadow for nodes in the displayed PPI network (STRING/BIOGRID)
#' @param ppi_show_drugs logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network
#' @param ppi_show_isolated_nodes logical indicating if targets/nodes without any interactions should be displayed in the protein-protein interaction network
#' @param bgset_description character indicating type of background (e.g. "All lipid-binding proteins (n = 200)")
#' @param bgset_id_type character indicating source of query (one of "uniprot_acc", "symbol",
#' "entrezgene", or "ensembl_gene","ensembl_mrna","refseq_mrna","ensembl_protein","refseq_protein")
#' @param enrichment_p_value_cutoff cutoff p-value for enrichment/over-representation analysis
#' @param enrichment_p_value_adj one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param enrichment_q_value_cutoff cutoff q-value for enrichment analysis
#' @param enrichment_min_geneset_size minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param enrichment_max_geneset_size maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
#' @param enrichment_plot_num_terms number of top enriched Gene Ontology terms (max) to show in enrichment barplot
#' @param enrichment_simplify_go remove highly similar GO terms in results from GO enrichment/over-representation analysis
#' @param subcellcomp_min_confidence minimun confidence level for subcellular compartment annotation in COMPARTMENTS (min = 3, max = 5)
#' @param subcellcomp_min_channels minimum number of channels that support a subcellular annotation in COMPARTMENTS
#' @param subcellcomp_show_cytosol logical indicating if subcellular heatmap should highlight cytosol as a subcellular protein location or not
#' @param regulatory_min_confidence minimum confidence level for regulatory interactions (TF-target) retrieved from DoRothEA ('A','B','C', or 'D')
#' @param fitness_max_score maximum loss-of-fitness score (scaled Bayes factor from BAGEL) for genes retrieved from Project Score
#' @param show_ppi logical indicating if report should contain protein-protein interaction data (STRING)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform, association_score >= 0.05, support from at least two data types)
#' @param show_top_diseases_only logical indicating if report should contain top (n = 20) disease associations only pr. query gene (Open Targets Platform)
#' @param show_cancer_hallmarks logical indicating if report should contain annotations/evidence of cancer hallmarks per query gene (COSMIC/Open Targets Platform)
#' @param show_drug logical indicating if report should contain targeted cancer drug information
#' @param show_enrichment logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME, NetPath, WikiPathways)
#' @param show_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of query set with oncogenes/tumor suppressor genes
#' @param show_cell_tissue logical indicating if report should contain tissue-specificity and single cell-type specificity assessments (Human Protein Atlas)
#' of target genes, using data from the Human Protein Atlas
#' @param show_ligand_receptor logical indicating if report should contain ligand-receptor interactions (CellChatDB)
#' @param show_regulatory logical indicating if report should contain data on transcription factor (TF) - target interactions relevant for the query set (DoRothEA)
#' @param show_unknown_function logical indicating if report should highlight target genes with unknown or poorly defined functions (GO/Uniprot KB/NCBI)
#' @param show_prognostic logical indicating if mRNA-based (single-gene) prognostic associations to cancer types should be listed (Human Protein Atlas/TCGA)
#' @param show_subcell_comp logical indicating if report should provide subcellular compartment annotations (COMPARTMENTS)
#' @param show_synleth logical indicating if report should list overlap with predicted synthetic lethality interactions (gene paralogs only, De Kegel et al., Cell Systems, 2021)
#' @param show_fitness logical indicating if report should provide fitness scores and target priority scores from CRISPR/Cas9 loss-of-fitness screens (Project Score)
#' @param show_complex logical indicating if report should provide target memberships in known protein complexes (ComplexPortal/Compleat/PDB/CORUM)
#' @param show_domain logical indicating if report should provide target memberships in known protein domains (Pfam)
#'
#' @keywords internal
#'

init_report <- function(oeDB,
                        project_title = "_Project title_",
                        project_owner = "_Project owner_",
                        html_floating_toc = T,
                        html_report_theme = "default",
                        query_id_type = "symbol",
                        ignore_id_err = TRUE,
                        project_description = "_Project description_",
                        ppi_string_min_score = 0.9,
                        ppi_string_network_type = "functional",
                        ppi_biogrid_min_evidence = 3,
                        ppi_add_nodes = 30,
                        ppi_show_isolated_nodes = FALSE,
                        ppi_node_shadow = TRUE,
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
                        subcellcomp_min_confidence = 3,
                        subcellcomp_min_channels = 1,
                        subcellcomp_show_cytosol = F,
                        regulatory_min_confidence = "D",
                        fitness_max_score = -2,
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
                        show_synleth = T,
                        show_complex = T,
                        show_domain = T) {

  ## report object
  rep <- list()

  stopifnot(!is.null(oeDB))
  stopifnot(!is.null(oeDB$release_notes))
  stopifnot(!is.null(oeDB$tcgadb))
  stopifnot(!is.null(oeDB$subcelldb$gganatogram_legend))

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
  rep[["config"]][["show"]][["protein_domain"]] <- show_domain
  rep[["config"]][["show"]][["aberration"]] <- show_aberration
  rep[["config"]][["show"]][["coexpression"]] <- show_coexpression
  rep[["config"]][["show"]][["subcellcomp"]] <- show_subcell_comp
  rep[["config"]][["show"]][["fitness"]] <- show_fitness
  rep[["config"]][["show"]][["cell_tissue"]] <- show_cell_tissue
  rep[["config"]][["show"]][["regulatory"]] <- show_regulatory
  rep[["config"]][["show"]][["ligand_receptor"]] <- show_ligand_receptor
  rep[["config"]][["show"]][["cancer_hallmark"]] <- show_cancer_hallmarks
  rep[["config"]][["show"]][["unknown_function"]] <- show_unknown_function
  rep[["config"]][["show"]][["synleth"]] <- show_synleth
  rep[["config"]][["show"]][["cancer_prognosis"]] <-
    show_prognostic

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
    regulatory_min_confidence

  ## config - color codes and thresholds for
  ## synthetic lethality prediction percentiles
  rep[["config"]][["synleth"]] <- list()
  rep[["config"]][["synleth"]][["breaks"]] <-
    c(2, 5, 10)
  rep[["config"]][["synleth"]][["colors"]] <-
    c("#08306b","#08519c","#2171b5","#b8b8ba")

  ## config - disease - color codes and
  ## thresholds for quantitative target-disease associations
  rep[["config"]][["disease"]] <- list()
  rep[["config"]][["disease"]][["breaks"]] <-
     c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  rep[["config"]][["disease"]][["colors"]] <-
    c("#b8b8ba",
      "#deebf7",
      "#c6dbef",
      "#9ecae1",
      "#6baed6",
      "#4292c6",
      "#2171b5",
      "#08519c",
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

  ## config - fitness/loss_of_function - plot height for hits in CRISPR screens
  rep[["config"]][["fitness"]] <- list()
  rep[["config"]][["fitness"]][["plot_height_fitness"]] <- 10
  rep[['config']][["fitness"]][["max_BF_score"]] <- fitness_max_score

  ## config/ppi - protein-protein interaction settings
  rep[["config"]][["ppi"]] <- list()
  rep[["config"]][["ppi"]][["string"]] <- list()
  rep[["config"]][["ppi"]][["string"]][["minimum_score"]] <- ppi_string_min_score
  rep[["config"]][["ppi"]][["string"]][["visnetwork_shape"]] <- "dot"
  rep[["config"]][["ppi"]][["string"]][["visnetwork_shadow"]] <- ppi_node_shadow
  rep[["config"]][["ppi"]][["string"]][["show_drugs"]] <- ppi_show_drugs
  rep[["config"]][["ppi"]][["string"]][["add_nodes"]] <- ppi_add_nodes
  rep[["config"]][["ppi"]][["string"]][["query_type"]] <- "network"
  rep[["config"]][["ppi"]][["string"]][["network_type"]] <- ppi_string_network_type
  rep[["config"]][["ppi"]][["string"]][["show_isolated_nodes"]] <- ppi_show_isolated_nodes


  rep[["config"]][["ppi"]][["biogrid"]] <- list()
  rep[["config"]][["ppi"]][["biogrid"]][['minimum_evidence']] <- ppi_biogrid_min_evidence
  rep[["config"]][["ppi"]][["biogrid"]][["visnetwork_shape"]] <- "dot"
  rep[["config"]][["ppi"]][["biogrid"]][["visnetwork_shadow"]] <- ppi_node_shadow
  rep[["config"]][["ppi"]][["biogrid"]][["show_drugs"]] <- ppi_show_drugs
  rep[["config"]][["ppi"]][["biogrid"]][["add_nodes"]] <- ppi_add_nodes
  rep[["config"]][["ppi"]][["biogrid"]][["show_isolated_nodes"]] <- ppi_show_isolated_nodes

  ## config/enrichment - general functional enrichment settings
  rep[["config"]][["enrichment"]] <- list()
  rep[["config"]][["enrichment"]][["p_value_cutoff"]] <-
    enrichment_p_value_cutoff
  rep[["config"]][["enrichment"]][["q_value_cutoff"]] <-
    enrichment_q_value_cutoff
  rep[["config"]][["enrichment"]][["p_adjust_method"]] <-
    enrichment_p_value_adj
  rep[["config"]][["enrichment"]][["min_gs_size"]] <-
    enrichment_min_geneset_size
  rep[["config"]][["enrichment"]][["max_gs_size"]] <-
    enrichment_max_geneset_size
  rep[["config"]][["enrichment"]][["simplify_go"]] <-
    enrichment_simplify_go
  rep[["config"]][["enrichment"]][["bgset_description"]] <-
    bgset_description
  rep[["config"]][["enrichment"]][["bgset_size"]] <- 0
  rep[["config"]][["enrichment"]][["enrichment_plot_num_terms"]] <-
    enrichment_plot_num_terms


  ## specifiy plot height (tile/heatmap) - genes along y-axis
  ## tumor types / tissues / celltypes along x-axis
  rep[["config"]][["aberration"]] <- list()
  rep[["config"]][["aberration"]][["plot_height"]] <- 14

  ## config/prognosis - settings for prognostic associations
  rep[["config"]][["prognosis"]] <- list()
  rep[["config"]][["prognosis"]][["breaks"]] <-
    c(3, 5.2, 7.4, 9.5, 11.7, 13.9)
  rep[["config"]][["prognosis"]][["colors_unfavorable"]] <-
    c("#fee5d9",
      "#fcbba1",
      "#fc9272",
      "#fb6a4a",
      "#ef3b2c",
      "#cb181d",
      "#99000d")
  rep[["config"]][["prognosis"]][["colors_favorable"]] <-
    c("#edf8e9",
      "#c7e9c0",
      "#a1d99b",
      "#74c476",
      "#41ab5d",
      "#238b45",
      "#005a32")

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
    c("#084594",
      "#2171B5",
      "#4292C6",
      "#6BAED6",
      "#9ECAE1",
      "#b8b8ba")

  rep[["config"]][["unknown_function"]] <- list()
  rep[["config"]][["unknown_function"]][["rank"]] <-
    c(1, 2, 3, 4, 5, 6)
  rep[["config"]][["unknown_function"]][["colors"]] <-
    c("#99000d",
      "#cb181d",
      "#ef3b2c",
      "#fb6a4a",
      "#fc9272",
      "#fcbba1")

  rep[["config"]][["unknown_function"]][['num_candidates']] <-
    oeDB[['genedb']]$all |>
    dplyr::filter(.data$gene_biotype == "protein-coding") |>
    dplyr::filter(!is.na(.data$unknown_function_rank)) |>
    dplyr::filter(.data$unknown_function_rank <= 6) |>
    nrow()

  rep[['config']][['subcellcomp']] <- list()
  rep[['config']][['subcellcomp']][['minimum_confidence']] <-
    subcellcomp_min_confidence
  rep[['config']][['subcellcomp']][['minimum_channels']] <-
    subcellcomp_min_channels
  rep[['config']][['subcellcomp']][['show_cytosol']] <-
    subcellcomp_show_cytosol
  rep[['config']][['subcellcomp']][['gganatogram_legend']] <-
    oeDB$subcelldb$gganatogram_legend

  rep[['config']][['complex']] <- list()
  rep[['config']][['complex']][['breaks']] <-
    rep[['config']][['disease']][['breaks']]
  rep[['config']][['complex']][['colors']] <-
    rep[['config']][['disease']][['colors']]


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
                     "protein_domain",
                     "ligand_receptor",
                     "subcellcomp",
                     "fitness",
                     "synleth",
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

  ## synthetic lethality
  rep[["data"]][["synleth"]][['both_in_pair']] <- data.frame()
  rep[["data"]][["synleth"]][['single_pair_member']] <- data.frame()

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

  ## protein domains
  rep[["data"]][["cancer_hallmark"]][["protein_domain"]] <- data.frame()

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
  for (db in c('string','biogrid')) {
    rep[["data"]][["ppi"]][[db]] <- list()
    rep[["data"]][["ppi"]][[db]][["source"]] <- toupper(db)
    rep[["data"]][["ppi"]][[db]][["complete_network"]] <- NULL
    rep[["data"]][["ppi"]][[db]][["hubscores"]] <- data.frame()
    rep[["data"]][["ppi"]][[db]][["community_network"]] <- NULL
  }

  ## functional enrichment
  for (c in c("go","msigdb","wikipathways","kegg","netpath")) {
    rep[["data"]][["enrichment"]][[c]] <- data.frame()
  }

  ## unknown function
  rep[["data"]][["unknown_function"]][["hits_df"]] <- data.frame()

  ## Fitness scores and target priority scores from CRISPR/Cas9 screens
  ## (Project Score (ps))
  rep[["data"]][["fitness"]][['fitness_scores']] <- list()
  rep[["data"]][["fitness"]][['fitness_scores']][["targets"]] <- data.frame()
  rep[["data"]][["fitness"]][['fitness_scores']][["n_targets"]] <- 0
  rep[["data"]][["fitness"]][['target_priority_scores']] <- list()
  rep[["data"]][["fitness"]][['target_priority_scores']][['targets']] <-
    data.frame()
  rep[["data"]][["fitness"]][['target_priority_scores']][['n_pri_targets']] <-
    0


  ## TCGA co-expression
  rep[["data"]][["tcga"]][["coexpression"]] <- data.frame()

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
    if (v != "snv_indel") {
      rep[["data"]][["tcga"]][["aberration"]][["matrix"]][[v]] <- NULL
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- data.frame()
    }
    else{
      rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <- list()
      i <- 1
      while(i <= nrow(oeDB$tcgadb$maf_codes)) {
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
#' Function that interrogates a list of human gene identifiers for
#' cancer relevance. Multiple perspectives are offered, including
#' tumor aberration and co-expression patterns, druggability,
#' protein-protein interactions, gene fitness effects, regulatory
#' interactions, subcellular compartment enrichment, pathway enrichment,
#' synthetic lethality interactions, prognostic associations, and more.
#'
#' @param query character vector with gene/query identifiers
#' @param oeDB oncoEnrichR data repository object - as returned from `load_db()`
#' @param query_id_type character indicating source of query (one of
#' "uniprot_acc", "symbol","entrezgene", or "ensembl_gene", "ensembl_mrna",
#' "refseq_mrna", "ensembl_protein", "refseq_protein")
#' @param html_floating_toc logical - float the table of contents to the left
#' of the main document content (HTML report). The floating table of contents
#' will always be visible even when the document is scrolled
#' @param html_report_theme Bootswatch theme for HTML report (any of
#' "bootstrap", "cerulean", "cosmo", "default",
#' "flatly", "journal", "lumen", "paper", "sandstone", "simplex",
#' "spacelab", "united", "yeti")
#' @param ignore_id_err logical indicating if analysis should
#' continue when uknown query identifiers are encountered
#' @param project_title project title (title of report)
#' @param project_owner name of project owner
#' @param project_description project background information
#' @param bgset character vector with gene identifiers, used as
#' reference/background for enrichment/over-representation analysis
#' @param bgset_id_type character indicating source of background
#' ("uniprot_acc", "symbol", "entrezgene",
#' "ensembl_gene", "ensembl_mrna", "refseq_mrna", "ensembl_protein",
#' "refseq_protein"), default: "symbol"
#' @param bgset_description character indicating type of background
#' (e.g. "All lipid-binding proteins (n = 200)")
#' @param enrichment_p_value_cutoff cutoff p-value for enrichment/
#' over-representation analysis (default: 0.05)
#' @param enrichment_p_value_adj one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none" (clusterProfiler, default: "BH")
#' @param enrichment_q_value_cutoff cutoff q-value for enrichment
#' analysis (clusterProfiler, default: 0.2)
#' @param enrichment_min_geneset_size minimal size of geneset annotated by term
#' for testing in enrichment/over-representation analysis (clusterProfiler,
#' default: 10)
#' @param enrichment_max_geneset_size maximal size of geneset annotated by term
#' for testing in enrichment/over-representation analysis (clusterProfiler,
#' default: 500)
#' @param enrichment_plot_num_terms number of top enriched Gene Ontology terms
#' (max) to show in enrichment barplot (default: 15)
#' @param enrichment_simplify_go remove highly similar GO terms in results from
#' GO enrichment/over-representation analysis (default: TRUE)
#' @param subcellcomp_min_confidence minimum confidence level for subcellular
#' compartment annotation in COMPARTMENTS (min = 3, max = 5, default: 3)
#' @param subcellcomp_min_channels minimum number of channels that support a
#' subcellular compartment annotation in COMPARTMENTS (min = 1,
#' max = 3, default: 1)
#' @param subcellcomp_show_cytosol logical indicating if subcellular heatmap
#' should show highlight proteins located in the cytosol or not (default: FALSE)
#' @param regulatory_min_confidence minimum confidence level for regulatory
#' interactions (TF-target) retrieved from DoRothEA ('A','B','C', or 'D',
#' default: 'D')
#' @param fitness_max_score maximum loss-of-fitness score (scaled Bayes factor
#' from BAGEL) for genes retrieved from DepMap/Project Score, default:-2
#' @param ppi_add_nodes number of nodes to add to target set when computing the
#' protein-protein interaction network (STRING/BioGRID, default: 30)
#' @param ppi_string_min_score minimum score (between 0 and 1) for confidence of
#' retrieved protein-protein interactions (STRING, default: 0.9)
#' @param ppi_string_network_type type of network to show for interactions in
#' STRING ('functional' or 'physical', default: 'functional')
#' @param ppi_biogrid_min_evidence minimum number of evidence items required
#' for protein-protein interactions retrieved (BioGRID, default: 3)
#' @param ppi_node_shadow show shadow for nodes in the displayed PPI
#' network (default: TRUE)
#' @param ppi_show_drugs logical indicating if targeted drugs (>= phase 3)
#' should be displayed in protein-protein interaction networks (default: TRUE)
#' @param ppi_show_isolated_nodes logical indicating if targets/nodes without
#' any interactions should be displayed in the protein-protein
#' interaction networks (default: FALSE)
#' @param show_ppi logical indicating if report should contain protein-protein
#' interaction views of the query set (STRING and BioGRID, default: TRUE)
#' @param show_disease logical indicating if report should contain disease
#' associations (Open Targets Platform, association_score >= 0.05, support
#' from at least two data types), and tumor suppressor/oncogene annotations (
#' default: TRUE)
#' @param show_top_diseases_only logical indicating if report should contain
#' top (n = 20) disease associations only pr. query gene
#' default: TRUE)
#' @param show_cancer_hallmarks logical indicating if report should
#' contain annotations/evidence of cancer hallmarks per query gene
#' (COSMIC/Open Targets Platform, default: TRUE
#' @param show_drug logical indicating if report should contain targeted
#' cancer drug information (default: TRUE)
#' @param show_enrichment logical indicating if report should contain functional
#' enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME,
#' NetPath, WikiPathways, default: TRUE)
#' @param show_aberration logical indicating if report should contain TCGA
#' aberration plots (amplifications/deletions, default: TRUE)
#' @param show_coexpression logical indicating if report should contain TCGA
#' co-expression data (RNAseq) of query set with oncogenes/tumor
#' suppressor genes (default: TRUE)
#' @param show_cell_tissue logical indicating if report should contain
#' tissue-specificity and single cell-type specificity assessments
#' (Human Protein Atlas) of target genes (default: FALSE)
#' @param show_ligand_receptor logical indicating if report should contain
#' ligand-receptor interactions (CellChatDB, default: TRUE)
#' @param show_regulatory logical indicating if report should contain data on
#' transcription factor (TF) - target interactions relevant for the
#' query set (DoRothEA, default: TRUE)
#' @param show_unknown_function logical indicating if report should highlight
#' target genes with unknown or poorly defined functions
#' (GO/Uniprot KB/NCBI, default: TRUE)
#' @param show_prognostic logical indicating if mRNA-based (single-gene)
#' prognostic associations to cancer types should be listed
#' (Human Protein Atlas/TCGA, default: TRUE
#' @param show_subcell_comp logical indicating if report should provide
#' subcellular compartment annotations (COMPARTMENTS, default: TRUE)
#' @param show_synleth logical indicating if report should list overlap with
#' predicted synthetic lethality interactions (gene paralogs only,
#' De Kegel et al., Cell Systems, 2021). Default: TRUE
#' @param show_fitness logical indicating if report should provide fitness
#' scores and target priority scores from CRISPR/Cas9 loss-of-fitness
#' screens (DepMap/Project Score, default: TRUE)
#' @param show_complex logical indicating if report should provide target
#' memberships in known protein complexes
#' (ComplexPortal/Compleat/hu.MAP2/PDB/CORUM, default: TRUE)
#' @param show_domain logical indicating if report should provide target
#' memberships in known protein domains (Pfam, default: TRUE)

#' @param ... arguments for Galaxy/web-based processing
#'
#' @return
#' An oncoEnrichR report list object, with two main elements,
#' `data` and `config`. `data` contains data that goes into
#' each section of the output reports, `config` contains all
#' metadata for annotation resources used.
#'
#' @export
#'
onco_enrich <- function(query = NULL,
                        oeDB = NULL,
                        query_id_type = "symbol",
                        ignore_id_err = TRUE,
                        html_floating_toc = T,
                        html_report_theme = "default",
                        project_title = "_Project title_",
                        project_owner = "_Project owner_",
                        project_description = "_Project description_",
                        bgset = NULL,
                        bgset_id_type = "symbol",
                        bgset_description = "All protein-coding genes",
                        enrichment_p_value_cutoff = 0.05,
                        enrichment_p_value_adj = "BH",
                        enrichment_q_value_cutoff = 0.2,
                        enrichment_min_geneset_size = 10,
                        enrichment_max_geneset_size = 500,
                        enrichment_plot_num_terms = 20,
                        enrichment_simplify_go = TRUE,
                        subcellcomp_min_confidence = 3,
                        subcellcomp_min_channels = 1,
                        subcellcomp_show_cytosol = FALSE,
                        regulatory_min_confidence = "D",
                        fitness_max_score = -2,
                        ppi_add_nodes = 30,
                        ppi_string_min_score = 0.9,
                        ppi_string_network_type = "functional",
                        ppi_biogrid_min_evidence = 3,
                        ppi_node_shadow = TRUE,
                        ppi_show_drugs = TRUE,
                        ppi_show_isolated_nodes = FALSE,
                        show_ppi = TRUE,
                        show_disease = TRUE,
                        show_top_diseases_only = TRUE,
                        show_cancer_hallmarks = TRUE,
                        show_drug = TRUE,
                        show_enrichment = TRUE,
                        show_aberration = TRUE,
                        show_coexpression = TRUE,
                        show_cell_tissue = FALSE,
                        show_ligand_receptor = TRUE,
                        show_regulatory = TRUE,
                        show_unknown_function = TRUE,
                        show_prognostic = TRUE,
                        show_subcell_comp = TRUE,
                        show_synleth = TRUE,
                        show_fitness = TRUE,
                        show_complex = TRUE,
                        show_domain = TRUE,
                        ...) {


  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  dot_args <- list(...)

  if (is.null(oeDB)) {
    lgr::lgr$info( paste0(
      "ERROR: mandatory argument 'oeDB' cannot be NULL"))
    return()
  }
  oedb_val <- validate_db(oeDB)
  if (oedb_val != 0) {
    return()
  }

  if (is.null(query)) {
    lgr::lgr$info( paste0(
      "ERROR: mandatory argument 'query' cannot be NULL"))
    return()
  }
  if (!is.character(query)) {
    lgr::lgr$info( paste0(
      "ERROR: mandatory argument 'query' is of wrong type (not character)"))
    return()
  }

  if (!is.null(bgset)) {
    val <- assertthat::validate_that(
      is.character(bgset))
    if (!is.logical(val)) {
      lgr::lgr$info( paste0(
        "ERROR: ", val))
      return()
    }
  }

  val <- assertthat::validate_that(length(query) >= 2)
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: query set must contain at least two entries - length of query is: ", length(query)))
    return()
  }

  val <- assertthat::validate_that(
    enrichment_p_value_adj %in% c("holm", "hochberg",
                                     "hommel", "bonferroni",
                                     "BH", "BY",
                                     "fdr", "none")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'enrichment_p_value_adj' must take on any of the following values: ",
      "'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'",
      " (value provided was '", enrichment_p_value_adj,"')"))
    return()
  }

  val <- assertthat::validate_that(
    regulatory_min_confidence %in% c("A", "B",
                                          "C", "D")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'regulatory_min_confidence' must take on any of the following values: 'A', 'B', 'C', 'D'",
      " (value provided was '", regulatory_min_confidence,"')"))
    return()
  }

  val <- assertthat::validate_that(
    html_report_theme %in% c("bootstrap","cerulean","cosmo","default",
                        "flatly","journal","lumen","paper","sandstone",
                        "simplex","spacelab","united","yeti")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'html_report_theme' must take on any of the following values: ",
      "'bootstrap', 'cerulean', 'cosmo', 'default', 'flatly', 'journal', 'lumen',",
      "'paper', 'sandstone', 'simplex', 'spacelab', 'united', 'yeti'",
      " (value provided was '", enrichment_p_value_adj,"')"))
    return()
  }

  ## Number of allowed query genes
  oncoenrichr_query_limit <- 1000

  if (length(names(dot_args)) > 0) {
    if ("galaxy" %in% names(dot_args))
      lgr::lgr$info(
                 "NOTE: Running oncoEnrichR workflow in Galaxy mode")
      if (is.logical(dot_args$galaxy)) {
        if (dot_args$galaxy == T) {
          oncoenrichr_query_limit <- 1000
        }
      }
  }

  if (length(query) > oncoenrichr_query_limit) {
    lgr::lgr$info(
      paste0("WARNING: oncoEnrichR is limited to the analysis of n = ",
             oncoenrichr_query_limit," entries. Query contained n = ",
             length(query), " identifiers, limiting to top ",
             oncoenrichr_query_limit," genes"))
    query <- utils::head(unique(query),
                         oncoenrichr_query_limit)
  }

  val <- assertthat::validate_that(
    query_id_type %in%
      c("symbol", "entrezgene",
        "refseq_mrna", "ensembl_mrna",
        "refseq_protein", "ensembl_protein",
        "uniprot_acc",
        "ensembl_gene")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'query_id_type' must take on of the following values: ",
      "'symbol', 'entrezgene', 'refseq_mrna', 'ensembl_mrna', ",
      "'refseq_protein', 'ensembl_protein', 'uniprot_acc', 'ensembl_gene'",
      " (value provided was '", query_id_type,"')"))
    return()
  }

  val <- assertthat::validate_that(
    bgset_id_type %in%
      c("symbol", "entrezgene", "refseq_mrna", "ensembl_mrna",
        "refseq_protein", "ensembl_protein", "uniprot_acc",
        "ensembl_gene")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'bgset_id_type' must take on any of the following values: ",
      "'symbol', 'entrezgene', 'refseq_mrna', 'ensembl_mrna', ",
      "'refseq_protein', 'ensembl_protein', 'uniprot_acc', 'ensembl_gene'",
      " (value provided was '", bgset_id_type,"')"))
    return()
  }

  val <- fitness_max_score <= 0 & is.numeric(fitness_max_score)

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'fitness_max_score' must be a value (scaled Bayes factor from BAGEL) less than zero ",
      "(current type and value: '",typeof(fitness_max_score),"' - ",
      fitness_max_score,")")
    )
    return()
  }


  val <-
    (is.numeric(ppi_string_min_score)) & ## check that number is numeric
      (ppi_string_min_score > 0) &
      (ppi_string_min_score <= 1)

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'ppi_string_min_score' must take a numeric value - greater than 0 and less than 1 ",
      "(current type and value: '",typeof(ppi_string_min_score),"' - ",
      ppi_string_min_score,")")
    )
    return()
  }

  val <-
    (ppi_biogrid_min_evidence %% 1 == 0) & ## check that number is whole integer
    (ppi_biogrid_min_evidence >= 2) &
    (ppi_biogrid_min_evidence <= 10)

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'ppi_biogrid_min_evidence' must be an integer/whole number and take a value from 2 to 10 ",
      "(current type and value: '",typeof(ppi_biogrid_min_evidence),"' - ",
      ppi_biogrid_min_evidence,")")
    )
    return()
  }

  val <-
    enrichment_plot_num_terms %% 1 == 0 &
      enrichment_plot_num_terms >= 10 &
      enrichment_plot_num_terms <= 30

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'enrichment_plot_num_terms' must be an integer/whole number and take a value from 10 to 30 ",
      "(current type and value: '",typeof(enrichment_plot_num_terms),"' - ",
      enrichment_plot_num_terms,")")
    )
    return()
  }

  val <-
    subcellcomp_min_confidence %% 1 == 0 &
    subcellcomp_min_confidence >= 3 &
    subcellcomp_min_confidence <= 5

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'subcellcomp_min_confidence' must be an integer/whole number and take a value from 3 to 5 ",
      "(current type and value: '",typeof(subcellcomp_min_confidence),"' - ",
      subcellcomp_min_confidence,")")
    )
    return()
  }

  # val <- T
  # subcellcomp_num_channels_selected <- 0
  # for (e in subcellcomp_channel_types) {
  #   if (!(e %in% c('Experimental', 'Text mining', 'Knowledge'))) {
  #     val <- F
  #   } else {
  #     subcellcomp_num_channels_selected <-
  #       subcellcomp_num_channels_selected + 1
  #   }
  # }
  #
  # if (val == F) {
  #   lgr::lgr$info( paste0(
  #     "ERROR: 'subcellcomp_channel_types' must be any selection of 'Knowledge' ",
  #     ", 'Experimental', and 'Text mining', (current type and value: '",
  #     typeof(subcellcomp_channel_types),"' - ",
  #     paste(subcellcomp_channel_types, collapse = ', '),")"
  #   ))
  #   return()
  # }

  val <-
    subcellcomp_min_channels %% 1 == 0 &
    subcellcomp_min_channels > 0 &
    subcellcomp_min_channels <= 3

  if (val == F) {
    lgr::lgr$info( paste0(
      "ERROR: 'subcellcomp_min_channels' must be an integer/whole number and not larger than the selected number of channel types",
      "(current type and value: '",typeof(subcellcomp_min_channels),"' - ",
      subcellcomp_min_channels,")")
    )
    return()
  }


  val <- assertthat::validate_that(
    is.numeric(enrichment_p_value_cutoff) &
      enrichment_p_value_cutoff > 0 & enrichment_p_value_cutoff < 1)

  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'enrichment_p_value_cutoff' must be of type numeric and be greater than 0 and less than 1 ",
      "(current type and value: '",typeof(enrichment_p_value_cutoff),"' - ",
      enrichment_p_value_cutoff,")")
    )
    return()
  }

  val <- assertthat::validate_that(
    is.numeric(enrichment_q_value_cutoff) &
      enrichment_q_value_cutoff > 0 & enrichment_q_value_cutoff < 1)

  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: 'enrichment_q_value_cutoff' must be of type numeric and be greater than 0 and less than 1 ",
      "(current type and value: '",typeof(enrichment_q_value_cutoff),"' - ",
      enrichment_q_value_cutoff,")")

    )
    return()
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
    ppi_string_min_score = ppi_string_min_score,
    ppi_string_network_type = ppi_string_network_type,
    ppi_biogrid_min_evidence = ppi_biogrid_min_evidence,
    ppi_add_nodes = ppi_add_nodes,
    ppi_show_isolated_nodes = ppi_show_isolated_nodes,
    ppi_node_shadow = ppi_node_shadow,
    bgset_description =
      bgset_description,
    bgset_id_type = bgset_id_type,
    enrichment_p_value_cutoff = enrichment_p_value_cutoff,
    enrichment_p_value_adj = enrichment_p_value_adj,
    enrichment_q_value_cutoff = enrichment_q_value_cutoff,
    enrichment_min_geneset_size = enrichment_min_geneset_size,
    enrichment_max_geneset_size = enrichment_max_geneset_size,
    enrichment_simplify_go = enrichment_simplify_go,
    enrichment_plot_num_terms = enrichment_plot_num_terms,
    subcellcomp_show_cytosol = subcellcomp_show_cytosol,
    subcellcomp_min_confidence = subcellcomp_min_confidence,
    subcellcomp_min_channels = subcellcomp_min_channels,
    regulatory_min_confidence = regulatory_min_confidence,
    fitness_max_score = fitness_max_score,
    show_ppi = show_ppi,
    ppi_show_drugs = ppi_show_drugs,
    show_disease = show_disease,
    show_top_diseases_only = show_top_diseases_only,
    show_cancer_hallmarks = show_cancer_hallmarks,
    show_drug = show_drug,
    show_enrichment = show_enrichment,
    show_aberration = show_aberration,
    show_coexpression = show_coexpression,
    show_subcell_comp = show_subcell_comp,
    show_fitness = show_fitness,
    show_cell_tissue = show_cell_tissue,
    show_ligand_receptor = show_ligand_receptor,
    show_regulatory = show_regulatory,
    show_unknown_function = show_unknown_function,
    show_synleth = show_synleth,
    show_prognostic =
      show_prognostic,
    show_domain = show_domain,
    show_complex = show_complex)

  ## validate query gene set
  qgenes_match <-
    validate_query_genes(
      qgenes = query,
      q_id_type = query_id_type,
      ignore_id_err = ignore_id_err,
      genedb = oeDB[['genedb']][['all']],
      transcript_xref = oeDB[['genedb']][['transcript_xref']])

  val <- assertthat::validate_that(NROW(qgenes_match$found) >= 2)
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: query set must contain at least two valid entries - number of validated entries: ", NROW(qgenes_match$found)))
    return()
  }

  ## assign validation result to report object
  onc_rep[['data']][['query']][['target']] <-
    qgenes_match[['all']]
  onc_rep[["data"]][["query"]][["validation_status"]] <-
    qgenes_match[["match_status"]]

  ## if query list contains non-validated entries and 'ignore_id_err' is
  ## set to FALSE, oncoEnrichR analysis will halt (no analysis modules
  ## will be included in report)
  if (qgenes_match[["match_status"]] == "imperfect_stop") {
    for (e in c("disease",
               "drug",
               "enrichment",
               "protein_complex",
               "protein_domain",
               "ppi",
               "aberration",
               "coexpression",
               "cell_tissue",
               "cancer_hallmark",
               "fitness",
               "regulatory_interactions",
               "ligand_receptor",
               "subcellcomp",
               "synleth",
               "unknown_function",
               "cancer_prognosis")) {
      onc_rep[['config']][['show']][[e]] <- F
    }
    return(onc_rep)
  }

  ## If list of genes are less than 5, pathway and GO enrichment
  ## is not performed (too few genes)
  if (NROW(qgenes_match[["found"]]) < 5) {
    onc_rep[['config']][['show']][['enrichment']] <- F
    lgr::lgr$info(
      paste0("WARNING: Function (GO) and pathway enrichment is NOT performed for gene sets of size < 5. Query contained n = ",
             NROW(qgenes_match[["found"]])," valid entries"))
  }


  ## validate background gene set (if provided)
  background_entrezgene <- NULL
  background_genes_match <- NULL
  if (!is.null(bgset)) {
    background_genes_match <-
      validate_query_genes(
        bgset,
        q_id_type = bgset_id_type,
        ignore_id_err = ignore_id_err,
        genedb = oeDB[['genedb']][['all']],
        transcript_xref = oeDB[['genedb']][['transcript_xref']],
        qtype = "background")
    if (background_genes_match[["match_status"]] == "imperfect_stop" |
        NROW(background_genes_match[['found']]) <= 1) {
      lgr::lgr$info( paste0("WARNING: Background geneset not defined properly - ",
                        "using all protein-coding genes instead"))
      background_entrezgene <- as.character(
        unique(oeDB[['genedb']][['all']]$entrezgene)
      )
    } else {
      background_entrezgene <- as.character(
        unique(background_genes_match[["found"]]$entrezgene)
      )
      if (onc_rep[['config']][['enrichment']][['bgset_description']] ==
         "All protein-coding genes") {
        lgr::lgr$info( paste0("WARNING: Description of background set is not set"))
        onc_rep[['config']][['enrichment']][['bgset_description']] <-
          "Undefined"

      }
    }
  } else {
    bg <- dplyr::select(oeDB[['genedb']][['all']],
                        c("entrezgene",
                        "gene_biotype")) |>
      dplyr::filter(!is.na(.data$entrezgene) &
                      .data$gene_biotype == "protein-coding") |>
      dplyr::distinct()
    background_entrezgene <- as.character(bg$entrezgene)
  }
  onc_rep[['config']][['enrichment']][['bgset_size']] <-
    length(background_entrezgene)


  query_entrezgene <- unique(qgenes_match[["found"]]$entrezgene)
  query_symbol <- unique(qgenes_match[["found"]]$symbol)

  if (length(query_symbol) > 20) {
    onc_rep[["config"]][["aberration"]][["plot_height"]] <-
      onc_rep[["config"]][["aberration"]][["plot_height"]] +
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
        otdb_all = oeDB[['otdb']][['all']],
        otdb_gene_rank = oeDB[['otdb']][['gene_rank']],
        min_association_score = 0.05)
  }

  ## Include synthetic lethality interactions
  if (show_synleth == T) {
    onc_rep[["data"]][["synleth"]] <- annotate_synleth_paralog_pairs(
      qgenes = query_symbol,
      genedb = oeDB[['genedb']][['all']],
      slparalogdb = oeDB[['slparalogdb']])
  }

  ## Include gene-drug annotations in the report (targeted cancer drugs)
  if (show_drug == T) {
    onc_rep[["data"]][["drug"]] <-
      target_drug_associations(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']])
  }

  ## Include ligand-receptor interactions in the report
  if (show_ligand_receptor == T) {
    onc_rep[["data"]][["ligand_receptor"]] <-
      annotate_ligand_receptor_interactions(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        ligand_receptor_db =
          oeDB[['ligandreceptordb']][['cellchatdb']][['db']],
        ligand_receptor_xref =
          oeDB[['ligandreceptordb']][['cellchatdb']][['xref']])
  }

  ## Include enrichment analyses in the report (pathway, GO, MSigDB)
  if (onc_rep[['config']][['show']][['enrichment']] == T) {
    for (c in names(oeDB[['pathwaydb']][["msigdb"]][["COLLECTION"]])) {
      for (subcat in names(oeDB[['pathwaydb']][["msigdb"]][["COLLECTION"]][[c]])) {
        if (c == "C5" & subcat != "HPO") {
          enr <- get_go_enrichment(
            query_entrez = as.character(query_entrezgene),
            background_entrez = background_entrezgene,
            bgset_description =
              onc_rep[["config"]][["enrichment"]][["bgset_description"]],
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
            genedb = oeDB[['genedb']][['all']])
          if (!is.null(enr)) {
            onc_rep[["data"]][["enrichment"]][["go"]] <-
              dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["go"]], enr) |>
              dplyr::distinct()
          }
        } else {
          if (c != "C5") {
            db = paste0("MSIGdb/", c, "/",subcat)
            enr <- get_universal_enrichment(
              query_entrez = as.character(query_entrezgene),
              genedb = oeDB[['genedb']][['all']],
              background_entrez = background_entrezgene,
              bgset_description =
                onc_rep[["config"]][["enrichment"]][["bgset_description"]],
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
              dbsource = db)
            if (!is.null(enr)) {
              onc_rep[["data"]][["enrichment"]][["msigdb"]] <-
                dplyr::bind_rows(onc_rep[["data"]][["enrichment"]][["msigdb"]], enr) |>
                dplyr::distinct() |>
                dplyr::arrange(.data$qvalue)
            }
          }
        }
      }
    }

    for (pwaydb in c('wikipathways','netpath','kegg')) {
      db <- "NetPath"
      if (pwaydb == "wikipathways") {
        db <- "WikiPathways"
      }
      if (pwaydb == "kegg") {
        db <- "KEGG"
      }
      onc_rep[["data"]][["enrichment"]][[pwaydb]] <-
        get_universal_enrichment(
          as.character(query_entrezgene),
          genedb = oeDB[['genedb']][['all']],
          background_entrez = background_entrezgene,
          bgset_description =
            onc_rep[["config"]][["enrichment"]][["bgset_description"]],
          min_geneset_size = onc_rep[["config"]][["enrichment"]][["min_gs_size"]],
          max_geneset_size = onc_rep[["config"]][["enrichment"]][["max_gs_size"]],
          q_value_cutoff = onc_rep[["config"]][["enrichment"]][["q_value_cutoff"]],
          p_value_cutoff = onc_rep[["config"]][["enrichment"]][["p_value_cutoff"]],
          p_value_adjustment_method =
            onc_rep[["config"]][["enrichment"]][["p_adjust_method"]],
          TERM2GENE = oeDB[['pathwaydb']][[pwaydb]]$TERM2GENE,
          TERM2NAME = oeDB[['pathwaydb']][[pwaydb]]$TERM2NAME,
          #TERM2SOURCE = oeDB[['pathwaydb']][[pwaydb]]$TERM2SOURCE,
          dbsource = db)

    }
  }

  ## Include protein-protein interactions in the report
  if (show_ppi == T) {
    ## TEST TO CHECK OMNIPATHDB IS LIVE IS NOT WORKING
    # service_is_down <- unique(is.na(pingr::ping("string-db.org")))
    onc_rep[["data"]][["ppi"]][['string']] <-
      get_ppi_network(
        qgenes = as.integer(query_entrezgene),
        ppi_source = "string",
        genedb = oeDB[['genedb']][['all']],
        cancerdrugdb = oeDB[['cancerdrugdb']],
        biogrid = oeDB[['biogrid']],
        settings = onc_rep[["config"]][["ppi"]],
        ppi_source_release = oeDB$release_notes$string$version)

    onc_rep[["data"]][["ppi"]][['biogrid']] <-
      get_ppi_network(
        qgenes = as.integer(query_entrezgene),
        ppi_source = "biogrid",
        genedb = oeDB[['genedb']][['all']],
        cancerdrugdb = oeDB[['cancerdrugdb']],
        biogrid = oeDB[['biogrid']],
        settings = onc_rep[["config"]][["ppi"]],
        ppi_source_release = oeDB$release_notes$biogrid$version)
  }

  ## include protein complex annotations in the report
  if (show_complex == T) {
    onc_rep[["data"]][["protein_complex"]] <-
      annotate_protein_complex(
        query_entrez = as.integer(query_entrezgene),
        genedb = oeDB[['genedb']][['all']],
        complex_db = oeDB[['genedb']][['proteincomplexdb']][['db']],
        complex_up_xref = oeDB[['genedb']][['proteincomplexdb']][['up_xref']],
        transcript_xref = oeDB[['genedb']][['transcript_xref']],
        otdb_gene_rank = oeDB[['otdb']][['gene_rank']])
  }

  if (show_domain == T) {
    onc_rep[["data"]][["protein_domain"]][['target']] <-
      annotate_protein_domain(
        query_entrez = as.integer(query_entrezgene),
        genedb = oeDB[['genedb']][['all']],
        pfamdb = oeDB[['pfamdb']],
        transcript_xref = oeDB[['genedb']][['transcript_xref']])
  }

  ## include annotations regarding genes of unknown/poorly defined function
  if (show_unknown_function == T) {
    onc_rep[["data"]][["unknown_function"]][["hits_df"]] <-
      get_genes_unknown_function(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']])
  }

  ## include subcellular enrichment and annotations
  if (show_subcell_comp == T) {
     subcellcomp_annotations <-
      annotate_subcellular_compartments(
        query_entrez = as.integer(query_entrezgene),
        compartments_min_confidence = subcellcomp_min_confidence,
        compartments_min_channels = subcellcomp_min_channels,
        show_cytosol = subcellcomp_show_cytosol,
        genedb = oeDB[['genedb']][['all']],
        compartments = oeDB[['subcelldb']][['compartments']],
        go_gganatogram_map = oeDB[['subcelldb']][['go_gganatogram_map']])

     onc_rep[["data"]][["subcellcomp"]][["all"]] <-
       subcellcomp_annotations[["all"]]
     onc_rep[["data"]][["subcellcomp"]][["grouped"]] <-
       subcellcomp_annotations[["grouped"]]
     onc_rep[["data"]][["subcellcomp"]][["anatogram"]] <-
       subcellcomp_annotations[["anatogram"]]

  }

  ##
  if (show_fitness == T) {
    onc_rep[["data"]][["fitness"]][["fitness_scores"]] <-
      get_fitness_lof_scores(
        qgenes = query_symbol,
        depmapdb = oeDB[['depmapdb']])

    if (onc_rep[["data"]][["fitness"]][["fitness_scores"]][["n_targets"]] <= 10) {
      onc_rep[["config"]][["fitness"]][["plot_height_fitness"]] <- 5
    }

    if (onc_rep[["data"]][["fitness"]][["fitness_scores"]][["n_targets"]] >= 20) {
      onc_rep[["config"]][["fitness"]][["plot_height_fitness"]]  <-
        onc_rep[["config"]][["fitness"]][["plot_height_fitness"]] +
        as.integer((
          onc_rep[["data"]][["fitness"]][["fitness_scores"]][["n_targets"]] - 20)/8.5)
    }

    onc_rep[["data"]][["fitness"]][["target_priority_scores"]] <-
      get_target_priority_scores(
        qgenes = query_symbol,
        depmapdb = oeDB[['depmapdb']])
  }


  if (show_aberration == T) {
    for (v in c("cna_homdel","cna_ampl")) {
      onc_rep[["data"]][["tcga"]][["aberration"]][["matrix"]][[v]] <-
        tcga_aberration_matrix(
          qgenes = as.integer(query_entrezgene),
          qsource = "entrezgene",
          genedb = oeDB[['genedb']][['all']],
          tcgadb = oeDB[['tcgadb']],
          vtype = v)
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][[v]] <-
        tcga_aberration_table(
          qgenes = as.integer(query_entrezgene),
          qsource = "entrezgene",
          genedb = oeDB[['genedb']][['all']],
          tcgadb = oeDB[['tcgadb']],
          vtype = v)
    }

    onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
      oeDB$tcgadb[["recurrent_variants"]] |>
      dplyr::inner_join(
        dplyr::select(qgenes_match$found, c("symbol")),
        by = c("SYMBOL" = "symbol"), multiple = "all") |>
      dplyr::distinct()

    if (nrow(onc_rep[["data"]][["tcga"]][["recurrent_variants"]]) > 0) {
      cosmic_variants <-
        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] |>
        dplyr::select(c("VAR_ID", "COSMIC_MUTATION_ID")) |>
        dplyr::filter(!is.na(.data$COSMIC_MUTATION_ID)) |>
        dplyr::distinct()

      if (nrow(cosmic_variants) > 0) {

        cosmic_variants <- as.data.frame(
          cosmic_variants |>
          tidyr::separate_rows("COSMIC_MUTATION_ID", sep ="&") |>
          dplyr::mutate(
            COSMIC_MUTATION_ID = paste0(
              "<a href=\"https://cancer.sanger.ac.uk/cosmic/search?q=",
              .data$COSMIC_MUTATION_ID,"\" target='_blank'>",
              .data$COSMIC_MUTATION_ID,"</a>"
            )) |>
          dplyr::group_by(.data$VAR_ID) |>
          dplyr::summarise(
            COSMIC_MUTATION_ID =
              paste(
                .data$COSMIC_MUTATION_ID, collapse = ", "
              ),
            .groups = "drop")
        )

        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
          onc_rep[["data"]][["tcga"]][["recurrent_variants"]] |>
          dplyr::select(-c("COSMIC_MUTATION_ID")) |>
          dplyr::left_join(
            cosmic_variants,
            by = c("VAR_ID"),
            multiple = "all")
      }

      onc_rep[["data"]][["tcga"]][["recurrent_variants"]] <-
        onc_rep[["data"]][["tcga"]][["recurrent_variants"]] |>
        dplyr::left_join(
          oeDB[['tcgadb']][['pfam']],
          by = "PFAM_ID",
          multiple = "all") |>
        dplyr::mutate(
          PROTEIN_DOMAIN = dplyr::if_else(
            !is.na(.data$PFAM_ID),
            paste0(
            "<a href=\"http://pfam.xfam.org/family/",
            .data$PFAM_ID,
            "\" target='_blank'>",
            .data$PFAM_DOMAIN_NAME,
            "</a>"),
            as.character(NA)
          )
        ) |>
        dplyr::select(
          -c("PFAM_DOMAIN_NAME", "PFAM_ID")) |>
        dplyr::left_join(
          dplyr::select(oeDB[['genedb']][['all']],
                        c("symbol", "ensembl_gene_id")),
          by = c("SYMBOL" = "symbol"),
          multiple = "all") |>
        dplyr::mutate(
          ENSEMBL_GENE_ID =
            paste0(
              "<a href='https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
              .data$ensembl_gene_id,"' target='_blank'>",
              .data$ensembl_gene_id,"</a>")) |>
        dplyr::mutate(
          ENSEMBL_TRANSCRIPT_ID =
            paste0(
              "<a href='https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",
              .data$ensembl_gene_id,
              ";t=",
              .data$ENSEMBL_TRANSCRIPT_ID,"' target='_blank'>",
              .data$ENSEMBL_TRANSCRIPT_ID,"</a>")) |>
        dplyr::select(-c("VAR_ID")) |>
        dplyr::rename(CONSEQUENCE_ALTERNATE = "VEP_ALL_CSQ") |>
        dplyr::select(c("SYMBOL",
                      "CONSEQUENCE",
                      "PROTEIN_CHANGE",
                      "MUTATION_HOTSPOT",
                      "PROTEIN_DOMAIN",
                      "LOSS_OF_FUNCTION",
                      "ENSEMBL_GENE_ID",
                      "ENSEMBL_TRANSCRIPT_ID",
                      "PRIMARY_SITE",
                      "SITE_RECURRENCE",
                      "TOTAL_RECURRENCE",
                      "COSMIC_MUTATION_ID",
                      "CONSEQUENCE_ALTERNATE"))
    }

    for (psite in names(onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]])) {
      onc_rep[["data"]][["tcga"]][["aberration"]][["table"]][["snv_indel"]][[psite]][['top_mutated_genes']] <-
        tcga_oncoplot_genes(
          qgenes = query_symbol,
          qsource = "symbol",
          genedb = oeDB[['genedb']][['all']],
          tcgadb = oeDB[['tcgadb']],
          site = psite)
    }
  }

  if (show_cancer_hallmarks == T) {

    lgr::lgr$info( "Open Targets Platform: Retrieving genes with evidence of cancer hallmark properties")
    onc_rep[["data"]][["cancer_hallmark"]][["target"]] <-
      oeDB[['genedb']][["cancer_hallmark"]][["short"]] |>
      dplyr::inner_join(
        dplyr::select(qgenes_match$found, c("entrezgene")),
        by = "entrezgene", multiple = "all") |>
      dplyr::select(-c("ensembl_gene_id","entrezgene")) |>
      dplyr::distinct()

    n_genes_cancerhallmarks <- 0

    if (nrow(onc_rep[["data"]][["cancer_hallmark"]][["target"]]) > 0) {
      n_genes_cancerhallmarks <-
        length(unique(onc_rep[["data"]][["cancer_hallmark"]][["target"]]$symbol))
    }
    lgr::lgr$info(
      paste0("Number of query genes attributed with cancer hallmark properties: ",
             n_genes_cancerhallmarks))

  }

  if (show_coexpression == T) {
    onc_rep[["data"]][["tcga"]][["coexpression"]] <-
      tcga_coexpression(
        qgenes = query_symbol,
        genedb = oeDB[['genedb']][['all']],
        tcgadb = oeDB[['tcgadb']])
  }

  if (show_regulatory == T) {

    for (collection in c('global','pancancer')) {
      onc_rep[["data"]][["regulatory"]][["interactions"]][[collection]] <-
        annotate_tf_targets(
          query_symbol,
          genedb = oeDB[['genedb']][['all']],
          tf_target_interactions = oeDB[['tftargetdb']],
          collection = collection,
          regulatory_min_confidence =
            regulatory_min_confidence)
    }

    if (NROW(onc_rep[["data"]][["regulatory"]][["interactions"]][["pancancer"]]) > 0) {
      onc_rep[["data"]][["regulatory"]][["network"]] <-
        retrieve_tf_target_network(
          tf_target_interactions =
            onc_rep[["data"]][["regulatory"]][["interactions"]][["pancancer"]]
        )
    }
  }

  if (show_prognostic == T) {
    onc_rep[["data"]][["cancer_prognosis"]][['hpa']][['assocs']] <-
      hpa_prognostic_genes(
        query_symbol,
        genedb = oeDB[['genedb']][['all']],
        hpadb = oeDB[['hpa']])

    for (feature in c('exp', 'mut', 'cna', 'meth')) {
      onc_rep[["data"]][["cancer_prognosis"]][['km_cshl']][['assocs']][[feature]] <-
        km_cshl_survival_genes(
          query_symbol,
          survivaldb = oeDB[['survivaldb']][[feature]],
          genetic_feature = feature)
    }

  }

  if (show_cell_tissue == T) {
    onc_rep[["data"]][["cell_tissue"]][['tissue_overview']] <-
      gene_tissue_cell_spec_cat(
        qgenes = query_symbol,
        q_id_type = "symbol",
        resolution = "tissue",
        genedb = oeDB[['genedb']][['all']],
        hpa_enrichment_db_df =
          oeDB[['tissuecelldb']][['tissue']][['te_df']],
        hpa_expr_db_df =
          oeDB[['tissuecelldb']][['tissue']][['expr_df']])

    onc_rep[["data"]][["cell_tissue"]][['tissue_enrichment']] <-
      gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(query_entrezgene),
        resolution = "tissue",
        background_entrez = as.integer(background_entrezgene),
        genedb = oeDB[['genedb']][['all']],
        hpa_enrichment_db_df =
          oeDB[['tissuecelldb']][['tissue']][['te_df']],
        hpa_enrichment_db_SE =
          oeDB[['tissuecelldb']][['tissue']][['te_SE']])

    onc_rep[["data"]][["cell_tissue"]][['scRNA_overview']] <-
      gene_tissue_cell_spec_cat(
        qgenes = as.integer(query_entrezgene),
        q_id_type = "entrezgene",
        resolution = "single_cell",
        genedb = oeDB[['genedb']][['all']],
        hpa_enrichment_db_df =
          oeDB[['tissuecelldb']][['single_cell']][['te_df']],
        hpa_expr_db_df =
          oeDB[['tissuecelldb']][['single_cell']][['expr_df']])

    onc_rep[["data"]][["cell_tissue"]][['scRNA_enrichment']] <-
      gene_tissue_cell_enrichment(
        qgenes_entrez = as.integer(query_entrezgene),
        background_entrez = as.integer(background_entrezgene),
        resolution = "single_cell",
        genedb = oeDB[['genedb']][['all']],
        hpa_enrichment_db_df =
          oeDB[['tissuecelldb']][['single_cell']][['te_df']],
        hpa_enrichment_db_SE =
          oeDB[['tissuecelldb']][['single_cell']][['te_SE']])

  }

  return(onc_rep)

}

#' Write oncoEnrichR report object to output file
#'
#' Function that writes an oncoEnrichR report object to file,
#' either as an interactive HTML report or as an Excel workbook.

#'
#' @param report object with oncoEnrichR report data (returned by oeDB$onco_enrich)
#' @param oeDB oncoEnrichR data repository object - as returned from `load_db()`
#' @param file full filename for report output (e.g. "oe_report.html" or "oe_report.xlsx")
#' @param ignore_file_extension logical to accept any type of filaname extensions (for Galaxy integration)
#' @param overwrite logical indicating if existing output files may be overwritten
#' @param format file format of output (html/excel)
#' @param ... options for Galaxy/non self-contained HTML. Only applicable for use in Galaxy
#'
#'
#' @export

write <- function(report,
                  oeDB,
                  file = "testReport.html",
                  ignore_file_extension = F,
                  overwrite = F,
                  format = "html",
                  ...) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  selfcontained <- T
  galaxy_run <- T
  html_extern_path <- NA
  # logger <- log4r::logger(threshold = "INFO",
  #                         appenders = log4r::console_appender(log4r_layout))
  dot_args <- list(...)
  if (length(names(dot_args)) > 0) {

    for (arg in names(dot_args)) {
      if (!(arg == "selfcontained_html" |
            arg == "galaxy" | arg == "extra_files_path")) {
        lgr::lgr$info( paste0(
          "ERROR: argument '",arg,"' does not exist for oeDB$write()"))
        return()
      }
    }

    if ("selfcontained_html" %in% names(dot_args)) {
      if (is.logical(dot_args$selfcontained_html)) {
        selfcontained <- dot_args$selfcontained_html
      }
    }
    if ("galaxy" %in% names(dot_args)) {
      if (is.logical(dot_args$galaxy)) {
        galaxy_run <- dot_args$galaxy
      }
    }
    if ("extra_files_path" %in% names(dot_args)) {
      if (is.character(dot_args$extra_files_path)) {
        html_extern_path <- dot_args$extra_files_path

        if (!dir.exists(html_extern_path)) {
          dir.create(html_extern_path)
        }
      }
    }
  }

  val <- assertthat::validate_that(
    format %in% c("html","excel")
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: ",val))
    return()
  }

  val <- assertthat::validate_that(
    is.character(file)
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: ",val))
    return()
  }

  val <- assertthat::validate_that(
    is.logical(overwrite)
  )
  if (!is.logical(val)) {
    lgr::lgr$info( paste0(
      "ERROR: ",val))
    return()
  }

  output_directory <- dirname(file)
  file_basename <- basename(file)
  file_basename_prefix <- stringr::str_replace(
    file_basename,"\\.(html|xlsx)$","")

  if (output_directory != ".") {
    val <- assertthat::validate_that(
      dir.exists(output_directory)
    )
    if (!is.logical(val)) {
      lgr::lgr$info(paste0("ERROR: ",val))
      return()
    }
  }

  if (overwrite == F) {
    val <- assertthat::validate_that(
      file.exists(file)
    )
    if (is.logical(val)) {
      lgr::lgr$info(
        paste0("ERROR: output file (",file, ") exists, (overwrite = F)"))
      return()
    }
  }

  if (ignore_file_extension == F) {
    if (format == "html" & tools::file_ext(file) != "html") {
      lgr::lgr$info( paste0(
        "ERROR: oncoEnrichR HTML output: File name must end with .html, not ",
        tools::file_ext(file))
      )
      return()
    }

    if (format == "excel" & tools::file_ext(file) != "xlsx") {
      lgr::lgr$info( paste0(
        "ERROR: oncoEnrichR Excel output: File name must end with .xlsx, not ",
        tools::file_ext(file))
      )
      return()
    }
  }

  ## Assign to env
  pos <- 1
  envir = as.environment(pos)
  #for (e in export) assign(e, get(e), envir = envir)

  ## TODO: check that report parameter is a valid oncoEnrichR result object
  if (!is.null(report)) {
    # assign("onc_enrich_report",
    #        report, envir = .GlobalEnv)
    assign("onc_enrich_report",
           report,
           envir = envir)
  } else {
    lgr::lgr$info(
      "ERROR: oncoEnrichR report object is NULL - cannot write report contents")
    return()
  }



  if (!is.null(oeDB[['tcgadb']][['maf']])) {
    # assign("tcga_maf_datasets",
    #        oeDB[['tcgadb']][['maf']], envir = .GlobalEnv)
    assign("tcga_maf_datasets",
           oeDB[['tcgadb']][['maf']],
           envir = envir)
  } else {
    lgr::lgr$info(
      "ERROR: oeDB$tcgadb$maf is NULL - cannot write report contents")
    return()
  }




  if (format == "html") {

    lgr::lgr$info( "------")
    lgr::lgr$info( "Writing HTML file with report contents")

    if (selfcontained == F) {

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
      dir.create(tmpdir)
      #system(paste0('mkdir ', tmpdir))
      system(paste0('cp ',
                    oe_rmarkdown_template_dir,
                    .Platform$file.sep,
                    "* ",
                    tmpdir))

      floating_toc <- "false"
      if (report$config$rmarkdown$floating_toc == T) {
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
          paste0("    theme: ", report$config$rmarkdown$theme),
          "    includes:",
          "      in_header: _header.html",
          "      after_body: _disclaimer.md", sep ="\n")
      sink()


      rmdown_html <- file.path(tmpdir,"_site", "index.html")
      rmdown_supporting1 <- file.path(tmpdir,"_site","index_files")
      rmdown_supporting2 <- file.path(tmpdir,"_site","site_libs")

      if (galaxy_run == T) {

        if (dir.exists(file.path(
          html_extern_path, "site_libs"))) {
          if (overwrite == F) {
            lgr::lgr$info( paste0(
              "ERROR: Cannot create HTML since 'site_libs' exist in output_directory",
              " and 'overwrite' is FALSE")
            )
            return()
          } else {
            system(paste0('rm -rf ',file.path(
              html_extern_path, "site_libs"
            )))
          }
        }
        if (dir.exists(file.path(
          html_extern_path, "index_files"))) {
          if (overwrite == F) {
            lgr::lgr$info( paste0(
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

        suppressWarnings(
          rmarkdown::render_site(
            input = tmpdir,
            quiet = T
          )
        )

        # target_html <- file.path(output_directory, paste0(
        #   file_basename_prefix, ".html")
        # )

        if (file.exists(rmdown_html) & dir.exists(rmdown_supporting1) &
           dir.exists(rmdown_supporting2)) {
          system(paste0('mv ', rmdown_html, ' ',
                        file))
          system(paste0('mv ', rmdown_supporting1," ", html_extern_path))
          system(paste0('mv ', rmdown_supporting2," ", html_extern_path))
          system(paste0('rm -rf ',tmpdir))

          lgr::lgr$info( paste0("Output file: ",
                                          file))
          lgr::lgr$info( "------")
        }
      }

    } else {

      disclaimer <- system.file(
        "templates",
        "_disclaimer.md",
        package = "oncoEnrichR")

      header <- system.file(
        "templates",
        "_header.html",
        package = "oncoEnrichR")

      report_theme <- report$config$rmarkdown$theme
      toc_float <- report$config$rmarkdown$floating_toc

      markdown_input <- system.file(
        "templates",
        "index.Rmd",
        package = "oncoEnrichR")

      suppressWarnings(
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
            includes = rmarkdown::includes(
              in_header = header,
              after_body = disclaimer)),
          output_file = file_basename,
          output_dir = output_directory,
          clean = T,
          intermediates_dir = output_directory,
          quiet = T)
      )

      lgr::lgr$info( paste0("Output file (self-contained HTML): ",
                                      file))
      lgr::lgr$info( "------")
    }


  }
  if (format == "excel") {

    wb <- openxlsx::createWorkbook()
    lgr::lgr$info( "------")
    lgr::lgr$info( "Writing Excel workbook with report contents")

    table_style_index <- 15
    for (elem in c("settings",
                  "query",
                  "unknown_function",
                  "cancer_association",
                  "cancer_hallmark",
                  "drug_known",
                  "drug_tractability",
                  "synthetic_lethality",
                  "fitness_scores",
                  "fitness_prioritized",
                  "protein_complex",
                  "protein_domain",
                  "ppi_string",
                  "ppi_biogrid",
                  "enrichment",
                  "regulatory",
                  "ligand_receptor",
                  "subcellcomp",
                  "cell_tissue",
                  "aberration",
                  "recurrent_variants",
                  "coexpression",
                  "prognostic_association_I",
                  "prognostic_association_II"
                  )) {

      show_elem <- elem
      if (elem == "cancer_association") {
        show_elem <- "disease"
      }

      if (elem == "recurrent_variants") {
        show_elem <- "aberration"
      }
      if (elem == "ppi_string"){
        show_elem <- "ppi"
      }
      if(elem == "ppi_biogrid"){
        show_elem <- "ppi"
      }

      if (elem == "prognostic_association_I") {
        show_elem <- "cancer_prognosis"
      }
      if (elem == "prognostic_association_II") {
        show_elem <- "cancer_prognosis"
      }
      if (elem == "drug_known") {
        show_elem <- "drug"
      }
      if (elem == "synthetic_lethality") {
        show_elem <- "synleth"
      }
      if (elem == "fitness_scores" | elem == "fitness_prioritized") {
        show_elem <- "fitness"
      }

      if (elem != "settings") {
        if (report[['config']][['show']][[show_elem]] == FALSE) {
          next
        }
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
      if (table_style_index == 21) {
        table_style_index <- 14
      }
      table_style_index <- table_style_index + 1

    }

    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    lgr::lgr$info( paste0("Output file: ",file))
    lgr::lgr$info( "------")
  }
  else{
    if (format == "json") {
      lgr::lgr$info( "JSON output not yet implemented")
    }
  }
}
