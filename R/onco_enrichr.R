
init_report <- function(title = "Project Title",
                        project_owner = "Project Owner",
                        ppi_min_string_score = 900,
                        ppi_add_nodes = 50,
                        background_enrichment_description = "All protein-coding genes",
                        p_value_cutoff_enrichment = 0.05,
                        q_value_cutoff_enrichment = 0.01,
                        show_ppi = T,
                        show_disease = T,
                        show_enrichment = T,
                        show_tcga_aberration = T,
                        show_tcga_coexpression = T,
                        #show_gtex_coexp = F,
                        show_complex = T){

  ## report object
  rep <- list()

  ## two main elements
  # 1. data - contains all annotations and enrichment results
  # 2. config - contains all settings and underlying data sources
  for(e in c('data','config')){
    rep[[e]] <- list()
  }

  ## release notes - software and database versions
  rep[['config']][['resources']] <-  oncoEnrichR::release_notes

  ## logicals indicating which sections/analyses of the report to include
  rep[['config']][['show']] <- list()
  rep[['config']][['show']][['ppi']] <- show_ppi
  rep[['config']][['show']][['target_disease']] <- show_disease
  rep[['config']][['show']][['enrichment']] <- show_enrichment
  rep[['config']][['show']][['protein_complex']] <- show_complex
  rep[['config']][['show']][['tcga_aberration']] <- show_tcga_aberration
  rep[['config']][['show']][['tcga_coexpression']] <- show_tcga_coexpression
  #rep[['config']][['show']][['gtex_coexp']] <- show_gtex_coexp

  ## report metadata - project owner, background, title
  rep[['config']][['title']] <- title
  rep[['config']][['project_background']] <- list()
  rep[['config']][['project_background']][['owner']] <- project_owner
  rep[['config']][['project_background']][['items']] <- data.frame()

  ## settings for showing disease associations
  rep[['config']][['target_disease']] <- list()
  rep[['config']][['target_disease']][['breaks']] <- c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  rep[['config']][['target_disease']][['colors']] <- c("#b8b8ba","#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#084594")

  ## protein-protein interaction settings
  rep[['config']][['ppi']] <- list()
  rep[['config']][['ppi']][['stringdb']] <- list()
  rep[['config']][['ppi']][['stringdb']][['minimum_score']] <- ppi_min_string_score
  rep[['config']][['ppi']][['stringdb']][['visnetwork_shape']] <- 'dot'
  rep[['config']][['ppi']][['stringdb']][['visnetwork_shadow']] <- T
  rep[['config']][['ppi']][['stringdb']][['add_nodes']] <- ppi_add_nodes
  rep[['config']][['ppi']][['stringdb']][['query_type']] <- 'network'

  ## enrichment settings
  rep[['config']][['enrichment']] <- list()
  rep[['config']][['enrichment']][['p_value_cutoff']] <- p_value_cutoff_enrichment
  rep[['config']][['enrichment']][['q_value_cutoff']] <- q_value_cutoff_enrichment
  rep[['config']][['enrichment']][['p_adjust_method']] <- 'BH'
  rep[['config']][['enrichment']][['min_gs_size']] <- 5
  rep[['config']][['enrichment']][['background_set']] <- background_enrichment_description

  ## TCGA/GTex
  rep[['config']][['co_expression_gtex']] <- list()
  rep[['config']][['co_expression_gtex']][['plot_height']] <- 5
  rep[['config']][['tcga_aberration']] <- list()
  rep[['config']][['tcga_aberration']][['plot_height']] <- 13

  rep[['data']][['target_disease']] <- data.frame()
  rep[['data']][['ppi']] <- list()
  rep[['data']][['ppi']][['stringdb']] <- list()
  rep[['data']][['ppi']][['stringdb']][['network']] <- NULL
  rep[['data']][['enrichment']] <- list()
  for(c in c('go','msigdb','wikipathwaydb','keggdb')){
    rep[['data']][['enrichment']][[c]] <- data.frame()
  }

  rep[['data']][['co_expression_gtex']] <- list()
  rep[['data']][['co_expression_gtex']][['plots']] <- NULL
  rep[['data']][['co_expression_gtex']][['df']] <- data.frame()
  rep[['data']][['co_expression_tcga']] <- data.frame()
  rep[['data']][['protein_complex']] <- data.frame()
  rep[['data']][['tcga_aberration']] <- list()
  rep[['data']][['tcga_aberration']][['table']] <- list()
  rep[['data']][['tcga_aberration']][['plot']] <- list()

  for(v in c('cna_homdel','cna_ampl')){
    rep[['data']][['tcga_aberration']][['plot']][[v]] <- NULL
    rep[['data']][['tcga_aberration']][['table']][[v]] <- data.frame()
  }

  return(rep)
}

#' Function that interrogates a list of human protein-coding genes for known disease associations, protein-protein interactions,
#' enrichment in gene ontology (GO), pathway databases and curated signatures, and aberration frequencies and co-expression patterns in tumor samples.
#'
#' @param query character vector with gene identifiers
#' @param query_source character indicating source of query ('uniprot_acc','symbol','entrezgene','ensembl_gene_id')
#' @param p_title title of report
#' @param p_owner name of project owner
#' @param report_fname filename prefix for HTML report
#' @param background_fname filename for simple text file with project background information, one line per background item
#' @param ppi_add_nodes number of nodes to add to protein-protein interaction network
#' @param background_enrichment character vector with gene identifiers, used as reference/background for enrichment analysis
#' @param background_enrichment_source character indicating source of background ('uniprot_acc','symbol','entrezgene','ensembl_gene_id')
#' @param background_enrichment_description character indicating type of background (e.g. 'All lipid-binding proteins (n = 200)')
#' @param p_value_cutoff_enrichment cutoff p-value for enrichment analysis
#' @param q_value_cutoff_enrichment cutoff q-value for enrichment analysis
#' @param show_ppi logical indicating if report should contain protein-protein interaction network (STRING)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform)
#' @param show_enrichment logical indicating if report should contain functional enrichment analysis (MSigDB, GO, KEGG, REACTOME etc.)
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexpression logical indicating if report should contain TCGA co-expression data (RNAseq) of queryset with oncogenes/tumor suppressor genes
#' @param show_complex logical indicating if report should list proteins in known protein complexes
#' @param report_name filename for report
#' @param format file format of output (html/json)
#' @export
#'
generate_report_data <- function(query,
                            query_source = "symbol",
                            p_title = "Project Title",
                            p_owner = "Project Owner",
                            report_fname = "oncoEnrichR_Report",
                            background_fname = NULL,
                            ppi_add_nodes = 50,
                            background_enrichment = NULL,
                            background_enrichment_source = "symbol",
                            background_enrichment_description = "All protein-coding genes",
                            p_value_cutoff_enrichment = 0.05,
                            q_value_cutoff_enrichment = 0.01,
                            show_ppi = T,
                            show_disease = T,
                            show_enrichment = T,
                            show_tcga_aberration = T,
                            show_tcga_coexpression = T,
                            #show_gtex_coexp = F,
                            show_complex = T){
                            #gtex_atlasassay_groups = c("g32","g9","g29","g10","g28","g44","g33","g50","g37","g38","g42","g35")){
  stopifnot(is.character(query))
  stopifnot(query_source == "symbol" | query_source == "entrezgene" |
              query_source == "uniprot_acc" | query_source == "ensembl_gene_id" |
              query_source == "any")

  qgenes_match <-
    oncoEnrichR::verify_query_genes(query,
                                    qsource = query_source,
                                    genedb = oncoEnrichR::genedb,
                                    uniprot_acc = oncoEnrichR::uniprot_xref)


  background_entrez <- NULL
  background_genes_match <- NULL
  if(!is.null(background_enrichment)){
    background_genes_match <-
      oncoEnrichR::verify_query_genes(background_enrichment,
                                      qsource = background_enrichment_source,
                                      genedb = oncoEnrichR::genedb,
                                      uniprot_acc = oncoEnrichR::uniprot_xref)
    if(background_genes_match[['success']] == -1){
      return(-1)
    }
    background_entrez <- background_genes_match[['found']]$entrezgene
  }

  if(qgenes_match[['success']] == -1){
    return(-1)
  }

  query_entrezgene <- qgenes_match[['found']]$entrezgene
  query_symbol <- qgenes_match[['found']]$symbol

  onc_rep <- init_report(title = p_title,
                                        project_owner = p_owner,
                                        ppi_add_nodes = ppi_add_nodes,
                                        show_ppi = show_ppi,
                                        show_disease = show_disease,
                                        show_enrichment = show_enrichment,
                                        show_tcga_aberration = show_tcga_aberration,
                                        show_tcga_coexpression = show_tcga_coexpression,
                                        #show_gtex_coexp = show_gtex_coexp,
                                        show_complex = show_complex,
                                        background_enrichment_description = background_enrichment_description,
                                        p_value_cutoff_enrichment = p_value_cutoff_enrichment,
                                        q_value_cutoff_enrichment = q_value_cutoff_enrichment)

  if(length(query_symbol) > 30){
    onc_rep[['config']][['tcga_aberration']][['plot_height']] <-
      onc_rep[['config']][['tcga_aberration']][['plot_height']] + as.integer((length(query_symbol) - 30)/ 10)
    onc_rep[['config']][['co_expression_gtex']][['plot_height']] <-
      onc_rep[['config']][['co_expression_gtex']][['plot_height']] + as.integer((length(query_symbol) - 30)/ 10)
  }

  if(!is.null(background_fname)){
    if(file.exists(background_fname)){
      project_bg <- read.table(file=background_fname,stringsAsFactors = F,
                             header=F,sep="\n",na.strings="",comment.char="")
      colnames(project_bg) <- 'text'
      onc_rep[['config']][['project_background']][['items']] <- project_bg
    }
  }


  if(show_disease == T){
    onc_rep[['data']][['target_disease']] <-
      oncoEnrichR::target_disease_associations(query_symbol, genedb = oncoEnrichR::genedb)
  }

  for(c in names(oncoEnrichR::msigdb[['COLLECTION']])){
    for(subcat in names(oncoEnrichR::msigdb[['COLLECTION']][[c]])){
      if(c == "C5"){
        enr <- oncoEnrichR::get_go_enrichment(query_entrezgene,
                                  background_entrez = background_entrez,
                                  minGSSize = onc_rep[['config']][['enrichment']][['min_gs_size']],
                                  q_value_cutoff = onc_rep[['config']][['enrichment']][['q_value_cutoff']],
                                  p_value_cutoff = onc_rep[['config']][['enrichment']][['p_value_cutoff']],
                                  ontology = subcat,
                                  genedb = oncoEnrichR::genedb)
        if(!is.null(enr)){
          onc_rep[['data']][['enrichment']][['go']] <-
            dplyr::bind_rows(onc_rep[['data']][['enrichment']][['go']], enr) %>%
            dplyr::distinct()
        }
      }else{
        db = paste0("MSIGdb/",c,"/",subcat)
        enr <- oncoEnrichR::get_universal_enrichment(query_entrezgene,
                                         genedb = oncoEnrichR::genedb,
                                         background_entrez = background_entrez,
                                         minGSSize = onc_rep[['config']][['enrichment']][['min_gs_size']],
                                         q_value_cutoff = onc_rep[['config']][['enrichment']][['q_value_cutoff']],
                                         p_value_cutoff = onc_rep[['config']][['enrichment']][['p_value_cutoff']],
                                         TERM2GENE = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2GENE,
                                         TERM2NAME = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2NAME,
                                         TERM2SOURCE = oncoEnrichR::msigdb$TERM2SOURCE,
                                         dbsource = db)
        if(!is.null(enr)){
          onc_rep[['data']][['enrichment']][['msigdb']] <-
            dplyr::bind_rows(onc_rep[['data']][['enrichment']][['msigdb']], enr) %>%
            dplyr::distinct()
        }
      }
    }
  }


  db <- 'WikiPathways'
  onc_rep[['data']][['enrichment']][['wikipathwaydb']] <-
    oncoEnrichR::get_universal_enrichment(query_entrezgene,
                              genedb = oncoEnrichR::genedb,
                              background_entrez = background_entrez,
                              minGSSize = onc_rep[['config']][['enrichment']][['min_gs_size']],
                              q_value_cutoff = onc_rep[['config']][['enrichment']][['q_value_cutoff']],
                              p_value_cutoff = onc_rep[['config']][['enrichment']][['p_value_cutoff']],
                              TERM2GENE = oncoEnrichR::wikipathwaydb$TERM2GENE,
                              TERM2NAME = oncoEnrichR::wikipathwaydb$TERM2NAME,
                              TERM2SOURCE = oncoEnrichR::wikipathwaydb$TERM2SOURCE,
                              dbsource = db)

  db <- 'KEGG'
  onc_rep[['data']][['enrichment']][['keggdb']] <-
    oncoEnrichR::get_universal_enrichment(query_entrezgene,
                                          genedb = oncoEnrichR::genedb,
                                          background_entrez = background_entrez,
                                          minGSSize = onc_rep[['config']][['enrichment']][['min_gs_size']],
                                          q_value_cutoff = onc_rep[['config']][['enrichment']][['q_value_cutoff']],
                                          p_value_cutoff = onc_rep[['config']][['enrichment']][['p_value_cutoff']],
                                          TERM2GENE = oncoEnrichR::keggdb$TERM2GENE,
                                          TERM2NAME = oncoEnrichR::keggdb$TERM2NAME,
                                          TERM2SOURCE = oncoEnrichR::keggdb$TERM2SOURCE,
                                          dbsource = db)


  # if(show_gtex_coexp == T){
  #   gtex_results <-
  #     oncoEnrichR::gtex_co_expression(query_symbol, genedb = oncoEnrichR::genedb, gids = gtex_atlasassay_groups)
  #   onc_rep[['data']][['co_expression_gtex']][['plots']] <- gtex_results[['plots']]
  #   onc_rep[['data']][['co_expression_gtex']][['df']] <- gtex_results[['df']]
  # }

  if(show_ppi == T){
    onc_rep[['data']][['ppi']][['stringdb']][['network']] <-
      oncoEnrichR::get_string_ppi_network(query_entrezgene,
                                          genedb = oncoEnrichR::genedb,
                                          settings = onc_rep[['config']][['ppi']][['stringdb']])
  }

  if(show_complex == T){
    onc_rep[['data']][['protein_complex']] <-
      oncoEnrichR::annotate_protein_complex(query_symbol,
                                            genedb = oncoEnrichR::genedb,
                                            corum_db = oncoEnrichR::corumdb,
                                            uniprot_acc = oncoEnrichR::uniprot_xref)
  }

  if(show_tcga_aberration == T){
    for(v in c('cna_homdel','cna_ampl')){
      onc_rep[['data']][['tcga_aberration']][['plot']][[v]] <-
        oncoEnrichR::tcga_aberration_plot(query_entrezgene, qsource = "entrezgene", genedb = oncoEnrichR::genedb, vtype = v)
      onc_rep[['data']][['tcga_aberration']][['table']][[v]] <-
        oncoEnrichR::tcga_aberration_table(query_entrezgene, qsource = "entrezgene", genedb = oncoEnrichR::genedb, vtype = v)
    }
  }

  if(show_tcga_coexpression == T){
    onc_rep[['data']][['co_expression_tcga']] <-
      oncoEnrichR::tcga_co_expression(query_symbol, genedb = oncoEnrichR::genedb)
  }

  return(onc_rep)

}

#' Function that writes the contents in the oncoEnrichR report object to an interactive HTML report
#'
#' @param report object with oncoEnrichR report data (returned by oncoEnrichR::generate_report)
#' @param project_directory working directory
#' @param report_name filename for report
#' @param format file format of output (html/json)
#' @export

write_report <- function(report, project_directory, report_name, format = 'html'){

  outfname <- list()
  outfname[['html']] <- paste(report_name,"html",sep=".")
  disclaimer <- system.file("templates", "disclaimer.md", package = "oncoEnrichR")
  report_theme <- "default"

  assign("onc_enrich_report", report, envir = .GlobalEnv)

  if(format == "html"){
    rlogging::message("------")
    rlogging::message("Writing HTML file with report contents")
    markdown_input <- system.file("templates", "onco_enrich_report.Rmd", package = "oncoEnrichR")
    rmarkdown::render(markdown_input,
                      output_format = rmarkdown::html_document(theme = report_theme, toc = T, toc_depth = 3,
                                                               toc_float = T, number_sections = F,
                                                               includes = rmarkdown::includes(after_body = disclaimer)),
                      output_file = outfname[['html']],
                      output_dir = project_directory,
                      clean = T,
                      intermediates_dir = project_directory,
                      quiet = T)
  }else{
    rlogging::message('JSON output not yet implemented')
  }
}
