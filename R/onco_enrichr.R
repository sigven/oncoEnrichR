
init_report <- function(title = "Project Title",
                        project_owner = "Project Owner",
                        ppi_min_string_score = 900,
                        ppi_add_nodes = 50,
                        background_enrichment_description = "All protein-coding genes",
                        p_value_cutoff = 0.05,
                        q_value_cutoff = 0.01,
                        show_ppi = T,
                        show_disease = T,
                        show_enrichment = T,
                        show_tcga_aberration = T,
                        show_tcga_coexp = T,
                        show_gtex_coexp = F,
                        show_complex = T){

  oncr_report <- list()
  oncr_report[['show_section']] <- list()
  oncr_report[['show_section']][['ppi']] <- show_ppi
  oncr_report[['show_section']][['target_disease']] <- show_disease
  oncr_report[['show_section']][['enrichment']] <- show_enrichment
  oncr_report[['show_section']][['complex']] <- show_complex
  oncr_report[['show_section']][['tcga_aberration']] <- show_tcga_aberration
  oncr_report[['show_section']][['tcga_coexp']] <- show_tcga_coexp
  oncr_report[['show_section']][['gtex_coexp']] <- show_gtex_coexp

  oncr_report[['title']] <- title
  oncr_report[['project_background']] <- list()
  oncr_report[['project_background']][['owner']] <- project_owner
  oncr_report[['project_background']][['items']] <- data.frame()
  oncr_report[['target_disease']] <- data.frame()
  oncr_report[['disease_rank']] <- list()
  oncr_report[['disease_rank']][['breaks']] <- c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  oncr_report[['disease_rank']][['colors']] <- c("#b8b8ba","#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#084594")
  oncr_report[['ppi']] <- list()
  oncr_report[['ppi']][['stringdb']] <- list()
  oncr_report[['ppi']][['stringdb']][['settings']] <- list()
  oncr_report[['ppi']][['stringdb']][['network']] <- NULL
  oncr_report[['ppi']][['stringdb']][['settings']][['minimum_score']] <- ppi_min_string_score
  oncr_report[['ppi']][['stringdb']][['settings']][['visnetwork_shape']] <- 'dot'
  oncr_report[['ppi']][['stringdb']][['settings']][['visnetwork_shadow']] <- T
  oncr_report[['ppi']][['stringdb']][['settings']][['add_nodes']] <- ppi_add_nodes
  oncr_report[['ppi']][['stringdb']][['settings']][['query_type']] <- 'network'

  oncr_report[['enrichment']] <- list()
  oncr_report[['enrichment']][['settings']] <- list()
  oncr_report[['enrichment']][['settings']][['p_value_cutoff']] <- p_value_cutoff
  oncr_report[['enrichment']][['settings']][['q_value_cutoff']] <- q_value_cutoff
  oncr_report[['enrichment']][['settings']][['p_adjust_method']] <- 'BH'
  oncr_report[['enrichment']][['settings']][['min_gs_size']] <- 5
  oncr_report[['enrichment']][['settings']][['background_set']] <- background_enrichment_description
  oncr_report[['enrichment']][['results']] <- list()
  for(c in c('go','msigdb','wikipathwaydb','acsn')){
    oncr_report[['enrichment']][['results']][[c]] <- data.frame()
  }
  oncr_report[['co_expression_gtex']] <- list()
  oncr_report[['co_expression_gtex']][['plot_height']] <- 5
  oncr_report[['co_expression_gtex']][['plots']] <- NULL
  oncr_report[['co_expression_gtex']][['df']] <- data.frame()
  oncr_report[['protein_complex']] <- data.frame()
  oncr_report[['tcga_aberration']] <- list()
  oncr_report[['tcga_aberration']][['table']] <- list()
  oncr_report[['tcga_aberration']][['plot']] <- list()
  oncr_report[['tcga_aberration']][['plot_height']] <- 13
  for(v in c('cna_homdel','cna_ampl')){
    oncr_report[['tcga_aberration']][['plot']][[v]] <- NULL
    oncr_report[['tcga_aberration']][['table']][[v]] <- data.frame()
  }

  return(oncr_report)
}

#' Function that interrogates a list of human protein-coding genes for known disease associations, protein-protein interactions,
#' enrichment in gene ontology (GO), pathway databases and curated signatures, and aberration frequencies and co-expression patterns in tumor samples.
#'
#' @param query character vector with gene identifiers
#' @param query_source character indicating source of query ('uniprot_acc','symbol','entrezgene','ensembl_gene_id')
#' @param p_title title of report
#' @param p_owner name of project owner
#' @param report_fname filename prefix for HTML report
#' @param background_fname filename for text file with project background information
#' @param ppi_add_nodes number of nodes to add to protein-protein interaction network
#' @param background_enrichment character vector with gene identifiers, used as reference/background for enrichment analysis
#' @param background_enrichment_source character indicating source of background ('uniprot_acc','symbol','entrezgene','ensembl_gene_id')
#' @param background_enrichment_description character indicating type of background (e.g. 'All lipid-binding proteins (n = 200)')
#' @param show_ppi logical indicating if report should contain protein-protein interaction network (STRING)
#' @param show_disease logical indicating if report should contain disease associations (Open Targets Platform)
#' @param show_enrichment logical indicating if report should contain functional enrichment analysis (MSigDB (GO, KEGG, REACTOME etc.))
#' @param show_tcga_aberration logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
#' @param show_tcga_coexp logical indicating if report should contain co-expression data of queryset with oncogenes/tumor suppressor genes (TCGA)
#' @param show_complex logical indicating if report should list proteins in known protein complexes
#' @param report_name filename for report
#' @param format file format of output (html/json)
#' @export
#'
generate_report <- function(query,
                            query_source = "symbol",
                            p_title = "Project Title",
                            p_owner = "Project Owner",
                            report_fname = "OncoEnrichR_TestReport",
                            background_fname = NULL,
                            ppi_add_nodes = 50,
                            background_enrichment = NULL,
                            background_enrichment_source = "symbol",
                            background_enrichment_description = "All protein-coding genes",
                            show_ppi = T,
                            show_disease = T,
                            show_enrichment = T,
                            show_tcga_aberration = T,
                            show_tcga_coexp = T,
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
                                    uniprot_acc = oncoEnrichR::uniprot_acc)


  background_entrez <- NULL
  background_genes_match <- NULL
  if(!is.null(background_enrichment)){
    background_genes_match <-
      oncoEnrichR::verify_query_genes(background_enrichment,
                                      qsource = background_enrichment_source,
                                      genedb = oncoEnrichR::genedb,
                                      uniprot_acc = oncoEnrichR::uniprot_acc)
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


  ge_report <- oncoEnrichR::init_report(title = p_title, project_owner = p_owner,
                                        ppi_add_nodes = ppi_add_nodes,
                                        show_ppi = show_ppi,
                                        show_disease = show_disease,
                                        show_enrichment = show_enrichment,
                                        show_tcga_aberration = show_tcga_aberration,
                                        show_tcga_coexp = show_tcga_coexp,
                                        show_gtex_coexp = show_gtex_coexp,
                                        show_complex = show_complex,
                                        background_enrichment_description = background_enrichment_description)

  if(length(query_symbol) > 30){
    ge_report[['tcga_aberration']][['plot_height']] <-
      ge_report[['tcga_aberration']][['plot_height']] + as.integer((length(query_symbol) - 30)/ 10)
    ge_report[['co_expression_gtex']][['plot_height']] <-
      ge_report[['co_expression_gtex']][['plot_height']] + as.integer((length(query_symbol) - 30)/ 10)
  }

  if(!is.null(background_fname) & file.exists(background_fname)){
    project_bg <- read.table(file=background_fname,stringsAsFactors = F,
                             header=F,sep="\n",na.strings="",comment.char="")
    colnames(project_bg) <- 'text'
    ge_report$project_background$items <- project_bg
  }


  if(show_disease == T){
    ge_report[['target_disease']] <-
      oncoEnrichR::target_disease_associations(query_symbol, genedb = oncoEnrichR::genedb)
  }

  for(c in names(oncoEnrichR::msigdb[['COLLECTION']])){
    for(subcat in names(oncoEnrichR::msigdb[['COLLECTION']][[c]])){
      if(c == "C5"){
        enr <- oncoEnrichR::get_go_enrichment(query_entrezgene,
                                  background_entrez = background_entrez,
                                  minGSSize = ge_report[['enrichment']][['settings']][['min_gs_size']],
                                  q_value_cutoff = ge_report[['enrichment']][['settings']][['q_value_cutoff']],
                                  p_value_cutoff = ge_report[['enrichment']][['settings']][['p_value_cutoff']],
                                  ontology = subcat,
                                  genedb = oncoEnrichR::genedb)
        if(!is.null(enr)){
          ge_report[['enrichment']][['results']][['go']] <-
            dplyr::bind_rows(ge_report[['enrichment']][['results']][['go']], enr) %>%
            dplyr::distinct()
        }
      }else{
        db = paste0("MSIGdb/",c,"/",subcat)
        enr <- oncoEnrichR::get_universal_enrichment(query_entrezgene,
                                         genedb = oncoEnrichR::genedb,
                                         background_entrez = background_entrez,
                                         minGSSize = ge_report[['enrichment']][['settings']][['min_gs_size']],
                                         q_value_cutoff = ge_report[['enrichment']][['settings']][['q_value_cutoff']],
                                         p_value_cutoff = ge_report[['enrichment']][['settings']][['p_value_cutoff']],
                                         TERM2GENE = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2GENE,
                                         TERM2NAME = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2NAME,
                                         TERM2SOURCE = oncoEnrichR::msigdb$TERM2SOURCE,
                                         dbsource = db)
        if(!is.null(enr)){
          ge_report[['enrichment']][['results']][['msigdb']] <-
            dplyr::bind_rows(ge_report[['enrichment']][['results']][['msigdb']], enr) %>%
            dplyr::distinct()
        }
      }
    }
  }


  db <- 'WikiPathways'
  ge_report[['enrichment']][['results']][['wikipathwaydb']] <-
    oncoEnrichR::get_universal_enrichment(query_entrezgene,
                              genedb = oncoEnrichR::genedb,
                              background_entrez = background_entrez,
                              minGSSize = ge_report[['enrichment']][['settings']][['min_gs_size']],
                              q_value_cutoff = ge_report[['enrichment']][['settings']][['q_value_cutoff']],
                              p_value_cutoff = ge_report[['enrichment']][['settings']][['p_value_cutoff']],
                              TERM2GENE = oncoEnrichR::wikipathwaydb$TERM2GENE,
                              TERM2NAME = oncoEnrichR::wikipathwaydb$TERM2NAME,
                              TERM2SOURCE = oncoEnrichR::wikipathwaydb$TERM2SOURCE,
                              dbsource = db)


  if(show_gtex_coexp == T){
    gtex_results <-
      oncoEnrichR::gtex_co_expression(query_symbol, genedb = oncoEnrichR::genedb, gids = gtex_atlasassay_groups)
    ge_report[['co_expression_gtex']][['plots']] <- gtex_results[['plots']]
    ge_report[['co_expression_gtex']][['df']] <- gtex_results[['df']]
  }

  if(show_ppi == T){
    ge_report[['ppi']][['stringdb']][['network']] <-
      oncoEnrichR::get_string_ppi_network(query_entrezgene,
                                          genedb = oncoEnrichR::genedb,
                                          settings = ge_report[['ppi']][['stringdb']][['settings']])
  }

  if(show_complex == T){
    ge_report[['protein_complex']] <-
      oncoEnrichR::annotate_protein_complex(query_symbol,
                                            genedb = oncoEnrichR::genedb,
                                            corum_db = oncoEnrichR::corum,
                                            uniprot_acc = oncoEnrichR::uniprot_acc)
  }

  if(show_tcga_aberration == T){
    for(v in c('cna_homdel','cna_ampl')){
      ge_report[['tcga_aberration']][['plot']][[v]] <-
        oncoEnrichR::tcga_aberration_plot(query_entrezgene, qsource = "entrezgene", genedb = oncoEnrichR::genedb, vtype = v)
      ge_report[['tcga_aberration']][['table']][[v]] <-
        oncoEnrichR::tcga_aberration_table(query_entrezgene, qsource = "entrezgene", genedb = oncoEnrichR::genedb, vtype = v)
    }
  }

  return(ge_report)

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

  assign("gene_enrich_report", report, envir = .GlobalEnv)
  #gene_enrich_report <- report

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
