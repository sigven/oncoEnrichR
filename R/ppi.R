
#' Annotate geneset with protein complexes they are involved with
#'
#' @param qgenes character vector of human gene symbols
#' @param genedb data frame with gene annotations
#' @param settings oncoEnrichR report settings
#' @return data frame with unique protein complexes and their associated members
#'
get_string_ppi_network <- function(qgenes, genedb = NULL, settings = NULL){

  stopifnot(!is.null(settings))
  stopifnot(!is.null(genedb))
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")
  for(var in c('symbol','entrezgene','p_oncogene','tsgene','cdriver','tcga_driver','ensembl_gene_id','name',
             'gencode_gene_biotype','corum_id','ot_tractability_compound','signaling_pw','genename')){
    stopifnot(var %in% colnames(genedb))
  }
  network <- list()
  network$nodes <- NULL
  network$edges <- NULL
  if(settings$query_type != 'interaction_partners' & settings$query_type != 'network'){
    return(network)
  }
  query_list <- paste(qgenes, collapse="%0d")
  rlogging::message("STRINGdb: retrieving protein-protein interaction network from (v11)")
  rlogging::message(paste0("STRINGdb: Settings -  required_score = ",settings$minimum_score,", add_nodes = ",settings$add_nodes))

  query_nodes <- data.frame('entrezgene' = qgenes, stringsAsFactors = F) %>%
    dplyr::distinct() %>%
    dplyr::left_join(genedb, by=c("entrezgene" = "entrezgene")) %>%
    dplyr::mutate(id = paste0("s",entrezgene)) %>%
    dplyr::select(-c(corum_id,ensembl_gene_id)) %>%
    dplyr::mutate(query_node = T) %>%
    dplyr::distinct()

  all_links <- jsonlite::fromJSON(paste0('https://string-db.org/api/json/',settings$query_type,'?identifiers=',query_list,'&required_score=',settings$minimum_score,'&add_nodes=',settings$add_nodes)) %>%
    dplyr::left_join(dplyr::select(genedb,entrezgene,symbol),by=c("preferredName_A" = "symbol")) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::rename(entrezgene_a = entrezgene) %>%
    dplyr::mutate(entrezgene_a = as.character(entrezgene_a)) %>%
    dplyr::mutate(from = paste0("s",entrezgene_a)) %>%
    dplyr::left_join(dplyr::select(genedb, entrezgene, p_oncogene, tsgene, cdriver, tcga_driver),by=c("entrezgene_a" = "entrezgene")) %>%
    dplyr::rename(p_oncogene_A = p_oncogene, tsgene_A = tsgene, cdriver_A = cdriver, tcga_driver_A = tcga_driver) %>%
    dplyr::left_join(dplyr::select(genedb,entrezgene,symbol),by=c("preferredName_B" = "symbol")) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::rename(entrezgene_b = entrezgene) %>%
    dplyr::mutate(entrezgene_b = as.character(entrezgene_b)) %>%
    dplyr::mutate(to = paste0("s",entrezgene_b)) %>%
    dplyr::left_join(dplyr::select(genedb, entrezgene, p_oncogene, tsgene, cdriver, tcga_driver),by=c("entrezgene_a" = "entrezgene")) %>%
    dplyr::rename(p_oncogene_B = p_oncogene, tsgene_B = tsgene, cdriver_B = cdriver, tcga_driver_B = tcga_driver) %>%
    dplyr::mutate(interaction_symbol = paste0(preferredName_A,"_",preferredName_B)) %>%
    dplyr::left_join(dplyr::select(query_nodes, symbol, query_node), by=c("preferredName_A" = "symbol")) %>%
    dplyr::rename(query_node_A = query_node) %>%
    dplyr::left_join(dplyr::select(query_nodes, symbol, query_node), by=c("preferredName_B" = "symbol")) %>%
    dplyr::rename(query_node_B = query_node) %>%
    dplyr::mutate(weight = score) %>%
    dplyr::distinct()

  network_nodes <- data.frame('symbol' = c(all_links$preferredName_A, all_links$preferredName_B), stringsAsFactors = F) %>%
    dplyr::distinct() %>%
    dplyr::left_join(genedb, by=c("symbol" = "symbol")) %>%
    dplyr::select(-c(corum_id,ensembl_gene_id)) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::left_join(dplyr::select(query_nodes, symbol, query_node),by=c("symbol")) %>%
    dplyr::mutate(id = paste0("s",entrezgene)) %>%
    dplyr::distinct()

  all_nodes <- dplyr::bind_rows(query_nodes, network_nodes) %>%
    dplyr::distinct() %>%
    dplyr::mutate(query_node = dplyr::if_else(is.na(query_node),FALSE,as.logical(query_node)))


  all_nodes <- all_nodes %>% dplyr::mutate(shape = dplyr::if_else(query_node == T,settings$visnetwork_shape,as.character("box")))
  #all_nodes$shape  <- shape
  all_nodes$shadow <- settings$visnetwork_shadow

  all_nodes$title <- unlist(lapply(stringr::str_match_all(all_nodes$genename,">.+<"),paste,collapse=","))
  all_nodes$title <- stringr::str_replace_all(all_nodes$title,">|<", "")
  all_nodes$label  <- all_nodes$symbol # Node label

  all_links$width <- all_links$weight

  all_nodes$gene_category <- 'protein_coding'
  all_nodes$size <- 25
  all_nodes <- all_nodes %>%
    dplyr::mutate(color.background  = dplyr::if_else(query_node == T,"lightblue","mistyrose")) %>%
    dplyr::mutate(color.background = dplyr::if_else(tsgene == T & p_oncogene == F,"firebrick",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.background = dplyr::if_else(p_oncogene == T & tsgene == F,"darkolivegreen",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.background = dplyr::if_else(p_oncogene == T & tsgene == T,"black",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.border = 'black', color.highlight.background = 'orange', color.highlight.border = 'darkred', font.color = 'white') %>%
    dplyr::mutate(font.color = dplyr::if_else(query_node == T,"black",as.character(font.color),as.character(font.color)))

  network <- list()
  network$nodes <- all_nodes
  network$edges <- all_links
  return(network)

}
