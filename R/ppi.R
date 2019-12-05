
get_network_hubs <- function(edges = NULL, nodes = NULL, genedb = NULL){

  stopifnot(!is.null(edges) & !is.null(nodes))
  oncoEnrichR::validate_db_df(nodes, dbtype = "ppi_nodes")
  oncoEnrichR::validate_db_df(edges, dbtype = "ppi_edges")

  edges <- dplyr::select(edges, preferredName_A, preferredName_B, from, to, weight)
  d <- igraph::graph_from_data_frame(d = edges, directed = F)

  ## hub score (Kleinberg's hub centrality)
  hscore <- igraph::hub_score(d)
  hub_scores <- data.frame(symbol = names(sort(hscore$vector,decreasing = T)),
                           hub_score = round(sort(hscore$vector,decreasing = T), digits = 3),
                           stringsAsFactors = F) %>%
    dplyr::left_join(dplyr::select(genedb, symbol,name),by=c("symbol")) %>%
    dplyr::left_join(dplyr::select(nodes, symbol, query_node),by=c("symbol")) %>%
    dplyr::filter(query_node == T) %>%
    dplyr::select(symbol, name, hub_score) %>%
    dplyr::distinct()

  #closeness_score <- igraph::closeness(d, mode="all")
  #degree.cent <- centr_degree(d, mode = "all")


  return(hub_scores)
}


get_network_communities <- function(edges = NULL, nodes = NULL){

  stopifnot(!is.null(edges) & !is.null(nodes))
  oncoEnrichR::validate_db_df(nodes, dbtype = "ppi_nodes")
  oncoEnrichR::validate_db_df(edges, dbtype = "ppi_edges")

  edges <- dplyr::select(edges, preferredName_A, preferredName_B, from, to, weight)

  d <- igraph::graph_from_data_frame(d = edges, directed = F)


  ## communities, fast greedy modularity optimization algorithm for finding community structure,
  cties <- igraph::fastgreedy.community(d)
  edges_communities <- data.frame()
  if(length(cties) > 0){
    n <- 1
    while(n <= length(cties)){
      members_community_n <- cties[[n]]
      community <- dplyr::filter(edges, preferredName_A %in% members_community_n & preferredName_B %in% members_community_n) %>%
        dplyr::mutate(community = n)

      edges_communities <- edges_communities %>%
        dplyr::bind_rows(community)
      n <- n + 1
    }
  }

  nodes_communities <- data.frame('id' = unique(c(unique(edges_communities$from),unique(edges_communities$to))), stringsAsFactors = F) %>%
    dplyr::inner_join(nodes, by=c("id")) %>%
    dplyr::distinct()

  community_structure <- list()
  community_structure[['edges']] <- edges_communities
  community_structure[['nodes']] <- nodes_communities

  return(community_structure)
}

get_ppi_network <- function(qgenes, ppi_source = "STRING", genedb = NULL, cancerdrugdb = NULL, settings = NULL){

  stopifnot(!is.null(settings))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(cancerdrugdb))
  stopifnot(settings$query_type == 'interaction_partners' | settings$query_type == 'network')
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")

  query_list <- paste(qgenes, collapse="%0d")

  query_nodes <- data.frame('entrezgene' = qgenes, stringsAsFactors = F) %>%
    dplyr::distinct() %>%
    dplyr::left_join(genedb, by=c("entrezgene" = "entrezgene")) %>%
    dplyr::mutate(id = paste0("s",entrezgene)) %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(query_node = T) %>%
    dplyr::distinct()

  rlogging::message("STRINGdb: retrieving protein-protein interaction network from (v11)")
  rlogging::message(paste0("STRINGdb: Settings -  required_score = ",settings$minimum_score,", add_nodes = ",settings$add_nodes))

  all_edges <- jsonlite::fromJSON(paste0('https://string-db.org/api/json/',settings$query_type,'?species=9606&identifiers=',query_list,'&required_score=',settings$minimum_score,'&add_nodes=',settings$add_nodes)) %>%
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
    dplyr::distinct() %>%
    dplyr::select(-c(ncbiTaxonId,stringId_A,stringId_B))

  network_nodes <- data.frame('symbol' = unique(c(all_edges$preferredName_A, all_edges$preferredName_B)), stringsAsFactors = F) %>%
    dplyr::distinct() %>%
    dplyr::left_join(genedb, by=c("symbol" = "symbol")) %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::left_join(dplyr::select(query_nodes, symbol, query_node),by=c("symbol")) %>%
    dplyr::mutate(id = paste0("s",entrezgene)) %>%
    dplyr::distinct()

  all_nodes <- dplyr::bind_rows(query_nodes, network_nodes) %>%
    dplyr::distinct() %>%
    dplyr::mutate(query_node = dplyr::if_else(is.na(query_node),FALSE,as.logical(query_node)))


  all_nodes <- all_nodes %>%
    dplyr::mutate(shape = dplyr::if_else(query_node == T,settings$visnetwork_shape,as.character("box")))
  all_nodes$shadow <- settings$visnetwork_shadow

  all_nodes$title <- unlist(lapply(stringr::str_match_all(all_nodes$genename,">.+<"),paste,collapse=","))
  all_nodes$title <- stringr::str_replace_all(all_nodes$title,">|<", "")
  all_nodes$label  <- all_nodes$symbol # Node label

  all_edges$width <- all_edges$weight

  all_nodes$gene_category <- 'protein_coding'
  all_nodes$size <- 25
  all_nodes <- all_nodes %>%
    dplyr::mutate(color.background  = dplyr::if_else(query_node == T,"lightblue","mistyrose")) %>%
    dplyr::mutate(color.background = dplyr::if_else(tsgene == T & p_oncogene == F,"firebrick",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.background = dplyr::if_else(p_oncogene == T & tsgene == F,"darkolivegreen",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.background = dplyr::if_else(p_oncogene == T & tsgene == T,"black",as.character(color.background),as.character(color.background))) %>%
    dplyr::mutate(color.border = 'black', color.highlight.background = 'orange', color.highlight.border = 'darkred', font.color = 'white') %>%
    dplyr::mutate(font.color = dplyr::if_else(query_node == T | color.background == "mistyrose","black",as.character(font.color),as.character(font.color)))

  nodes_with_drugs <- all_nodes
  edges_with_drugs <- all_edges

  if(settings$show_drugs == T){
    drug_target_ids <- dplyr::inner_join(dplyr::select(all_nodes, id), dplyr::select(cancerdrugdb$edges, from), by = c("id" = "from")) %>%
      dplyr::distinct()

    if(nrow(drug_target_ids) > 0){
      drug_edges <- dplyr::inner_join(cancerdrugdb$edges, drug_target_ids, by = c("from" = "id"))
      drug_nodes <- drug_edges %>%
        dplyr::select(-from) %>%
        dplyr::inner_join(cancerdrugdb$nodes, by = c("to" = "id")) %>%
        dplyr::rename(id = to) %>%
        dplyr::distinct()

      nodes_with_drugs <- dplyr::bind_rows(all_nodes, drug_nodes)
      edges_with_drugs <- dplyr::bind_rows(all_edges, drug_edges)

    }
  }

  network <- list()
  network[['source']] <- ppi_source
  network[['complete_network']] <- list()
  network[['complete_network']][['nodes']] <- all_nodes
  network[['complete_network']][['edges']] <- all_edges
  if(settings$show_drugs == T){
    network[['complete_network']][['nodes']] <- nodes_with_drugs
    network[['complete_network']][['edges']] <- edges_with_drugs
  }
  network[['community_network']] <- oncoEnrichR::get_network_communities(edges = all_edges, nodes = all_nodes)
  network[['hubscores']] <- oncoEnrichR::get_network_hubs(edges = all_edges, nodes = all_nodes, genedb = genedb)

  return(network)

}
