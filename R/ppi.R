
get_network_hubs <- function(edges = NULL,
                             nodes = NULL,
                             genedb = NULL){

  stopifnot(!is.null(edges) & !is.null(nodes))
  validate_db_df(nodes, dbtype = "ppi_nodes")
  validate_db_df(edges, dbtype = "ppi_edges")
  validate_db_df(genedb, dbtype = "genedb")

  edges <- dplyr::select(edges, .data$preferredName_A,
                         .data$preferredName_B,
                         .data$from, .data$to, .data$weight)
  d <- igraph::graph_from_data_frame(d = edges, directed = F)

  ## hub score (Kleinberg"s hub centrality)
  hscore <- igraph::hub_score(d)
  hub_scores <- data.frame(symbol = names(sort(hscore$vector,decreasing = T)),
                           hub_score = round(sort(hscore$vector,decreasing = T), digits = 3),
                           stringsAsFactors = F) |>
    dplyr::left_join(dplyr::select(genedb, .data$symbol, .data$name),by=c("symbol")) |>
    dplyr::left_join(dplyr::select(nodes, .data$symbol, .data$query_node),by=c("symbol")) |>
    dplyr::filter(.data$query_node == T) |>
    dplyr::select(.data$symbol, .data$name, .data$hub_score) |>
    dplyr::distinct()

  #closeness_score <- igraph::closeness(d, mode="all")
  #degree.cent <- centr_degree(d, mode = "all")


  return(hub_scores)
}


get_network_communities <- function(edges = NULL, nodes = NULL){

  stopifnot(!is.null(edges) & !is.null(nodes))
  validate_db_df(nodes, dbtype = "ppi_nodes")
  validate_db_df(edges, dbtype = "ppi_edges")

  edges <- dplyr::select(edges,
                         .data$preferredName_A,
                         .data$preferredName_B,
                         .data$from,
                         .data$to, .data$weight)

  d <- igraph::graph_from_data_frame(d = edges, directed = F)


  ## communities, fast greedy modularity optimization algorithm for finding community structure,
  cties <- igraph::fastgreedy.community(d)
  edges_communities <- data.frame()
  if(length(cties) > 0){
    n <- 1
    while(n <= length(cties)){
      members_community_n <- cties[[n]]
      community <- dplyr::filter(
        edges, .data$preferredName_A %in% members_community_n &
          .data$preferredName_B %in% members_community_n) |>
        dplyr::mutate(community = n)

      edges_communities <- edges_communities |>
        dplyr::bind_rows(community)
      n <- n + 1
    }
  }

  nodes_communities <- data.frame(
    "id" = unique(c(unique(edges_communities$from),
                    unique(edges_communities$to))),
    stringsAsFactors = F) |>
    dplyr::inner_join(nodes, by = c("id")) |>
    dplyr::distinct()

  community_structure <- list()
  community_structure[["edges"]] <- edges_communities
  community_structure[["nodes"]] <- nodes_communities

  return(community_structure)
}


get_string_network_nodes_edges <-
  function(qgenes,
           all_query_nodes = NULL,
           settings = NULL,
           genedb = NULL){

  query_list <- paste(qgenes, collapse="%0d")

  edges <- jsonlite::fromJSON(
    paste0("https://string-db.org/api/json/",
           settings$query_type, "?species=9606&identifiers=",
           query_list, "&required_score=",
           settings$minimum_score, "&add_nodes=", settings$add_nodes))

  if(NROW(edges) > 0){

    edges <- edges |>
      dplyr::left_join(dplyr::select(genedb, .data$entrezgene, .data$symbol),
                       by = c("preferredName_A" = "symbol")) |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::rename(entrezgene_a = .data$entrezgene) |>
      #dplyr::mutate(entrezgene_a = as.character(.data$entrezgene_a)) |>
      dplyr::mutate(from = paste0("s", .data$entrezgene_a)) |>
      dplyr::left_join(dplyr::select(genedb, .data$entrezgene, .data$oncogene,
                                     .data$tumor_suppressor,
                                     .data$cancer_driver),by=c("entrezgene_a" = "entrezgene")) |>
      dplyr::rename(oncogene_A = .data$oncogene,
                    tsgene_A = .data$tumor_suppressor,
                    cdriver_A = .data$cancer_driver) |>
      dplyr::left_join(dplyr::select(genedb, .data$entrezgene, .data$symbol),
                       by = c("preferredName_B" = "symbol")) |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::rename(entrezgene_b = .data$entrezgene) |>
      #dplyr::mutate(entrezgene_b = as.character(.data$entrezgene_b)) |>
      dplyr::mutate(to = paste0("s", .data$entrezgene_b)) |>
      dplyr::left_join(dplyr::select(genedb, .data$entrezgene,
                                     .data$oncogene, .data$tumor_suppressor,
                                     .data$cancer_driver),
                       by = c("entrezgene_a" = "entrezgene")) |>
      dplyr::rename(oncogene_B = .data$oncogene,
                    tsgene_B = .data$tumor_suppressor,
                    cdriver_B = .data$cancer_driver) |>
      dplyr::mutate(interaction_symbol =
                      paste0(.data$preferredName_A, "_",
                             .data$preferredName_B)) |>
      dplyr::left_join(dplyr::select(all_query_nodes, .data$symbol, .data$query_node),
                       by = c("preferredName_A" = "symbol")) |>
      dplyr::rename(query_node_A = .data$query_node) |>
      dplyr::left_join(dplyr::select(all_query_nodes, .data$symbol, .data$query_node),
                       by = c("preferredName_B" = "symbol")) |>
      dplyr::rename(query_node_B = .data$query_node) |>
      dplyr::mutate(weight = .data$score) |>
      dplyr::mutate(title = paste0("Interaction score:", .data$score)) |>
      dplyr::distinct() |>
      dplyr::select(-c(.data$ncbiTaxonId, .data$stringId_A, .data$stringId_B))

    nodes <- data.frame(
      "symbol" = unique(c(edges$preferredName_A, edges$preferredName_B)),
      stringsAsFactors = F) |>
      dplyr::distinct() |>
      dplyr::left_join(dplyr::select(
        genedb, .data$symbol, .data$entrezgene,
        .data$tumor_suppressor, .data$oncogene,
        .data$cancer_driver, .data$genename,
        .data$targeted_cancer_drugs_ep,
        .data$targeted_cancer_drugs_lp),
        by = "symbol") |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::left_join(dplyr::select(all_query_nodes,
                                     .data$symbol,
                                     .data$query_node),
                       by = "symbol") |>
      dplyr::mutate(id = paste0("s", .data$entrezgene)) |>
      dplyr::distinct() |>
      dplyr::bind_rows(all_query_nodes) |>
      dplyr::distinct() |>
      dplyr::mutate(query_node = dplyr::if_else(
        is.na(.data$query_node), FALSE, as.logical(.data$query_node)))

    network <- list()
    network$edges <- edges
    network$nodes <- nodes

  }
  else{
    network <- NULL
  }

  return(network)


}

get_ppi_network <- function(qgenes,
                            ppi_source = "STRING",
                            genedb = NULL,
                            cancerdrugdb = NULL,
                            settings = NULL,
                            logger = NULL){

  stopifnot(!is.null(logger))
  stopifnot(!is.null(settings))
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(cancerdrugdb))
  stopifnot(!is.null(cancerdrugdb[['network']]))
  stopifnot("query_type" %in% names(settings))
  stopifnot(identical(names(settings),
                      c("minimum_score",
                        "visnetwork_shape",
                        "visnetwork_shadow",
                        "show_drugs",
                        "add_nodes",
                        "query_type"))
  )

  stopifnot(settings$query_type == "interaction_partners" |
              settings$query_type == "network")
  validate_db_df(genedb, dbtype = "genedb")
  stopifnot(is.integer(qgenes))

  query_nodes <- data.frame("entrezgene" = qgenes, stringsAsFactors = F) |>
    dplyr::distinct() |>

    dplyr::left_join(dplyr::select(
      genedb, .data$symbol, .data$entrezgene, .data$tumor_suppressor,
      .data$cancer_driver, .data$oncogene, .data$genename,
      .data$targeted_cancer_drugs_ep, .data$targeted_cancer_drugs_lp),
      by = c("entrezgene" = "entrezgene")) |>
    dplyr::mutate(id = paste0("s", .data$entrezgene)) |>
    dplyr::mutate(query_node = T) |>
    dplyr::distinct()

  ## Occasionally, some entrezgene identifiers are ambiguous (wrt. gene symbols)
  ## if this is the case, pick random entry from duplicate records, so that
  ## all entrezgene identifiers point to a unique row

  if(nrow(query_nodes) > length(qgenes)){


    duplicate_records <- as.data.frame(
      query_nodes |>
      dplyr::group_by(.data$entrezgene) |>
      dplyr::summarise(n = dplyr::n(),
                       .groups = "drop") |>
      dplyr::filter(.data$n > 1)
    )

    log4r_info(logger, paste0(
      "Resolving ambiguous Entrez gene identifiers: ",
      paste(duplicate_records$entrezgene, collapse=", ")))

    query_nodes_clean <-
      as.data.frame(
        query_nodes |>
          dplyr::group_by(.data$entrezgene) |>
          dplyr::summarise(n = dplyr::n(),
                           .groups = "drop") |>
          dplyr::filter(.data$n == 1)
      )

    for(i in 1:nrow(duplicate_records)){
      entrezgene_id = duplicate_records[i,]$entrezgene
      sample_row <- dplyr::sample_n(
        dplyr::filter(query_nodes, .data$entrezgene == entrezgene_id),
        1
      )
      query_nodes_clean <- query_nodes_clean |>
        dplyr::bind_rows(sample_row)

      i <- i + 1
    }

    if(length(qgenes) == nrow(query_nodes_clean)){
      query_nodes <- query_nodes_clean
    }

  }

  log4r_info(logger, "STRINGdb: retrieving protein-protein interaction network from (v11.5)")
  log4r_info(logger, paste0("STRINGdb: Settings -  required_score = ",
                                  settings$minimum_score,", add_nodes = ",settings$add_nodes))

  all_edges <- data.frame()
  all_nodes <- data.frame()

  if(length(qgenes) > 200){
    i <- 1
    omnipathdb <- data.frame()
    while(i <= length(qgenes)){
      qgenes_set <- qgenes[i:min(length(qgenes),i + 199)]

      ppi_network_data <- get_string_network_nodes_edges(
        qgenes = qgenes_set,
        all_query_nodes = query_nodes,
        settings = settings,
        genedb = genedb)

      all_edges <- all_edges |>
        dplyr::bind_rows(ppi_network_data$edges) |>
        dplyr::distinct()
      all_nodes <- all_nodes |>
        dplyr::bind_rows(ppi_network_data$nodes) |>
        dplyr::distinct()
      i <- i + 199
    }

  }else{
    ppi_network_data <- get_string_network_nodes_edges(
      qgenes = qgenes,
      all_query_nodes = query_nodes,
      settings = settings,
      genedb = genedb)

    if(!is.null(ppi_network_data)){

      all_edges <- all_edges |>
        dplyr::bind_rows(ppi_network_data$edges) |>
        dplyr::distinct()
      all_nodes <- all_nodes |>
        dplyr::bind_rows(ppi_network_data$nodes) |>
        dplyr::distinct()
    }
  }


  if(nrow(all_edges) > 0 & nrow(all_nodes) > 0){

    all_nodes <- all_nodes |>
      dplyr::mutate(shape = dplyr::if_else(
        .data$query_node == T,
        settings$visnetwork_shape,
        as.character("box")))
    all_nodes$shadow <- settings$visnetwork_shadow

    all_nodes$title <- unlist(
      lapply(stringr::str_match_all(
        all_nodes$genename,">.+<"), paste, collapse=","))
    all_nodes$title <- stringr::str_replace_all(all_nodes$title,">|<", "")
    all_nodes$label <- all_nodes$symbol # Node label

    all_edges$width <- all_edges$weight

    all_nodes$gene_category <- "protein_coding"
    all_nodes$size <- 25
    all_nodes <- all_nodes |>
      dplyr::mutate(color.background  = dplyr::if_else(
        .data$query_node == T, "lightblue", "mistyrose")) |>
      dplyr::mutate(color.background = dplyr::if_else(
        .data$tumor_suppressor == T & .data$oncogene == F,
        "firebrick", as.character(.data$color.background),
        as.character(.data$color.background))) |>
      dplyr::mutate(color.background = dplyr::if_else(
        .data$oncogene == T & .data$tumor_suppressor == F,
        "darkolivegreen", as.character(.data$color.background),
        as.character(.data$color.background))) |>
      dplyr::mutate(color.background = dplyr::if_else(
        .data$oncogene == T & .data$tumor_suppressor == T,
        "black", as.character(.data$color.background),
        as.character(.data$color.background))) |>
      dplyr::mutate(color.border = "black",
                    color.highlight.background = "orange",
                    color.highlight.border = "darkred",
                    font.color = "white") |>
      dplyr::mutate(font.color = dplyr::if_else(
        .data$query_node == T | .data$color.background == "mistyrose",
        "black", as.character(.data$font.color),
        as.character(.data$font.color)))

    nodes_with_drugs <- all_nodes
    edges_with_drugs <- all_edges

    if(settings$show_drugs == T){
      drug_target_ids <-
        dplyr::inner_join(dplyr::select(all_nodes, .data$id),
                          dplyr::select(cancerdrugdb[['network']]$edges, .data$from),
                          by = c("id" = "from")) |>
        dplyr::distinct()

      if(nrow(drug_target_ids) > 0){
        drug_edges <- dplyr::inner_join(
          cancerdrugdb[['network']]$edges,
          drug_target_ids, by = c("from" = "id"))
        drug_nodes <- drug_edges |>
          dplyr::select(-.data$from) |>
          dplyr::inner_join(cancerdrugdb[['network']]$nodes,
                            by = c("to" = "id")) |>
          dplyr::rename(id = .data$to) |>
          dplyr::distinct()

        nodes_with_drugs <- dplyr::bind_rows(all_nodes, drug_nodes)
        edges_with_drugs <- dplyr::bind_rows(all_edges, drug_edges)

      }
    }

    network <- list()
    network[["source"]] <- ppi_source
    network[["complete_network"]] <- list()
    network[["complete_network"]][["nodes"]] <- all_nodes
    network[["complete_network"]][["edges"]] <- all_edges
    if(settings$show_drugs == T){
      network[["complete_network"]][["nodes"]] <- nodes_with_drugs
      network[["complete_network"]][["edges"]] <- edges_with_drugs
    }
    network[["community_network"]] <-
      get_network_communities(edges = all_edges,
                              nodes = all_nodes)
    network[["hubscores"]] <-
      get_network_hubs(edges = all_edges,
                       nodes = all_nodes,
                       genedb = genedb)
  }else{
    network <- list()
    network[["source"]] <- ppi_source
    network[["complete_network"]] <- list()
    network[["complete_network"]][["nodes"]] <- data.frame()
    network[["complete_network"]][["edges"]] <- data.frame()
    network[["hubscores"]] <- data.frame()
    network[['community_network']] <- NULL
  }

  return(network)

}
