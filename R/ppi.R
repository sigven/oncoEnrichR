
get_network_hubs <- function(edges = NULL,
                             nodes = NULL,
                             genedb = NULL) {

  stopifnot(!is.null(edges) & !is.null(nodes))
  validate_db_df(nodes, dbtype = "ppi_nodes")
  validate_db_df(edges, dbtype = "ppi_edges")
  validate_db_df(genedb, dbtype = "genedb")

  edges <- edges |>
    dplyr::filter(
      !is.na(.data$symbol_A) &
        !is.na(.data$symbol_B)) |>
    dplyr::select(
      c("symbol_A", "symbol_B", "from",
        "to", "weight"))

  d <- igraph::graph_from_data_frame(d = edges, directed = F)

  ## hub score (Kleinberg"s hub centrality)
  hscore <- igraph::hub_score(d)
  hub_scores <- data.frame(
    symbol = names(sort(hscore$vector,decreasing = T)),
    hub_score = round(sort(hscore$vector,decreasing = T), digits = 3),
    stringsAsFactors = F) |>
    dplyr::left_join(
      dplyr::select(
        genedb, c("symbol", "name")),
      by = c("symbol"), relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::select(
        nodes, c("symbol", "query_node")),
      by = c("symbol"), relationship = "many-to-many") |>
    dplyr::filter(.data$query_node == T) |>
    dplyr::select(c("symbol", "name", "hub_score")) |>
    dplyr::distinct()

  #closeness_score <- igraph::closeness(d, mode="all")
  #degree.cent <- centr_degree(d, mode = "all")


  return(hub_scores)
}


get_network_communities <- function(edges = NULL, nodes = NULL) {

  stopifnot(!is.null(edges) & !is.null(nodes))
  validate_db_df(nodes, dbtype = "ppi_nodes")
  validate_db_df(edges, dbtype = "ppi_edges")

  edges <- edges |>
    dplyr::filter(
      !is.na(.data$symbol_A) &
        !is.na(.data$symbol_B)) |>
    dplyr::select(
      c("symbol_A","symbol_B",
        "from","to","weight"))

  d <- igraph::graph_from_data_frame(d = edges, directed = F)

  ## communities, fast greedy modularity optimization algorithm
  ## for finding community structure,
  cties <- igraph::fastgreedy.community(d)
  edges_communities <- data.frame()
  if (length(cties) > 0) {
    n <- 1
    while(n <= length(cties)) {
      members_community_n <- cties[[n]]
      community <- dplyr::filter(
        edges, .data$symbol_A %in% members_community_n &
          .data$symbol_B %in% members_community_n) |>
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
    dplyr::inner_join(nodes, by = c("id"), relationship = "many-to-many") |>
    dplyr::distinct()

  community_structure <- list()
  community_structure[["edges"]] <- edges_communities
  community_structure[["nodes"]] <- nodes_communities

  return(community_structure)
}


get_biogrid_network_nodes_edges <-
  function(all_query_nodes = NULL,
           add_nodes = 50,
           minimum_evidence_items = 2,
           show_isolated_nodes = FALSE,
           genedb = NULL,
           biogrid = NULL) {


    stopifnot(!is.null(all_query_nodes))
    stopifnot(is.data.frame(all_query_nodes))
    stopifnot(NROW(all_query_nodes) > 1)

    ppi_edges <- data.frame()

    ## match edges (interactions) in BioGRID entrezgene_B
    edges_part1 <- biogrid |>
      dplyr::inner_join(
        dplyr::select(all_query_nodes, c("entrezgene","symbol")),
        by = c("entrezgene_B" = "entrezgene"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(symbol_B = "symbol") |>
      dplyr::left_join(
        dplyr::select(
          genedb, c("entrezgene", "symbol")),
        by = c("entrezgene_A" = "entrezgene"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(symbol_A = "symbol") |>
      dplyr::arrange(
        .data$symbol_A,
        .data$symbol_B,
        dplyr::desc(.data$pmid)) |>
      dplyr::mutate(
        evidence = paste(
          .data$method,
          .data$throughput,
          .data$pmid, sep ="|"
        )
      )

    ## match edges (interactions) in BioGRID against entrezgene_B
    edges_part2 <- biogrid |>
      dplyr::inner_join(
        dplyr::select(all_query_nodes, c("entrezgene","symbol")),
        by = c("entrezgene_A" = "entrezgene"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(symbol_A = "symbol") |>
      dplyr::left_join(
        dplyr::select(
          genedb, c("entrezgene", "symbol")),
        by = c("entrezgene_B" = "entrezgene"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(symbol_B = "symbol") |>
      dplyr::arrange(.data$symbol_A,
                     .data$symbol_B,
                     dplyr::desc(.data$pmid)) |>
      dplyr::mutate(
        evidence = paste(
          .data$method,
          .data$throughput,
          .data$pmid, sep ="|"
        )
      )

    all_edges <- as.data.frame(
      dplyr::bind_rows(
        edges_part1,
        edges_part2) |>
        dplyr::distinct() |>
        dplyr::filter(.data$symbol_A != .data$symbol_B) |>
        dplyr::arrange(.data$symbol_A, .data$symbol_B) |>
        dplyr::group_by(
          .data$entrezgene_A,
          .data$entrezgene_B,
          .data$symbol_A,
          .data$symbol_B) |>
        dplyr::summarise(
          evidence = paste(
            unique(.data$evidence), collapse = ","
          ),
          num_evidence_items = dplyr::n(),
          .groups = "drop") |>

        ## only keep interactions with a minimum number
        ## of evidence items
        dplyr::filter(
          .data$num_evidence_items >= minimum_evidence_items) |>
        dplyr::mutate(
          querynode_A = dplyr::if_else(
            .data$entrezgene_A %in% all_query_nodes$entrezgene,
            TRUE,
            FALSE
          ),
          querynode_B = dplyr::if_else(
            .data$entrezgene_B %in% all_query_nodes$entrezgene,
            TRUE,
            FALSE
          )
        )
    )

    ## network - query nodes only
    ppi_edges <- all_edges |>
      dplyr::filter(
        .data$querynode_A == T &
          .data$querynode_B == T)

    if (add_nodes > 0) {

      ## attach non-query nodes with most interactions with query nodes
      edges_single_node_inqueryset <- all_edges |>
        dplyr::filter(.data$querynode_A == F |
                        .data$querynode_B == F)

        if (NROW(edges_single_node_inqueryset) > 0) {

          set1 <- edges_single_node_inqueryset |>
            dplyr::filter(.data$querynode_A == F)

          if (NROW(set1) > 0) {
            set1 <- set1 |>
              dplyr::select(c("entrezgene_A", "entrezgene_B")) |>
              dplyr::distinct() |>
              dplyr::group_by(.data$entrezgene_A) |>
              dplyr::summarise(num_interactions = dplyr::n()) |>
              dplyr::arrange(dplyr::desc(.data$num_interactions)) |>
              dplyr::rename(entrezgene = "entrezgene_A")
          }

          set2 <- edges_single_node_inqueryset |>
            dplyr::filter(.data$querynode_B == F)

          if (NROW(set2) > 0) {
            set2 <- set2 |>
            dplyr::select(c("entrezgene_A", "entrezgene_B")) |>
            dplyr::distinct() |>
            dplyr::group_by(.data$entrezgene_B) |>
            dplyr::summarise(num_interactions = dplyr::n()) |>
            dplyr::arrange(dplyr::desc(.data$num_interactions)) |>
            dplyr::rename(entrezgene = "entrezgene_B")
          }

          if (NROW(set1) > 0 | NROW(set2) > 0) {

            all_non_querynodes_ranked <-
              dplyr::bind_rows(
                set1, set2) |>
              dplyr::arrange(
                dplyr::desc(.data$num_interactions)) |>
              utils::head(add_nodes)

            edges_added <-
              edges_single_node_inqueryset |>
              dplyr::filter(
                .data$entrezgene_A %in% all_non_querynodes_ranked$entrezgene |
                  .data$entrezgene_B %in% all_non_querynodes_ranked$entrezgene
              )

            if (NROW(edges_added) > 0) {
              ppi_edges <- dplyr::bind_rows(
                ppi_edges, edges_added
              )
            }
          }
        }

    }

    network <- NULL

    ## ignore proto-oncogenes/tumor-suppressors with weak
    ## support/confidence in the protein-protein interaction network
    genedb <- genedb |>
      dplyr::mutate(oncogene = dplyr::if_else(
        .data$oncogene_confidence_level == "MODERATE",
        FALSE,
        as.logical(.data$oncogene)
      )) |>
      dplyr::mutate(tumor_suppressor = dplyr::if_else(
        .data$tsg_confidence_level == "MODERATE",
        FALSE,
        as.logical(.data$tumor_suppressor)
      ))


    if(NROW(ppi_edges) > 0){

      edges <- ppi_edges |>
        dplyr::mutate(from = paste0("s", .data$entrezgene_A)) |>
        dplyr::left_join(
          dplyr::select(
            genedb,
            c("entrezgene","oncogene",
              "tumor_suppressor",
              "cancer_driver")),
          by = c("entrezgene_A" = "entrezgene"), relationship = "many-to-many") |>
        dplyr::rename(oncogene_A = "oncogene",
                      tsgene_A = "tumor_suppressor",
                      cdriver_A = "cancer_driver") |>
        dplyr::mutate(to = paste0("s", .data$entrezgene_B)) |>
        dplyr::left_join(
          dplyr::select(genedb,
                        c("entrezgene", "oncogene",
                          "tumor_suppressor","cancer_driver")),
          by = c("entrezgene_B" = "entrezgene"), relationship = "many-to-many") |>
        dplyr::rename(oncogene_B = "oncogene",
                      tsgene_B = "tumor_suppressor",
                      cdriver_B = "cancer_driver") |>
        dplyr::mutate(interaction_symbol =
                        paste0(.data$symbol_A, "_",
                               .data$symbol_B)) |>
        dplyr::mutate(weight = dplyr::if_else(
          .data$num_evidence_items >= 5,
          as.numeric(2),
          as.numeric(1)
        )) |>
        #dplyr::mutate(weight = .data$num_evidence_items) |>
        dplyr::mutate(title = paste0("Number of evidence items:",
                                     .data$num_evidence_items)) |>
        dplyr::distinct()

      nodes <- data.frame(
        "entrezgene" = unique(c(edges$entrezgene_A,
                            edges$entrezgene_B)),
        stringsAsFactors = F) |>
        dplyr::distinct() |>
        dplyr::left_join(dplyr::select(
          genedb,
          c("symbol","entrezgene",
            "tumor_suppressor","oncogene",
            "cancer_driver","genename",
            "targeted_cancer_drugs_ep",
            "targeted_cancer_drugs_lp")),
          by = "entrezgene", relationship = "many-to-many") |>
        dplyr::filter(!is.na(.data$entrezgene)) |>
        dplyr::left_join(
          dplyr::select(all_query_nodes,
                        c("entrezgene","query_node")),
          by = "entrezgene", relationship = "many-to-many") |>
        dplyr::mutate(id = paste0("s", .data$entrezgene)) |>
        dplyr::distinct() |>
        dplyr::mutate(query_node = dplyr::if_else(
          is.na(.data$query_node),
          FALSE,
          as.logical(.data$query_node)))


      if(show_isolated_nodes == T){
        nodes <- nodes |>
          dplyr::bind_rows(all_query_nodes) |>
          dplyr::distinct() |>
          dplyr::mutate(query_node = dplyr::if_else(
            is.na(.data$query_node),
            FALSE,
            as.logical(.data$query_node)))
      }

      network <- list()
      network$edges <- edges
      network$nodes <- nodes

    }

    return(network)

  }


get_string_network_nodes_edges <-
  function(qgenes,
           all_query_nodes = NULL,
           query_type = "network",
           minimum_score = 0.9,
           add_nodes = 50,
           show_isolated_nodes = FALSE,
           network_type = "physical",
           genedb = NULL) {

  query_list <- paste(qgenes, collapse = "%0d")

  edges <- jsonlite::fromJSON(
    paste0("https://string-db.org/api/json/",
           query_type,
           "?species=9606&identifiers=", query_list,
           "&required_score=", as.integer(minimum_score * 1000),
           "&network_type=", network_type,
           "&add_nodes=", add_nodes))


  ## ignore proto-oncogenes/tumor-suppressors with weak
  ## support/confidence in the protein-protein interaction network
  genedb <- genedb |>
    dplyr::mutate(oncogene = dplyr::if_else(
      .data$oncogene_confidence_level == "MODERATE",
      FALSE,
      as.logical(.data$oncogene)
    )) |>
    dplyr::mutate(tumor_suppressor = dplyr::if_else(
      .data$tsg_confidence_level == "MODERATE",
      FALSE,
      as.logical(.data$tumor_suppressor)
    ))

  if (NROW(edges) > 0) {

    edges <- edges |>
      dplyr::rename(symbol_A = "preferredName_A",
                    symbol_B = "preferredName_B") |>
      dplyr::left_join(
        dplyr::select(genedb, c("entrezgene", "symbol")),
        by = c("symbol_A" = "symbol"), relationship = "many-to-many") |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::rename(entrezgene_A = "entrezgene") |>
      dplyr::mutate(from = paste0("s", .data$entrezgene_A)) |>
      dplyr::left_join(
        dplyr::select(
          genedb,
          c("entrezgene","oncogene",
            "tumor_suppressor",
            "cancer_driver")),
        by = c("entrezgene_A" = "entrezgene"), relationship = "many-to-many") |>
      dplyr::rename(oncogene_A = "oncogene",
                    tsgene_A = "tumor_suppressor",
                    cdriver_A = "cancer_driver") |>
      dplyr::left_join(dplyr::select(genedb, c("entrezgene", "symbol")),
                       by = c("symbol_B" = "symbol"), relationship = "many-to-many") |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::rename(entrezgene_B = "entrezgene") |>
      dplyr::mutate(to = paste0("s", .data$entrezgene_B)) |>
      dplyr::left_join(
        dplyr::select(genedb,
                      c("entrezgene", "oncogene",
                        "tumor_suppressor","cancer_driver")),
        by = c("entrezgene_B" = "entrezgene"), relationship = "many-to-many") |>
      dplyr::rename(oncogene_B = "oncogene",
                    tsgene_B = "tumor_suppressor",
                    cdriver_B = "cancer_driver") |>
      dplyr::mutate(interaction_symbol =
                      paste0(.data$symbol_A, "_",
                             .data$symbol_B)) |>
      dplyr::left_join(
        dplyr::select(all_query_nodes, c("symbol","query_node")),
        by = c("symbol_A" = "symbol"), relationship = "many-to-many") |>
      dplyr::rename(querynode_A = "query_node") |>
      dplyr::left_join(
        dplyr::select(all_query_nodes, c("symbol","query_node")),
        by = c("symbol_B" = "symbol"), relationship = "many-to-many") |>
      dplyr::rename(querynode_B = "query_node") |>
      dplyr::mutate(
        querynode_A = dplyr::if_else(
          is.na(.data$querynode_A),
          FALSE,
          as.logical(.data$querynode_A)
        ),
        querynode_B = dplyr::if_else(
          is.na(.data$querynode_B),
          FALSE,
          as.logical(.data$querynode_B)
        )
      ) |>
      dplyr::mutate(weight = .data$score) |>
      dplyr::mutate(title = paste0("Interaction score:", .data$score)) |>
      dplyr::distinct() |>
      dplyr::select(-c("ncbiTaxonId", "stringId_A", "stringId_B"))

    nodes <- data.frame(
      "symbol" = unique(c(edges$symbol_A, edges$symbol_B)),
      stringsAsFactors = F) |>
      dplyr::distinct() |>
      dplyr::left_join(dplyr::select(
        genedb,
        c("symbol","entrezgene",
          "tumor_suppressor","oncogene",
          "cancer_driver","genename",
          "targeted_cancer_drugs_ep",
          "targeted_cancer_drugs_lp")),
        by = "symbol", relationship = "many-to-many") |>
      dplyr::filter(!is.na(.data$entrezgene)) |>
      dplyr::left_join(
        dplyr::select(all_query_nodes,
                      c("symbol","query_node")),
        by = "symbol", relationship = "many-to-many") |>
      dplyr::mutate(id = paste0("s", .data$entrezgene)) |>
      dplyr::distinct() |>
      dplyr::mutate(query_node = dplyr::if_else(
        is.na(.data$query_node),
        FALSE,
        as.logical(.data$query_node)))


    if(show_isolated_nodes == T){
      nodes <- nodes |>
        dplyr::bind_rows(all_query_nodes) |>
        dplyr::distinct() |>
        dplyr::mutate(query_node = dplyr::if_else(
          is.na(.data$query_node),
          FALSE,
          as.logical(.data$query_node)))
    }

    network <- list()
    network$edges <- edges
    network$nodes <- nodes

  }
  else{
    network <- NULL
  }

  return(network)


}

get_ppi_network <- function(qgenes = NULL,
                            ppi_source = "string",
                            genedb = NULL,
                            cancerdrugdb = NULL,
                            biogrid = NULL,
                            settings = NULL,
                            ppi_source_release = NA) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(!is.null(qgenes))
  stopifnot(!is.null(settings))
  stopifnot(!is.null(genedb))
  validate_db_df(genedb, dbtype = "genedb")
  stopifnot(is.integer(qgenes))
  stopifnot(!is.null(cancerdrugdb))
  stopifnot(!is.null(biogrid))
  validate_db_df(biogrid, dbtype = "biogrid")
  stopifnot(!is.null(cancerdrugdb[['network']]))
  stopifnot(ppi_source == "string" | ppi_source == "biogrid")

  if (ppi_source == "string") {

    stopifnot(identical(
      names(settings[[ppi_source]]),
      c("minimum_score",
        "visnetwork_shape",
        "visnetwork_shadow",
        "show_drugs",
        "add_nodes",
        "query_type",
        "network_type",
        "show_isolated_nodes"))
    )
    stopifnot(settings[[ppi_source]]$query_type == "interaction_partners" |
                settings[[ppi_source]]$query_type == "network")
    stopifnot(settings[[ppi_source]]$network_type == "physical" |
                settings[[ppi_source]]$network_type == "functional")
  }


  if (ppi_source == "biogrid") {

    stopifnot(identical(
      names(settings[[ppi_source]]),
      c("minimum_evidence",
        "visnetwork_shape",
        "visnetwork_shadow",
        "show_drugs",
        "add_nodes",
        "show_isolated_nodes"))
    )
  }

  query_nodes <- data.frame("entrezgene" = qgenes, stringsAsFactors = F) |>
    dplyr::distinct() |>

    dplyr::left_join(dplyr::select(
      genedb,
      c("symbol","entrezgene","tumor_suppressor",
        "cancer_driver", "oncogene", "genename",
        "targeted_cancer_drugs_ep",
        "targeted_cancer_drugs_lp")),
      by = c("entrezgene" = "entrezgene"), relationship = "many-to-many") |>
    dplyr::mutate(id = paste0("s", .data$entrezgene)) |>
    dplyr::mutate(query_node = T) |>
    dplyr::distinct()

  ## Occasionally, some entrezgene identifiers are ambiguous (wrt. gene symbols)
  ## if this is the case, pick random entry from duplicate records, so that
  ## all entrezgene identifiers point to a unique row

  if (nrow(query_nodes) > length(qgenes)) {


    duplicate_records <- as.data.frame(
      query_nodes |>
      dplyr::group_by(.data$entrezgene) |>
      dplyr::summarise(n = dplyr::n(),
                       .groups = "drop") |>
      dplyr::filter(.data$n > 1)
    )

    lgr::lgr$info( paste0(
      "Resolving ambiguous Entrez gene identifiers: ",
      paste(duplicate_records$entrezgene, collapse = ", ")))

    query_nodes_clean <-
      as.data.frame(
        query_nodes |>
          dplyr::group_by(.data$entrezgene) |>
          dplyr::summarise(n = dplyr::n(),
                           .groups = "drop") |>
          dplyr::filter(.data$n == 1)
      )

    for (i in 1:nrow(duplicate_records)) {
      entrezgene_id = duplicate_records[i,]$entrezgene
      sample_row <- dplyr::sample_n(
        dplyr::filter(query_nodes, .data$entrezgene == entrezgene_id),
        1
      )
      query_nodes_clean <- query_nodes_clean |>
        dplyr::bind_rows(sample_row)

      i <- i + 1
    }

    if (length(qgenes) == nrow(query_nodes_clean)) {
      query_nodes <- query_nodes_clean
    }

  }


  ppi_source_verbose <- "BioGRID"
  if (ppi_source == "string") {

    ppi_source_verbose <- "STRING"
    lgr::lgr$info( paste0(
      ppi_source_verbose,
      ": retrieving protein-protein interaction network from API (",
      ppi_source_release,")"))
    lgr::lgr$info( paste0(
      ppi_source_verbose,
      ": Settings -  required minimum score = ",
      settings[[ppi_source]]$minimum_score,
      ", add_nodes = ",settings[[ppi_source]]$add_nodes,
      ", show_isolated_nodes = ",settings[[ppi_source]]$show_isolated_nodes,
      ", network_type = ", settings[[ppi_source]]$network_type))

    all_edges <- data.frame()
    all_nodes <- data.frame()

    if (length(qgenes) > 200) {
      i <- 1
      while(i <= length(qgenes)) {
        qgenes_set <- qgenes[i:min(length(qgenes),i + 199)]

        ppi_network_data <- get_string_network_nodes_edges(
          qgenes = qgenes_set,
          all_query_nodes = query_nodes,
          query_type = settings[[ppi_source]]$query_type,
          network_type = settings[[ppi_source]]$network_type,
          minimum_score = settings[[ppi_source]]$minimum_score,
          add_nodes = settings[[ppi_source]]$add_nodes,
          show_isolated_nodes = settings[[ppi_source]]$show_isolated_nodes,
          genedb = genedb)

        all_edges <- all_edges |>
          dplyr::bind_rows(ppi_network_data$edges) |>
          dplyr::distinct()
        all_nodes <- all_nodes |>
          dplyr::bind_rows(ppi_network_data$nodes) |>
          dplyr::distinct()
        i <- i + 199
      }

    } else {
      ppi_network_data <- get_string_network_nodes_edges(
        qgenes = qgenes,
        all_query_nodes = query_nodes,
        query_type = settings[[ppi_source]]$query_type,
        network_type = settings[[ppi_source]]$network_type,
        minimum_score = settings[[ppi_source]]$minimum_score,
        show_isolated_nodes = settings[[ppi_source]]$show_isolated_nodes,
        add_nodes = settings[[ppi_source]]$add_nodes,
        genedb = genedb)

      if (!is.null(ppi_network_data)) {

        all_edges <- all_edges |>
          dplyr::bind_rows(ppi_network_data$edges) |>
          dplyr::distinct()
        all_nodes <- all_nodes |>
          dplyr::bind_rows(ppi_network_data$nodes) |>
          dplyr::distinct()
      }
    }

  }

  if (ppi_source == "biogrid") {

    lgr::lgr$info(paste0(
      ppi_source_verbose,
      ": retrieving physical protein-protein interaction network (",
      ppi_source_release,")"))
    lgr::lgr$info( paste0(
      ppi_source_verbose,
      ": Settings -  minimum evidence = ",
      settings[[ppi_source]]$minimum_evidence,
      ", show_isolated_nodes = ",settings[[ppi_source]]$show_isolated_nodes,
      ", add_nodes = ",settings[[ppi_source]]$add_nodes))

    all_edges <- data.frame()
    all_nodes <- data.frame()

    ppi_network_data <- get_biogrid_network_nodes_edges(
      all_query_nodes = query_nodes,
      minimum_evidence_items = settings[[ppi_source]]$minimum_evidence,
      add_nodes = settings[[ppi_source]]$add_nodes,
      show_isolated_nodes = settings[[ppi_source]]$show_isolated_nodes,
      biogrid = biogrid,
      genedb = genedb)

    if (!is.null(ppi_network_data)) {

      all_edges <- all_edges |>
        dplyr::bind_rows(ppi_network_data$edges) |>
        dplyr::distinct()
      all_nodes <- all_nodes |>
        dplyr::bind_rows(ppi_network_data$nodes) |>
        dplyr::distinct()
    }
  }


  if (nrow(all_edges) > 0 & nrow(all_nodes) > 0) {

    all_nodes <- all_nodes |>
      dplyr::mutate(shape = dplyr::if_else(
        .data$query_node == T,
        settings[[ppi_source]]$visnetwork_shape,
        as.character("box")))
    all_nodes$shadow <- settings[[ppi_source]]$visnetwork_shadow

    all_nodes$title <- unlist(
      lapply(stringr::str_match_all(
        all_nodes$genename,">.+<"), paste, collapse = ","))
    all_nodes$title <- stringr::str_replace_all(
      all_nodes$title,">|<", "")
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

    if (settings[[ppi_source]]$show_drugs == T) {
      drug_target_ids <-
        dplyr::inner_join(
          dplyr::select(all_nodes, c("id")),
          dplyr::select(cancerdrugdb[['network']]$edges, c("from")),
          by = c("id" = "from"), relationship = "many-to-many") |>
        dplyr::distinct()

      if (nrow(drug_target_ids) > 0) {
        drug_edges <- dplyr::inner_join(
          cancerdrugdb[['network']]$edges,
          drug_target_ids, by = c("from" = "id"))
        drug_nodes <- drug_edges |>
          dplyr::select(-c("from")) |>
          dplyr::inner_join(cancerdrugdb[['network']]$nodes,
                            by = c("to" = "id"), relationship = "many-to-many") |>
          dplyr::rename(id = "to") |>
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
    if (settings[[ppi_source]]$show_drugs == T) {
      network[["complete_network"]][["nodes"]] <- nodes_with_drugs
      network[["complete_network"]][["edges"]] <- edges_with_drugs
    }
    lgr::lgr$info(paste0(
      ppi_source_verbose,
      ": retrieving network communities (",
      ppi_source_release,")"))
    network[["community_network"]] <-
      get_network_communities(edges = all_edges,
                              nodes = all_nodes)

    lgr::lgr$info(paste0(
      ppi_source_verbose,
      ": retrieving hub scores (",
      ppi_source_release,")"))
    network[["hubscores"]] <-
      get_network_hubs(edges = all_edges,
                       nodes = all_nodes,
                       genedb = genedb)
  } else {
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
