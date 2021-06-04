get_go_enrichment <- function(query_entrez,
                              background_entrez = NULL,
                              ontology = "MF",
                              genedb = NULL,
                              p_value_cutoff = 0.05,
                              q_value_cutoff = 0.2,
                              p_value_adjustment_method = "BH",
                              min_geneset_size = 10,
                              max_geneset_size = 500,
                              simplify = F,
                              pool = F){

  oncoEnrichR:::log4r_info(
    paste0("Enrichment - GO: performing gene enrichment analysis of target set (subontology ",ontology,")"))
  oncoEnrichR:::log4r_info(paste0("Enrichment - GO: settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  oncoEnrichR:::log4r_info(paste0("Enrichment - GO: settings: p_value_adjustment_method = ",p_value_adjustment_method))
  oncoEnrichR:::log4r_info(paste0("Enrichment - GO: settings: minGSSize = ",min_geneset_size))
  oncoEnrichR:::log4r_info(paste0("Enrichment - GO: settings: maxGSSize = ",max_geneset_size))


  stopifnot(is.character(query_entrez))
  stopifnot(ontology == "MF" | ontology == "BP" |
              ontology == "CC" | ontology == "ALL")
  stopifnot(!is.null(genedb))
  stopifnot("symbol" %in% colnames(genedb) &
              "entrezgene" %in% colnames(genedb))
  if(is.null(background_entrez)){
    bg <- dplyr::select(genedb,entrezgene) %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::distinct()
    background_entrez <- bg$entrezgene
  }else{
    stopifnot(is.character(background_entrez))
  }
  ego <- clusterProfiler::enrichGO(gene          = query_entrez,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = ontology,
                                   minGSSize     = min_geneset_size,
                                   maxGSSize     = max_geneset_size,
                                   pAdjustMethod = p_value_adjustment_method,
                                   pvalueCutoff  = p_value_cutoff,
                                   qvalueCutoff  = q_value_cutoff,
                                   universe      = background_entrez,
                                   readable      = F,
                                   pool = pool)

  if(simplify == T){
    ego <- clusterProfiler::simplify(ego, cutoff=0.8,
                                     by="p.adjust", select_fun=min)
  }

  df <- as.data.frame(head(ego,5000))
  rownames(df) <- NULL
  if(ontology == "ALL" & "ONTOLOGY" %in% colnames(df)){
    df <- df %>% dplyr::rename(db = ONTOLOGY)
  }else{
    df <- df %>% dplyr::mutate(db = ontology)
  }
  if(nrow(df) > 0){
    df <- suppressWarnings(df %>%
      dplyr::mutate(db = paste0("GO_",db)) %>%
      dplyr::rename(go_id = ID, go_description = Description, count = Count,
                    gene_ratio = GeneRatio, background_ratio = BgRatio,
                    gene_id = geneID) %>%
      dplyr::mutate(
        go_description_link =
          paste0('<a href=\'http://amigo.geneontology.org/amigo/term/',
                 go_id,'\' target=\'_blank\'>',go_description,'</a>')) %>%
        tidyr::separate(gene_ratio,c('num_query_hits','num_query_all'),
                        sep='/',remove = F, convert = T) %>%
      dplyr::mutate(qvalue = as.numeric(qvalue)) %>%
      dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
      dplyr::mutate(
        qvalue =
          dplyr::if_else(!is.na(qvalue),
                         as.numeric(formatC(qvalue, format = "e",
                                            digits = 1)),
                         as.numeric(NA))) %>%
        dplyr::mutate(
          pvalue =
            dplyr::if_else(!is.na(pvalue),
                           as.numeric(formatC(pvalue, format = "e",
                                              digits = 1)),
                           as.numeric(NA))) %>%
        tidyr::separate(background_ratio,
                        c('num_background_hits','num_background_all'),
                        sep='/',remove = F, convert = T) %>%
      dplyr::mutate(
        enrichment_factor =
          round(as.numeric((num_query_hits/num_query_all) /
                             (num_background_hits/num_background_all)),
                digits = 1)) %>%
        dplyr::select(-c(num_query_hits,num_query_all,
                         num_background_hits,num_background_all))
    )

    gene2id <- NULL
    if(!is.null(genedb)){
      gene2id <- as.data.frame(
        df %>%
          dplyr::select(go_id,gene_id) %>%
          tidyr::separate_rows(gene_id,sep="/") %>%
          dplyr::mutate(
            gene_id = as.character(gene_id)) %>%
          dplyr::left_join(
            dplyr::select(genedb,entrezgene,symbol),
            by=c("gene_id" = "entrezgene")) %>%
          dplyr::mutate(
            entrez_url =
              paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                     gene_id,"' target=\'blank_\'>",symbol,"</a>")) %>%
          dplyr::group_by(go_id) %>%
          dplyr::summarise(
            gene_symbol_link =
              paste(unique(entrez_url),collapse=", "),
            gene_symbol = paste(unique(symbol),collapse=", "))
      )
    }
    if(!is.null(gene2id)){
      df <- df %>%
        dplyr::left_join(
          dplyr::select(gene2id, go_id,
                        gene_symbol_link,gene_symbol), by="go_id") %>%
        dplyr::rename(exact_source = go_id, description = go_description,
                      description_link = go_description_link) %>%
        dplyr::mutate(standard_name = exact_source)
    }

    if(!is.null(df)){
      df$setting_p_value_cutoff <- p_value_cutoff
      df$setting_q_value_cutoff <- q_value_cutoff
      df$setting_p_value_adj_method <- p_value_adjustment_method
      df$setting_min_geneset_size <- min_geneset_size
      df$setting_max_geneset_size <- max_geneset_size
    }
  }else{
    df <- NULL
  }

  return(df)

}


get_universal_enrichment <- function(query_entrez,
                                     background_entrez = NULL,
                                     genedb = NULL,
                                     p_value_cutoff = 0.05,
                                     p_value_adjustment_method = "BH",
                                     min_geneset_size = 10,
                                     max_geneset_size = 500,
                                     q_value_cutoff = 0.2,
                                     TERM2GENE = NULL,
                                     TERM2NAME = NULL,
                                     TERM2SOURCE = NULL,
                                     dbsource = ""){

  oncoEnrichR:::log4r_info(paste0("Enrichment - ",dbsource,": performing gene enrichment analysis of target set"))
  oncoEnrichR:::log4r_info(paste0("Enrichment - ",dbsource,": settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  oncoEnrichR:::log4r_info(paste0("Enrichment - ",dbsource,": settings: p_value_adjustment_method = ",p_value_adjustment_method))
  oncoEnrichR:::log4r_info(paste0("Enrichment - ",dbsource,": settings: minGSSize = ",min_geneset_size))
  oncoEnrichR:::log4r_info(paste0("Enrichment - ",dbsource,": settings: maxGSSize = ",max_geneset_size))

  stopifnot(is.character(query_entrez))
  stopifnot(!is.null(genedb) | !is.data.frame(genedb))
  #stopifnot(!is.null(TERM2SOURCE) | !is.data.frame(TERM2SOURCE))
  #stopifnot("standard_name" %in% colnames(TERM2SOURCE))
  stopifnot("symbol" %in% colnames(genedb) & "entrezgene" %in% colnames(genedb))
  if(is.null(background_entrez)){
    bg <- dplyr::select(genedb,entrezgene) %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::distinct()
    background_entrez <- bg$entrezgene
  }else{
    stopifnot(is.character(background_entrez))

  }
  enr <- clusterProfiler::enricher(gene          = query_entrez,
                                   universe      = background_entrez,
                                   pAdjustMethod = p_value_adjustment_method,
                                   minGSSize     = min_geneset_size,
                                   maxGSSize     = max_geneset_size,
                                   pvalueCutoff  = p_value_cutoff,
                                   qvalueCutoff  = q_value_cutoff,
                                   TERM2GENE     = TERM2GENE,
                                   TERM2NAME     = TERM2NAME)


  df <- as.data.frame(head(enr,5000))
  rownames(df) <- NULL
  if(nrow(df) > 0){
    df <- suppressWarnings(df %>%
      dplyr::rename(standard_name = ID,
                    description = Description,
                    count = Count,
                    gene_ratio = GeneRatio,
                    background_ratio = BgRatio,
                    gene_id = geneID)
    )

    if (dbsource == "WikiPathways"){
      df <- df %>%
        dplyr::mutate(description_link = paste0(
          "<a href=\"https://www.wikipathways.org/index.php/Pathway:",
          standard_name,"\" target='_blank'>",
          description,"</a>"
        )) %>%
        dplyr::mutate(exact_source = "https://wikipathways.org",
                      db = dbsource)
    }
    else if(dbsource == "KEGG"){
      df <- df %>%
        dplyr::mutate(description_link = paste0(
          "<a href=\"https://www.genome.jp/kegg-bin/show_pathway?",
          stringr::str_replace(standard_name,"hsa","map"),"\" target='_blank'>",
          description,"</a>"
        )) %>%
        dplyr::mutate(exact_source = "https://www.genome.jp/kegg/pathway.html",
                      db = dbsource)
    }
    else if(dbsource == "NetPath"){
      df <- df %>%
        dplyr::mutate(description_link = paste0(
          "<a href=\"http://netpath.org/pathways?path_id=",
          standard_name,"\" target='_blank'>",
          description,"</a>"
        )) %>%
        dplyr::mutate(exact_source = "http://netpath.org",
                      db = dbsource)
    }
    else{
      stopifnot(!is.null(TERM2SOURCE) | !is.data.frame(TERM2SOURCE))
      stopifnot("standard_name" %in% colnames(TERM2SOURCE))
      df <- df %>%
        dplyr::left_join(TERM2SOURCE, by="standard_name") %>%
        dplyr::mutate(description_link =
                        dplyr::if_else(!is.na(external_url),
                                       paste0("<a href='",external_url,
                                              "' target='_blank'>",
                                              description,"</a>"),
                                       description))


    }
    df <- suppressWarnings(df %>%
      dplyr::mutate(qvalue = as.numeric(qvalue)) %>%
      dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
      dplyr::mutate(
        qvalue =
          dplyr::if_else(!is.na(qvalue),
                         as.numeric(formatC(qvalue, format = "e",
                                            digits = 1)),
                         as.numeric(NA))) %>%
      dplyr::mutate(
        pvalue =
          dplyr::if_else(!is.na(pvalue),
                         as.numeric(formatC(pvalue, format = "e",
                                            digits = 1)),
                         as.numeric(NA))) %>%
      tidyr::separate(
        gene_ratio,c('num_query_hits','num_query_all'),
        sep='/',remove = F, convert = T) %>%
      tidyr::separate(
        background_ratio,c('num_background_hits','num_background_all'),
        sep='/',remove = F, convert = T) %>%
      dplyr::mutate(
        enrichment_factor =
          round(as.numeric((
            num_query_hits/num_query_all) /
              (num_background_hits/num_background_all)),digits = 1)) %>%
      dplyr::select(-c(num_query_hits,num_query_all,
                       num_background_hits,num_background_all)) %>%
      dplyr::mutate(
        db = dplyr::if_else(is.na(db) &
                              nchar(dbsource) > 0,
                            dbsource,as.character(db)))
    )

    gene2id <- NULL
    if(!is.null(genedb)){
      gene2id <- as.data.frame(
        df %>%
          dplyr::select(standard_name,gene_id) %>%
          tidyr::separate_rows(gene_id,sep="/") %>%
          dplyr::mutate(gene_id = as.character(gene_id)) %>%
          dplyr::left_join(
            dplyr::select(genedb,entrezgene,symbol),
            by=c("gene_id" = "entrezgene")) %>%
          dplyr::mutate(
            entrez_url =
              paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                     gene_id,"' target=\'blank_\'>",symbol,"</a>")) %>%
          dplyr::group_by(standard_name) %>%
          dplyr::summarise(
            gene_symbol_link =
              paste(unique(entrez_url),collapse=", "),
            gene_symbol = paste(unique(symbol),collapse=", "))
      )
    }
    if(!is.null(gene2id)){
      df <- df %>%
        dplyr::left_join(dplyr::select(gene2id,standard_name,
                                       gene_symbol_link,gene_symbol),
                         by="standard_name")
    }
    if(!is.null(df)){
      df$setting_p_value_cutoff <- p_value_cutoff
      df$setting_q_value_cutoff <- q_value_cutoff
      df$setting_p_value_adj_method <- p_value_adjustment_method
      df$setting_min_geneset_size <- min_geneset_size
      df$setting_max_geneset_size <- max_geneset_size
    }
  }else{
    df <- NULL
  }
  return(df)
}
