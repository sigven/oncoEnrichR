get_go_enrichment <- function(query_entrez,
                              background_entrez = NULL,
                              ontology = "MF",
                              genedb = NULL,
                              p_value_cutoff = 0.05,
                              q_value_cutoff = 0.2,
                              p_value_adjustment_method = "BH",
                              minGSSize = 10,
                              pool = F){

  rlogging::message(paste0("Enrichment - GO: performing gene enrichment analysis of query set (subontology ",ontology,")"))
  rlogging::message(paste0("Enrichment - GO: Settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  rlogging::message(paste0("Enrichment - GO: Settings: minGSSize = ",minGSSize))

  stopifnot(is.character(query_entrez))
  stopifnot(ontology == "MF" | ontology == "BP" | ontology == "CC" | ontology == "ALL")
  stopifnot(!is.null(genedb))
  stopifnot("symbol" %in% colnames(genedb) & "entrezgene" %in% colnames(genedb))
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
                                   minGSSize     = minGSSize,
                                   pAdjustMethod = p_value_adjustment_method,
                                   pvalueCutoff  = p_value_cutoff,
                                   qvalueCutoff  = q_value_cutoff,
                                   universe      = background_entrez,
                                   readable      = F,
                                   pool = pool)


  df <- as.data.frame(head(ego,5000))
  rownames(df) <- NULL
  if(ontology == "ALL" & "ONTOLOGY" %in% colnames(df)){
    df <- df %>% dplyr::rename(db = ONTOLOGY)
  }else{
    df <- df %>% dplyr::mutate(db = ontology)
  }
  df <- df %>%
    dplyr::mutate(db = paste0("GO_",db)) %>%
    dplyr::rename(go_id = ID, go_description = Description, count = Count,
                  gene_ratio = GeneRatio, background_ratio = BgRatio, gene_id = geneID) %>%
    dplyr::mutate(go_description_link = paste0('<a href=\'http://amigo.geneontology.org/amigo/term/',go_id,'\' target=\'_blank\'>',go_description,'</a>')) %>%
    tidyr::separate(gene_ratio,c('num_query_hits','num_query_all'),sep='/',remove = F, convert = T) %>%
    dplyr::mutate(qvalue = as.numeric(formatC(qvalue, format = "e", digits = 1))) %>%
    dplyr::mutate(pvalue = as.numeric(formatC(pvalue, format = "e", digits = 1))) %>%
    tidyr::separate(background_ratio,c('num_background_hits','num_background_all'),sep='/',remove = F, convert = T) %>%
    dplyr::mutate(fold_enrichment = round(as.numeric((num_query_hits/num_query_all) / (num_background_hits/num_background_all)),digits = 1)) %>%
    dplyr::select(-c(num_query_hits,num_query_all,num_background_hits,num_background_all))

  gene2id <- NULL
  if(!is.null(genedb)){
    gene2id <- as.data.frame(
      df %>%
        dplyr::select(go_id,gene_id) %>%
        tidyr::separate_rows(gene_id,sep="/") %>%
        dplyr::mutate(gene_id = as.character(gene_id)) %>%
        dplyr::left_join(dplyr::select(genedb,entrezgene,symbol), by=c("gene_id" = "entrezgene")) %>%
        dplyr::mutate(entrez_url = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",gene_id,"' target=\'blank_\'>",symbol,"</a>")) %>%
        dplyr::group_by(go_id) %>%
        dplyr::summarise(gene_symbol_link = paste(unique(entrez_url),collapse="|"), gene_symbol = paste(unique(symbol),collapse="|"))
    )
  }
  if(!is.null(gene2id)){
    df <- df %>% dplyr::left_join(dplyr::select(gene2id,go_id,gene_symbol_link,gene_symbol), by="go_id") %>%
      dplyr::rename(exact_source = go_id, description = go_description, description_link = go_description_link) %>%
      dplyr::mutate(standard_name = exact_source)
  }

  return(df)

}


get_universal_enrichment <- function(query_entrez,
                                     background_entrez = NULL,
                                     genedb = NULL,
                                     p_value_cutoff = 0.05,
                                     p_value_adjustment_method = "BH",
                                     minGSSize = 10,
                                     q_value_cutoff = 0.2,
                                     TERM2GENE = NULL,
                                     TERM2NAME = NULL,
                                     TERM2SOURCE = NULL,
                                     dbsource = ""){

  rlogging::message(paste0("Enrichment - ",dbsource,": performing gene enrichment analysis of query set"))
  rlogging::message(paste0("Enrichment - ",dbsource,": Settings: p_value_cutoff = ",p_value_cutoff,", q_value_cutoff = ",q_value_cutoff))
  rlogging::message(paste0("Enrichment - ",dbsource,": Settings: minGSSize = ",minGSSize))
  stopifnot(is.character(query_entrez))
  stopifnot(!is.null(genedb) | !is.data.frame(genedb))
  stopifnot(!is.null(TERM2SOURCE) | !is.data.frame(TERM2SOURCE))
  stopifnot("standard_name" %in% colnames(TERM2SOURCE))
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
                                   minGSSize     = minGSSize,
                                   pvalueCutoff  = p_value_cutoff,
                                   qvalueCutoff  = q_value_cutoff,
                                   TERM2GENE     = TERM2GENE,
                                   TERM2NAME     = TERM2NAME)


  df <- as.data.frame(head(enr,5000))
  rownames(df) <- NULL
  if(nrow(df) > 0){
    df <- df %>%
      dplyr::rename(standard_name = ID, description = Description, count = Count,
                    gene_ratio = GeneRatio, background_ratio = BgRatio, gene_id = geneID) %>%
      dplyr::left_join(TERM2SOURCE,by="standard_name") %>%
      dplyr::mutate(description_link = dplyr::if_else(!is.na(external_url),paste0("<a href='",external_url,"' target='_blank'>",description,"</a>"),description)) %>%
      dplyr::mutate(qvalue = as.numeric(formatC(qvalue, format = "e", digits = 1))) %>%
      dplyr::mutate(pvalue = as.numeric(formatC(pvalue, format = "e", digits = 1))) %>%
      tidyr::separate(gene_ratio,c('num_query_hits','num_query_all'),sep='/',remove = F, convert = T) %>%
      tidyr::separate(background_ratio,c('num_background_hits','num_background_all'),sep='/',remove = F, convert = T) %>%
      dplyr::mutate(fold_enrichment = round(as.numeric((num_query_hits/num_query_all) / (num_background_hits/num_background_all)),digits = 1)) %>%
      dplyr::select(-c(num_query_hits,num_query_all,num_background_hits,num_background_all)) %>%
      dplyr::mutate(db = dplyr::if_else(is.na(db) & nchar(dbsource) > 0,dbsource,as.character(db)))

    gene2id <- NULL
    if(!is.null(genedb)){
      gene2id <- as.data.frame(
        df %>%
          dplyr::select(standard_name,gene_id) %>%
          tidyr::separate_rows(gene_id,sep="/") %>%
          dplyr::mutate(gene_id = as.character(gene_id)) %>%
          dplyr::left_join(dplyr::select(genedb,entrezgene,symbol), by=c("gene_id" = "entrezgene")) %>%
          dplyr::mutate(entrez_url = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",gene_id,"' target=\'blank_\'>",symbol,"</a>")) %>%
          dplyr::group_by(standard_name) %>%
          dplyr::summarise(gene_symbol_link = paste(unique(entrez_url),collapse="|"), gene_symbol = paste(unique(symbol),collapse="|"))
      )
    }
    if(!is.null(gene2id)){
      df <- df %>% dplyr::left_join(dplyr::select(gene2id,standard_name,gene_symbol_link,gene_symbol), by="standard_name")
    }
  }else{
    df <- NULL
  }
  return(df)
}
