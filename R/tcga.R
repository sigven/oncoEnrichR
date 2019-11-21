
tcga_aberration_plot <- function(qgenes, qsource = "symbol", cstrata = "site", vtype = "snv_indel", genedb = NULL, percentile = FALSE){

  rlogging::message(paste0("TCGA: generating gene aberration plot, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == 'entrezgene'){
    query_genes_df <- data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "entrezgene") %>% dplyr::distinct()
  }else{
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "symbol") %>% dplyr::distinct()
  }

  title <- 'SNVs/InDels - TCGA'
  color <- 'steelblue'
  if(vtype == 'cna_ampl'){
    title <- 'Copy number amplifications (sCNA) - TCGA'
    color <- 'darkgreen'
  }
  if(vtype == 'cna_homdel'){
    title <- 'Homozygous deletions - (sCNA) - TCGA'
    color <- 'firebrick'
  }

  tcga_gene_stats <- oncoEnrichR::tcga_aberration_stats %>%
    dplyr::filter(clinical_strata == cstrata) %>%
    dplyr::filter(primary_site != "Other/Unknown") %>%
    dplyr::inner_join(dplyr::select(query_genes_df, symbol),by=c("symbol"))


  gene_candidates_init <- data.frame()
  tcga_ttypes <- sort(unique(tcga_gene_stats$primary_site))
  for(i in 1:length(sort(unique(tcga_gene_stats$symbol)))){
    init <- data.frame('primary_site' <- tcga_ttypes, 'primary_diagnosis_very_simplified' = NA,
                       'symbol' = sort(unique(tcga_gene_stats$symbol))[i], 'genomic_strata' = 'gene',
                       'clinical_strata' = cstrata, 'percent_mutated' = 0, 'percentile' = 0,
                       'variant_type' = vtype, 'consensus_calls' = F, 'fp_driver_gene' = as.logical(NA),
                       decile = 0, stringsAsFactors = F)
    colnames(init) <- c('primary_site','primary_diagnosis_very_simplified','symbol','genomic_strata',
                        'clinical_strata','percent_mutated','percentile','variant_type', 'consensus_calls',
                        'fp_driver_gene','decile')
    gene_candidates_init <- rbind(gene_candidates_init, init)
    i <- i + 1
  }
  gene_candidates_init <- gene_candidates_init %>%
    dplyr::filter(primary_site != 'Other/Unknown' & primary_site != "Pancancer")


  gene_aberrations <- tcga_gene_stats %>%
    dplyr::filter(variant_type == vtype) %>%
    dplyr::filter(primary_site != "Pancancer" & primary_site != "Other/Unknown")


  site_stats_zero <- tcga_gene_stats %>%
    dplyr::select(primary_site,tot_samples) %>%
    dplyr::distinct() %>%
    dplyr::mutate(samples_mutated = 0)

  pancan_order <- tcga_gene_stats %>%
    dplyr::filter(variant_type == vtype & primary_site == "Pancancer") %>%
    dplyr::mutate(pancancer_percent_mutated = percent_mutated) %>%
    dplyr::mutate(pancancer_percentile = percentile) %>%
    dplyr::select(symbol, pancancer_percent_mutated, pancancer_percentile)

  zero_frequency_genes <- dplyr::anti_join(gene_candidates_init, gene_aberrations, by=c("symbol","primary_site","variant_type")) %>%
    #dplyr::select(-c(tot_samples,samples_mutated)) %>%
    dplyr::left_join(site_stats_zero,by=c("primary_site"))
  gene_aberrations <- dplyr::left_join(dplyr::bind_rows(gene_aberrations, zero_frequency_genes), pancan_order,by=c("symbol"))


  gene_aberrations <- gene_aberrations %>%
    dplyr::arrange(pancancer_percent_mutated) %>%   # rearrange the df in the order we want
    dplyr::mutate(symbol = factor(symbol, unique(symbol)))

  p <- ggplot2::ggplot(gene_aberrations,ggplot2::aes(x=primary_site,y=symbol)) +
    ggplot2::geom_text(ggplot2::aes(label = round(percent_mutated,1)), color = "#E69F00") +
    ggplot2::geom_tile(ggplot2::aes(fill = percent_mutated), colour="black",size=0.40)

  if(percentile == T){
    p <- ggplot2::ggplot(gene_aberrations,ggplot2::aes(x=primary_site,y=symbol)) +
      ggplot2::geom_text(ggplot2::aes(label = round(percentile,1))) +
      ggplot2::geom_tile(ggplot2::aes(fill = percentile), colour="black",size=0.40)
  }
  p <- p +
    #add border white colour of line thickness 0.25
    #ggplot2::geom_tile(aes(fill = pct_mutated), colour="black",size=0.40)+
    ggplot2::scale_fill_gradient(low = "white", high = color) +
    #remove x and y axis labels
    ggplot2::labs(x="",y="")+
    ggplot2::ggtitle(title) +
    #remove extra space
    #scale_y_discrete(expand=c(0,0))+
    #ggplot2::coord_fixed(ratio = 1.3)+
    #set a base size for all fonts
    ggplot2::theme_grey(base_size=14)+
    #theme options
    ggplot2::theme(
      #bold font for both axis text
      #legend.title = ggplot2::element_text()
      axis.text = ggplot2::element_text(face="bold",family = "Helvetica", size = 14),
      #set thickness of axis ticks
      axis.ticks = ggplot2::element_line(size=0.2),
      #remove plot background
      plot.background = ggplot2::element_blank(),
      #remove plot border
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}

tcga_aberration_table <- function(qgenes, qsource = "entrezgene", genedb = NULL, vtype = "snv_indel"){

  rlogging::message(paste0("TCGA: collecting gene aberration data table, variant type =  ",vtype))
  stopifnot(!is.null(genedb))
  oncoEnrichR::validate_db_df(genedb, dbtype = "genedb")
  stopifnot(qsource == "symbol" | qsource == "entrezgene")
  stopifnot(is.character(qgenes))
  query_genes_df <- data.frame('symbol' = qgenes, stringsAsFactors = F)
  if(qsource == 'entrezgene'){
    query_genes_df <- data.frame('entrezgene' = qgenes, stringsAsFactors = F)
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "entrezgene") %>% dplyr::distinct()
  }else{
    query_genes_df <- dplyr::inner_join(genedb, query_genes_df, by = "symbol") %>% dplyr::distinct()
  }

  aberration_data <- oncoEnrichR::tcga_aberration_stats %>%
    dplyr::filter(clinical_strata == "site_diagnosis" & variant_type == vtype) %>%
    dplyr::filter(primary_site != "Pancancer") %>%
    dplyr::select(symbol, primary_site, primary_diagnosis_very_simplified, variant_type,
                  samples_mutated, tot_samples, percent_mutated) %>%
    dplyr::rename(primary_diagnosis = primary_diagnosis_very_simplified, cohort_size = tot_samples) %>%
    dplyr::left_join(dplyr::select(genedb, symbol, entrezgene),by=c("symbol")) %>%
    dplyr::inner_join(dplyr::select(query_genes_df, entrezgene),by=c("entrezgene")) %>%
    dplyr::distinct() %>%
    dplyr::mutate(gene = paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",entrezgene,"' target='_blank'>",symbol,"</a>")) %>%
    dplyr::select(-c(entrezgene,symbol)) %>%
    dplyr::select(gene,primary_site,primary_diagnosis,variant_type,dplyr::everything())

  return(aberration_data)
}
