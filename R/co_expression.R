
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

gtex_co_expression <- function(qgenes, genedb = NULL, gids = NULL,
                               dbdir = "/Volumes/sigven/research/DB/GTex/co_expression") {
  stopifnot(!is.null(gids))
  stopifnot(is.character(qgenes))
  stopifnot(!is.null(genedb))
  stopifnot("symbol" %in% colnames(genedb) & "entrezgene" %in% colnames(genedb))
  if (!dir.exists(dbdir)) {
    rlogging::warning(paste0("Database directory with GTex co-expression matrices (", dbdir,
                             ") does not exist - skipping"))
  }
  co_exp_res <- list()
  co_exp_res[["df"]] <- data.frame()
  co_exp_res[["plots"]] <- list()
  df <- data.frame("entrezgene" = qgenes, stringsAsFactors = F) %>%
    dplyr::left_join(dplyr::select(genedb, entrezgene, symbol), by = "entrezgene") %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::select(-entrezgene) %>%
    dplyr::distinct()
  for(gid in gids) {
    fname_rds <- paste0(dbdir, "/", gid, ".co_exp.rds")
    if(!file.exists(fname_rds)) {
      rlogging::warning(paste0("File with co-expression values for ",gid, " does not exist - continuing"))
    }else{
      co_exp_res[["plots"]][[gid]] <- NULL
      co_exp_data <- readRDS(fname_rds)

      part1 <- co_exp_data %>% dplyr::inner_join(df, by = c("symbol_A" = "symbol")) %>%
        dplyr::inner_join(df, by=c("symbol_B" = "symbol"))
      part2 <- co_exp_data %>% dplyr::inner_join(df, by = c("symbol_B" = "symbol")) %>%
        dplyr::inner_join(df, by=c("symbol_A" = "symbol"))

      coexp_res <- dplyr::bind_rows(part1, part2) %>%
        dplyr::distinct() %>%
        dplyr::left_join(oncoEnrichR::gtex_rnaseq_db$metadata$tissues, by = c("gid" = "AtlasAssayGroup"))

      rlogging::message(paste0("GTex: retrieved co-expressed genes from ", unique(coexp_res$tissue_name),
                               ": n = ", unique(coexp_res$n_samples)))

      all_dups <- data.frame()
      for (n in 1:nrow(coexp_res)) {
        entry <- coexp_res[n,]
        tissue_name <- coexp_res[n,]$tissue_name
        n_samples <- coexp_res[n,]$n_samples
        n_dup <- data.frame(symbol_A = entry$symbol_B, symbol_B = entry$symbol_A,
                            r = entry$r, p_value = entry$p_value, gid = entry$gid,
                            tissue_name = tissue_name, n_samples = n_samples, stringsAsFactors = F)
        n_dup2 <- data.frame(symbol_A = entry$symbol_B, symbol_B = entry$symbol_B, r = 1,
                             p_value = 0, gid = entry$gid,  tissue_name = tissue_name,
                             n_samples = n_samples, stringsAsFactors = F)
        n_dup3 <- data.frame(symbol_A = entry$symbol_A, symbol_B = entry$symbol_A, r = 1,
                             p_value = 0, gid = entry$gid,  tissue_name = tissue_name,
                             n_samples = n_samples, stringsAsFactors = F)
        all_dups <- dplyr::bind_rows(all_dups, n_dup, n_dup2, n_dup3)
      }

      coexp_res <- dplyr::bind_rows(coexp_res, all_dups) %>% dplyr::distinct()

      m <- as.matrix(reshape2::acast(coexp_res, symbol_A~symbol_B, value.var = "r"))
      melted_cormat <- reshape2::melt(get_upper_tri(m))

      title <- paste0(unique(coexp_res$tissue_name), ": n = ", unique(coexp_res$n_samples))
      co_exp_res[["plots"]][[gid]] <-
        ggplot2::ggplot(data = melted_cormat, ggplot2::aes(Var2, Var1, fill = value)) +
             ggplot2::geom_tile(color = "grey")+
             ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "grey",
                                          midpoint = 0, limit = c(-1, 1),
                                          name = "Spearman\nCorrelation") +
            ggplot2::theme_classic()+
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 9)) +
            ggplot2::xlab("") +
            ggplot2::ylab("") +
            ggplot2::ggtitle(title)


      co_exp_res[["df"]] <- dplyr::bind_rows(co_exp_res[["df"]], coexp_res)
      rm(co_exp_data)

    }
  }
  return(co_exp_res)
}

gtex_rnaseq_db <- NULL
# gtex_rnaseq_db <- list()
# gtex_rnaseq_db[['metadata']] <- list()
# gtex_rnaseq_db[['metadata']][['analysis']] <- metadata(sumexp)
# gtex_rnaseq_db[['metadata']][['samples']] <- as.data.frame(SummarizedExperiment::colData(sumexp))
# gtex_rnaseq_db[['metadata']][['samples']]$sample_id <- rownames(gtex_rnaseq_db[['metadata']][['samples']])
# rownames(gtex_rnaseq_db[['metadata']][['samples']]) <- NULL
#
# gtex_rnaseq_db[['metadata']][['tissues']] <- as.data.frame(
#   gtex_rnaseq_db$metadata$samples %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g9","Whole blood",as.character(NA))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g28","Skin (sun exposed)",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g44","Skin (non sun exposed)",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g10","Breast",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g45","Testis",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g33","Pancreas",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g29","Lung",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g50","Colon - transverse",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g37","Colon - sigmoid",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g41","Stomach",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g38","Muscle - skeletal",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g42","Adipose tissue",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g46","Thyroid",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g5","Adrenal gland",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g27","Liver",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g35","Prostate",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g32","Ovary",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g20","Espophagus - mucosa",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g19","Espophagus - gastric junction",as.character(tissue_name))) %>%
#     dplyr::mutate(tissue_name = dplyr::if_else(AtlasAssayGroup == "g42","Adipose tissue",as.character(tissue_name))) %>%
#     dplyr::filter(!is.na(tissue_name)) %>%
#     dplyr::group_by(tissue_name, AtlasAssayGroup) %>%
#     dplyr::summarise(n_samples = dplyr::n()) %>%
#     dplyr::arrange(desc(n_samples))
# )
#
#
# gene_info <- readRDS(file="data-raw/gene_info.rds")
#
# ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")
# ensembl_genes <- biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
# queryAttributes <- c('ensembl_gene_id','hgnc_id')
# ensembl_genes_xref <- as.data.frame(biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes) %>%
#                                       dplyr::mutate(hgnc_id = stringr::str_replace(hgnc_id,"HGNC:","")) %>%
#                                       dplyr::left_join(dplyr::select(gene_info, hgnc_id, symbol_entrez, gene_biotype)) %>%
#                                       dplyr::rename(symbol = symbol_entrez) %>%
#                                       dplyr::filter(gene_biotype == 'protein_coding') %>%
#                                       dplyr::distinct() %>%
#                                       dplyr::group_by(ensembl_gene_id) %>%
#                                       dplyr::summarise(symbol = paste(unique(symbol), collapse="&")))
#
# gtex_rnaseq_db[['counts']] <- list()
# for(gid in unique(gtex_rnaseq_db[['metadata']][['tissues']]$AtlasAssayGroup)) {
#   samples <- dplyr::filter(gtex_rnaseq_db[['metadata']][['samples']], AtlasAssayGroup == gid) %>%
#     dplyr::select(sample_id, AtlasAssayGroup)
#
#   gtex_rnaseq_db[['counts']][[gid]] <- list()
#   gtex_rnaseq_db[['counts']][[gid]][['counts']] <-
#     as.matrix(SummarizedExperiment::assays( sumexp )$counts)[,samples$sample_id]
#   cat(gid,'\n')
#
#   df <- data.frame('ensembl_gene_id' = rownames(gtex_rnaseq_db[['counts']][[gid]][['counts']]), stringsAsFactors = F)
#   df <- df %>% dplyr::left_join(ensembl_genes_xref,by="ensembl_gene_id")
#   rownames(gtex_rnaseq_db[['counts']][[gid]][['counts']]) <- df$symbol
#   gtex_rnaseq_db[['counts']][[gid]][['counts']] <-
#     as.matrix(gtex_rnaseq_db[['counts']][[gid]][['counts']][!is.na(rownames(gtex_rnaseq_db[['counts']][[gid]][['counts']])), ])
#
# }
#
# saveRDS(gtex_rnaseq_db,file="data-raw/gtex_rnaseq_db.rds")
# gtex_rnaseq_db$counts <- NULL
# if(!is.null(gtex_rnaseq_db)) {
#   usethis::use_data(gtex_rnaseq_db, overwrite = T)
#   rm(gtex_rnaseq_db)
# }

