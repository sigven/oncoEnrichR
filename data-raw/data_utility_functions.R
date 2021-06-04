
#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed <- function(pmid){

  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmid,ceiling(length(pmid)/400))
  j <- 0
  all_citations <- data.frame()
  cat('Retrieving PubMed citations for PMID list, total length', length(pmid))
  cat('\n')
  while(j < length(pmid_chunks)){
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    cat('Processing chunk ',j,' with ',length(pmid_chunk),'PMIDS')
    cat('\n')
    pmid_string <- paste(pmid_chunk,collapse = " ")
    res <- RISmed::EUtilsSummary(pmid_string, type="esearch", db="pubmed", retmax = 5000)
    year <- RISmed::YearPubmed(RISmed::EUtilsGet(res))
    authorlist <- RISmed::Author(RISmed::EUtilsGet(res))
    pmid_list <- RISmed::PMID(RISmed::EUtilsGet(res))
    i <- 1
    first_author <- c()
    while(i <= length(authorlist)){
      first_author <- c(first_author,paste(authorlist[[i]][1,]$LastName," et al.",sep=""))
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(RISmed::EUtilsGet(res))
    citations <- data.frame('pmid' = as.integer(pmid_list), 'citation' = paste(first_author,year,journal,sep=", "), stringsAsFactors = F)
    citations$link <- paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',citations$pmid,'\' target=\'_blank\'>',citations$citation,'</a>')
    all_citations <- dplyr::bind_rows(all_citations, citations)
    j <- j + 1
  }

  return(all_citations)

}

get_gene_info_ncbi <- function(basedir = NULL,
                               update = T){

  invisible(assertthat::assert_that(
    update == T | update == F,
    msg = "'update' is not of type logical (TRUE/FALSE)"))
  invisible(assertthat::assert_that(
    dir.exists(basedir),
    msg = paste0("Directory '",
                 basedir,"' does not exist")))

  rlogging::message("Retrieving gene_info from NCBI/Entrez")
  gene_info_fname <- paste0(basedir,'/data-raw/gene_info/Homo_sapiens.gene_info.gz')
  if(update == F){
    invisible(assertthat::assert_that(
      file.exists(gene_info_fname),
      msg = paste0("File ",gene_info_fname,
                   " does not exist")))
  }

  if(!file.exists(gene_info_fname) | update == T){
    remote_url <- "ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
    if(RCurl::url.exists(remote_url)){
      download.file(
        remote_url, destfile = gene_info_fname, quiet=T)
    }else{
      rlogging::message(
        paste0("Cannot update gene_info - download not available: ",
               remote_url))
    }
  }
  gene_info <- read.table(
    gzfile(gene_info_fname),
    sep = "\t", stringsAsFactors = F,na.strings = "-",
    skip = 1, comment.char = "#", quote = "", fill = T) %>%
    dplyr::filter(V1 == 9606) %>%
    dplyr::select(c(V2, V3, V5, V9, V6, V10, V11)) %>%
    dplyr::rename(
      entrezgene = V2, synonyms = V5, symbol = V3, name = V9,
      gene_biotype = V10, symbol_entrez = V11) %>%
    dplyr::mutate(ensembl_gene_id = stringr::str_replace(
      stringr::str_match(V6,"Ensembl:ENSG[0-9]{1,}"), "Ensembl:", "")) %>%
    dplyr::mutate(hgnc_id = stringr::str_replace(
      stringr::str_match(V6,"HGNC:HGNC:[0-9]{1,}"), "HGNC:HGNC:", "")) %>%
    dplyr::select(-V6) %>%
    dplyr::mutate(symbol_entrez = dplyr::if_else(
      is.na(symbol_entrez) | symbol_entrez == "-",
      symbol, as.character(symbol_entrez))) %>%
    dplyr::mutate(gene_biotype = dplyr::if_else(
      gene_biotype ==  "protin-coding",
      "protein_coding", as.character(gene_biotype))) %>%
    dplyr::filter(nchar(symbol_entrez) > 0)

  ### for genes annotated with the same ensembl gene ids, ignore this annotation (set to NA)
  ensgene_id_count <- as.data.frame(
    dplyr::filter(gene_info, !is.na(ensembl_gene_id)) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
  )

  gene_info <- gene_info %>%
    dplyr::left_join(ensgene_id_count,by=c("ensembl_gene_id")) %>%
    dplyr::mutate(ensembl_gene_id = dplyr::if_else(
      !is.na(n) & n > 1,
      as.character(NA),
      as.character(ensembl_gene_id))) %>%
    dplyr::select(-n)

  # if (update == T) {
  #   saveRDS(gene_info, file = paste0(basedir, "/data-raw/ncbi_gene/gene_info.rds"))
  # }

  return(gene_info)

}

get_gene_aliases <- function(gene_info = NULL){

  invisible(assertthat::assert_that(!is.null(gene_info)))
  assertable::assert_colnames(gene_info, colnames = c("synonyms","symbol"),
                              only_colnames = F, quiet = T)

  unique_synonyms <- as.data.frame(
    gene_info %>%
      dplyr::select(synonyms, symbol) %>%
      tidyr::separate_rows(synonyms,sep="\\|") %>%
      dplyr::group_by(synonyms) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n == 1)
  ) %>%
    dplyr::select(-n)

  alias2primary <- as.data.frame(
    gene_info %>%
      dplyr::select(synonyms, symbol) %>%
      tidyr::separate_rows(synonyms,sep="\\|") %>%
      dplyr::inner_join(unique_synonyms, by = "synonyms") %>%
      dplyr::rename(alias = synonyms)
  )

  return(alias2primary)
}

get_msigdb_signatures <- function(basedir = NULL,
                                  db_version = 'v7.4 (April 2021)'){

  ## get full dataset: Broad Institute's Molecular Signatures Database (v7.1)
  invisible(assertthat::assert_that(
    dir.exists(basedir),
    msg = paste0("Directory '",
                 basedir,"' does not exist")))
  msigdb_xml_fname <- file.path(
    basedir, "data-raw","msigdb","msigdb.xml")
  invisible(assertthat::assert_that(
    file.exists(msigdb_xml_fname),
    msg = paste0("File '",
                 msigdb_xml_fname,
                 "' does not exist")))

  msig_data_xml <- xml2::read_xml(msigdb_xml_fname)

  ## make data frame with signatures, one record pr. gene-signature association
  all_genesets <- msig_data_xml %>% xml2::xml_find_all("//GENESET")
  category_code <- all_genesets %>% xml2::xml_attr("CATEGORY_CODE")
  all_msigdb <- data.frame('category_code' = category_code, stringsAsFactors = F)
  all_msigdb$description <- all_genesets %>% xml2::xml_attr("DESCRIPTION_BRIEF")
  all_msigdb$standard_name <- all_genesets %>% xml2::xml_attr("STANDARD_NAME")
  all_msigdb$organism <- all_genesets %>% xml2::xml_attr("ORGANISM")
  all_msigdb$pmid <- all_genesets %>% xml2::xml_attr("PMID")
  all_msigdb$systematic_name <- all_genesets %>% xml2::xml_attr("SYSTEMATIC_NAME")
  all_msigdb$subcategory_code <- all_genesets %>% xml2::xml_attr("SUB_CATEGORY_CODE")
  all_msigdb$entrezgene <- all_genesets %>% xml2::xml_attr("MEMBERS_EZID")
  all_msigdb$contributor <- all_genesets %>% xml2::xml_attr("CONTRIBUTOR")
  all_msigdb$exact_source <- all_genesets %>% xml2::xml_attr("EXACT_SOURCE")
  all_msigdb$external_url <- all_genesets %>% xml2::xml_attr("EXTERNAL_DETAILS_URL")
  all_msigdb <- all_msigdb %>%
    tidyr::separate_rows(entrezgene,sep=",") %>%
    dplyr::filter(organism == "Homo sapiens") %>%
    dplyr::filter(subcategory_code != 'CP:KEGG') %>%
    dplyr::arrange(category_code) %>%
    dplyr::mutate(pmid = dplyr::if_else(nchar(pmid) == 0,as.character(NA),as.character(pmid))) %>%
    dplyr::mutate(subcategory_code = dplyr::if_else(nchar(subcategory_code) == 0,as.character("ALL"),as.character(subcategory_code))) %>%
    dplyr::filter(category_code != "ARCHIVED" & category_code != "C1") %>%
    dplyr::mutate(external_url = dplyr::if_else(nchar(external_url) == 0,paste0("http://software.broadinstitute.org/gsea/msigdb/cards/",standard_name),as.character(external_url))) %>%
    dplyr::mutate(description = stringr::str_replace_all(description,
                                                         "( \\[ICI [0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[GeneID=[0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[PubChem=[0-9]{1,}(;[0-9]{1,})*\\]( )?)",""))

  msigdb_category_description <- read.table(file="data-raw/msigdb/msigdb_collection_description.tsv",
                                   sep = "\t", header = T, stringsAsFactors = F)

  msigdb_complete <- as.data.frame(all_msigdb %>%
    dplyr::left_join(msigdb_category_description,
                     by = c("category_code", "subcategory_code")) %>%
    dplyr::mutate(db = "MSigDB", db_version = db_version) %>%
    dplyr::select(db, db_version, category_code, category_description,
                  subcategory_code, subcategory_description,
                  standard_name, description, organism,
                  entrezgene) %>%
    dplyr::rename(signature_description = description,
                  signature_name = standard_name)
  )


  msigdb <- list()
  msigdb[['VERSION']] <- version
  msigdb[['TERM2SOURCE']] <- all_msigdb %>%
    dplyr::filter(subcategory_code != "MIR:MIR_Legacy") %>%
    dplyr::filter(subcategory_code != "TFT:TFT_Legacy") %>%
    dplyr::filter(subcategory_code != "CP:WIKIPATHWAYS") %>%
    dplyr::filter(subcategory_code != "VAX") %>%
    dplyr::select(standard_name, exact_source, external_url, category_code, subcategory_code) %>%
    dplyr::distinct() %>%
    dplyr::mutate(external_url = stringr::str_replace(external_url,"\\|https://reactome.org/PathwayBrowser/","")) %>%
    dplyr::mutate(external_url = dplyr::if_else(nchar(external_url) == 0,as.character(NA),as.character(external_url))) %>%
    dplyr::mutate(db = dplyr::if_else(stringr::str_detect(standard_name,"^GO_"),subcategory_code,as.character(NA))) %>%
    dplyr::mutate(db = dplyr::if_else(stringr::str_detect(standard_name,"^HP_"),subcategory_code,as.character(NA))) %>%
    dplyr::mutate(db = dplyr::if_else(stringr::str_detect(standard_name,"^REACTOME_"),"REACTOME",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(stringr::str_detect(standard_name,"^BIOCARTA_"),"BIOCARTA",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(stringr::str_detect(standard_name,"^PID_"),"PATHWAY_INTERACTION_DB",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "H","HALLMARK",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CGP","CHEM_GEN_PERTURB",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CP","CANONICAL_PATHWAY",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "MIR:MIRDB","MICRORNA_TARGET_MIRDB",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "TFT:GTRD","TF_TARGET_GTRD",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CGN","CANCER_NEIGHBOURHOOD",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CM","CANCER_MODULE",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C6","ONCOGENIC",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C7","IMMUNESIGDB",as.character(db))) %>%
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C8","CELLTYPE_SIGNATURES",as.character(db))) %>%
    dplyr::select(-c(category_code,subcategory_code)) %>%
    dplyr::filter(!is.na(db))



  msigdb[['COLLECTION']] <- list()

  for(c in c('H','C2','C3','C4','C5','C6','C7','C8')){
    msigdb[['COLLECTION']][[c]] <- list()
    subcategories <- unique(all_msigdb[all_msigdb$category_code == c,]$subcategory_code)
    for(scat in subcategories){
      subcat <- stringr::str_replace(scat,"GO:","")
      if(subcat == "CP:WIKIPATHWAYS" |
         subcat == "VAX" |
         subcat == "MIR:MIR_Legacy" |
         subcat == "TFT:TFT_Legacy"){
        next
      }
      #cat(c,subcat,sep=" - ")
      #cat('\n')
      msigdb[['COLLECTION']][[c]][[subcat]] <- list()
      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2GENE']] <- data.frame()
      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2GENE']] <-
        dplyr::filter(all_msigdb, category_code == c & subcategory_code == scat) %>%
        dplyr::select(standard_name, entrezgene) %>%
        dplyr::rename(entrez_gene = entrezgene)
      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2NAME']] <-
        dplyr::filter(all_msigdb, category_code == c & subcategory_code == scat) %>%
        dplyr::select(standard_name, description) %>%
        dplyr::rename(name = description) %>%
        dplyr::distinct()
    }
  }
  rm(all_msigdb)
  return(list('db' = msigdb, 'df' = msigdb_complete))
}


get_ts_oncogene_annotations <- function(basedir = NULL,
                                        version = "34",
                                        gene_info = NULL){
  rlogging::message(
    "Retrieving proto-oncogenes/tumor suppressor status ",
    "according to classifications form CancerMine & Network of Cancer Genes (NCG)")

  fp_cancer_drivers <- get_curated_fp_cancer_genes(basedir = basedir, gene_info = gene_info)

  ## Tumor suppressor annotations or oncogene annotations from Network of cancer genes (NCG)
  ncg <- read.table(file = "data-raw/ncg/ncg_tsgene_oncogene.tsv",header = T,
                    stringsAsFactors = F, sep = "\t", quote = "", comment.char = "") %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      ncg_tumor_suppressor = dplyr::if_else(stringr::str_detect(cgc_annotation,"TSG") | ncg6_tsg == 1,
                                            as.logical(TRUE),as.logical(FALSE))) %>%
    dplyr::mutate(
      ncg_oncogene = dplyr::if_else(stringr::str_detect(cgc_annotation,"oncogene") | ncg6_oncogene == 1,
                                    as.logical(TRUE),as.logical(FALSE))) %>%
    dplyr::select(symbol, entrez, ncg_tumor_suppressor, ncg_oncogene) %>%
    dplyr::anti_join(dplyr::select(fp_cancer_drivers, symbol), by = "symbol") %>%
    dplyr::select(-entrez) %>%
    dplyr::distinct()

  pmids <- as.data.frame(
    read.table(gzfile(paste0(basedir,'/data-raw/cancermine/cancermine_sentences.tsv.gz')),
               header = T, comment.char = "", quote = "", sep="\t", stringsAsFactors = F) %>%
      dplyr::filter(predictprob >= 0.9) %>%
      dplyr::group_by(role, gene_normalized, pmid) %>%
      dplyr::summarise(doid = paste(unique(cancer_id), collapse=","), .groups = "drop") %>%
      dplyr::rename(symbol = gene_normalized) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, symbol), by = "symbol") %>%
      dplyr::distinct()
  )

  all_citations <- data.frame()
  if(file.exists(paste0(basedir,"/data-raw/cancermine/citations/cancermine_citations.v",version,".tsv"))){
    all_citations <- read.table(file=paste0(basedir,"/data-raw/cancermine/citations/cancermine_citations.v",version,".tsv"),
                                sep = "\t", header = F, quote = "", comment.char = "", stringsAsFactors = F) %>%
      magrittr::set_colnames(c('pmid','citation','citation_link')) %>%
      dplyr::distinct() %>%
      dplyr::arrange(desc(pmid))

  }

  pmids <- pmids %>%
    dplyr::inner_join(all_citations, by=c("pmid"))

  pmids_oncogene <- as.data.frame(
    pmids %>%
      dplyr::filter(role == "Oncogene") %>%
      dplyr::arrange(symbol,desc(pmid)) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(pmids_oncogene = paste(pmid,collapse=", "),
                       citations_oncogene = paste(citation,collapse="; "),
                       citation_links_oncogene = paste(head(citation_link,50),collapse=", "),
                       .groups = "drop") %>%
      dplyr::mutate(n_citations_oncogene = as.integer(stringr::str_count(pmids_oncogene,", ")) + 1) %>%
      dplyr::filter(n_citations_oncogene > 1) %>%
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(n_citations_oncogene > 4,
                                                         "MC",
                                                         as.character("LC"))) %>%
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(n_citations_oncogene >= 10,
                                                         "HC",
                                                         as.character(oncogene_cancermine)))
  )

  pmids_tsgene <- as.data.frame(
    pmids %>%
      dplyr::filter(role == "Tumor_Suppressor") %>%
      dplyr::arrange(symbol,desc(pmid)) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(pmids_tsgene = paste(unique(pmid),collapse=", "),
                       citations_tsgene = paste(citation,collapse="; "),
                       citation_links_tsgene = paste(head(citation_link,50),collapse=", "),
                       .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_citations_tsgene = as.integer(stringr::str_count(pmids_tsgene,", ")) + 1) %>%
      dplyr::filter(n_citations_tsgene > 1) %>%
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(n_citations_tsgene > 4,
                                                                 "MC",
                                                                 as.character("LC"))) %>%
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(n_citations_tsgene >= 10,
                                                                 "HC",
                                                                 as.character(tumor_suppressor_cancermine)))
  )

  pmids_cdriver <- as.data.frame(
    pmids %>%
      dplyr::filter(role == "Driver") %>%
      dplyr::arrange(symbol,desc(pmid)) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(pmids_cdriver = paste(pmid,collapse=", "),
                       citations_cdriver = paste(citation,collapse="; "),
                       citation_links_cdriver = paste(head(citation_link,50),collapse=", "),
                       .groups = "drop") %>%
      dplyr::mutate(n_citations_cdriver = as.integer(stringr::str_count(pmids_cdriver,", ")) + 1) %>%
      dplyr::filter(n_citations_cdriver > 1) %>%
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(n_citations_cdriver > 4,
                                                              "MC",
                                                              as.character("LC"))) %>%
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(n_citations_cdriver >= 10,
                                                              "HC",
                                                              as.character(cancer_driver_cancermine)))
  )


  rlogging::message("Retrieving known proto-oncogenes/tumor suppressor genes from CancerMine")
  oncogene <- as.data.frame(readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                                            col_names = T,na ="-", comment = "#", quote = "") %>%
                              dplyr::filter(role == "Oncogene") %>%
                              dplyr::rename(symbol = gene_normalized) %>%
                              dplyr::anti_join(dplyr::select(fp_cancer_drivers, symbol), by = "symbol") %>%
                              dplyr::group_by(symbol) %>%
                              dplyr::summarise(doid_oncogene = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
                              dplyr::inner_join(pmids_oncogene, by = "symbol") %>%
                              dplyr::distinct())

  n_hc_oncogene <- oncogene %>%
    dplyr::filter(oncogene_cancermine == "HC") %>%
    nrow()


  tsgene <- as.data.frame(
    readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                    col_names = T,na ="-", comment = "#", quote = "") %>%
      dplyr::filter(role == "Tumor_Suppressor") %>%
      dplyr::rename(symbol = gene_normalized) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, symbol), by = "symbol") %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(doid_tsgene = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
      dplyr::inner_join(pmids_tsgene, by = "symbol") %>%
      dplyr::distinct()
  )
  n_hc_tsgene <- tsgene %>%
    dplyr::filter(tumor_suppressor_cancermine == "HC") %>%
    nrow()


  cdriver <- as.data.frame(
    readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                    col_names = T,na ="-", comment = "#", quote = "") %>%
      dplyr::filter(role == "Driver") %>%
      dplyr::rename(symbol = gene_normalized) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, symbol), by = "symbol") %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(doid_cdriver = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
      dplyr::inner_join(pmids_cdriver, by = "symbol") %>%
      dplyr::distinct()
  )
  n_hc_cdriver <- cdriver %>%
    dplyr::filter(cancer_driver_cancermine == "HC") %>%
    nrow()

  tsgene_full <- dplyr::full_join(tsgene, oncogene, by = c("symbol")) %>%
    dplyr::full_join(cdriver, by="symbol") %>%
    dplyr::full_join(ncg, by = "symbol") %>%
    dplyr::mutate(ncg_oncogene = dplyr::if_else(is.na(ncg_oncogene) | ncg_oncogene == F,FALSE,TRUE)) %>%
    dplyr::mutate(ncg_tumor_suppressor = dplyr::if_else(is.na(ncg_tumor_suppressor) | ncg_tumor_suppressor == F,FALSE,TRUE)) %>%
    dplyr::mutate(n_citations_oncogene = dplyr::if_else(is.na(n_citations_oncogene),as.integer(0),as.integer(n_citations_oncogene))) %>%
    dplyr::mutate(n_citations_tsgene = dplyr::if_else(is.na(n_citations_tsgene),as.integer(0),as.integer(n_citations_tsgene))) %>%
    dplyr::mutate(n_citations_cdriver = dplyr::if_else(is.na(n_citations_cdriver),as.integer(0),as.integer(n_citations_cdriver))) %>%

    dplyr::mutate(oncogene = dplyr::if_else((ncg_oncogene == T & !is.na(oncogene_cancermine)) | oncogene_cancermine == "HC" | oncogene_cancermine == "MC",TRUE,FALSE)) %>%
    dplyr::mutate(oncogene = dplyr::if_else(is.na(oncogene),FALSE,as.logical(oncogene))) %>%
    dplyr::mutate(tumor_suppressor = dplyr::if_else((ncg_tumor_suppressor == T & !is.na(tumor_suppressor_cancermine)) | tumor_suppressor_cancermine == "HC" | tumor_suppressor_cancermine == "MC",TRUE,FALSE)) %>%
    dplyr::mutate(tumor_suppressor = dplyr::if_else(is.na(tumor_suppressor),FALSE,as.logical(tumor_suppressor))) %>%
    dplyr::mutate(cancer_driver = dplyr::if_else(cancer_driver_cancermine == "HC" |
                                                   cancer_driver_cancermine == "MC" |
                                                   cancer_driver_cancermine == "LC",TRUE,FALSE)) %>%
    dplyr::mutate(cancer_driver = dplyr::if_else(is.na(cancer_driver),FALSE,as.logical(cancer_driver))) %>%
    dplyr::filter(tumor_suppressor == T | oncogene == T | cancer_driver == T) %>%
    dplyr::mutate(ncg_links_tsgene = dplyr::if_else(ncg_tumor_suppressor == T,
                                                    paste0("NCG: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",symbol,"\" target=\"_blank\">Tumor Suppressor</a>"),
                                                    as.character("NCG: Not defined"))) %>%
    dplyr::mutate(ncg_links_tsgene = dplyr::if_else(is.na(ncg_links_tsgene),
                                                    as.character("NCG: Not defined"),
                                                    as.character(ncg_links_tsgene))) %>%
    dplyr::mutate(ncg_links_oncogene = dplyr::if_else(ncg_oncogene == T,
                                                      paste0("NCG: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",symbol,"\" target=\"_blank\">Oncogene</a>"),
                                                      as.character("NCG: Not defined"))) %>%
    dplyr::mutate(ncg_links_oncogene = dplyr::if_else(is.na(ncg_links_oncogene),
                                                      as.character("NCG: Not defined"),
                                                      as.character(ncg_links_oncogene))) %>%
    dplyr::mutate(ncg_link = paste(ncg_links_tsgene, ncg_links_oncogene, sep=", ")) %>%
    dplyr::mutate(citation_links_cdriver = dplyr::if_else(is.na(cancer_driver),as.character(NA),as.character(citation_links_cdriver)),
                  citations_cdriver = dplyr::if_else(is.na(cancer_driver),as.character(NA),as.character(citations_cdriver)),
                  citation_links_oncogene = dplyr::if_else(is.na(oncogene_cancermine),
                                                           "Oncogenic role (CancerMine): No records",
                                                           as.character(paste0("Oncogenic role (CancerMine): ",citation_links_oncogene))),
                  citations_oncogene = dplyr::if_else(is.na(oncogene),as.character(NA),
                                                      as.character(citations_oncogene)),
                  citation_links_tsgene = dplyr::if_else(is.na(tumor_suppressor_cancermine),
                                                         "Tumor suppressive role (CancerMine): No records",
                                                         as.character(paste0("Tumor suppressor role (CancerMine): ",citation_links_tsgene))),
                  citations_tsgene = dplyr::if_else(is.na(tumor_suppressor),as.character(NA),as.character(citations_tsgene))) %>%
    dplyr::mutate(cancergene_support = stringr::str_replace(paste(citation_links_oncogene,citation_links_tsgene,ncg_link, sep=", "),"^, ","")) %>%
    dplyr::mutate(tumor_suppressor_evidence = paste0("NCG:",ncg_tumor_suppressor,"&CancerMine:",tumor_suppressor_cancermine,":",n_citations_tsgene)) %>%
    dplyr::mutate(oncogene_evidence = paste0("NCG:",ncg_oncogene,"&CancerMine:",oncogene_cancermine,":",n_citations_oncogene)) %>%
    dplyr::mutate(cancer_driver_evidence = paste0("CancerMine:",cancer_driver,":",n_citations_cdriver)) %>%
    dplyr::mutate(onco_ts_ratio = dplyr::if_else(tumor_suppressor == T & oncogene == T,
                                                 as.numeric(n_citations_oncogene)/n_citations_tsgene,
                                                 as.numeric(NA))) %>%
    dplyr::mutate(tumor_suppressor = dplyr::if_else(tumor_suppressor == T & !is.na(onco_ts_ratio) & onco_ts_ratio > 3,
                                                    FALSE,
                                                    as.logical(tumor_suppressor))) %>%
    dplyr::mutate(oncogene = dplyr::if_else(!is.na(onco_ts_ratio) & oncogene == T & onco_ts_ratio <= 0.33,
                                            FALSE,
                                            as.logical(oncogene)))


  n_onc_ts <- dplyr::filter(tsgene_full, oncogene == T & tumor_suppressor == T)
  n_ts <- dplyr::filter(tsgene_full, tumor_suppressor == T & oncogene == F)
  n_onc <- dplyr::filter(tsgene_full, oncogene == T & tumor_suppressor == F)

  rlogging::message("A total of n = ",nrow(n_onc)," classified proto-oncogenes were retrieved from CancerMine/NCG")
  rlogging::message("A total of n = ",nrow(n_ts)," classified tumor suppressor genes were retrieved from CancerMine/NCG")
  rlogging::message("A total of n = ",nrow(n_onc_ts)," genes were annotated with dual roles as tumor suppressor genes and oncogenes from CancerMine/NCG")

  return(tsgene_full)

}


get_curated_fp_cancer_genes <- function(basedir = NULL,
                                        gene_info = NULL){

  invisible(assertthat::assert_that(!is.null(gene_info), msg = "'gene_info' object is NULL"))
  invisible(assertthat::assert_that(dir.exists(basedir), msg = paste0("Directory '",basedir,"' does not exist")))
  assertable::assert_colnames(gene_info,c("symbol","entrezgene"), only_colnames = F, quiet = T)
  assertable::assert_coltypes(gene_info, list(symbol = character(), entrezgene = integer()), quiet = T)
  xlsx_fname <- paste0(basedir,"/data-raw/cancer_driver/bailey_2018_cell.xlsx")
  invisible(assertthat::assert_that(file.exists(xlsx_fname),msg = paste0("File ",xlsx_fname," does not exist")))
  tmp <- openxlsx::read.xlsx(xlsx_fname,sheet = 8,startRow = 4)
  fp_cancer_drivers <- data.frame('symbol' = tmp[,2], stringsAsFactors = F) %>%
    dplyr::mutate(fp_driver_gene = TRUE) %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::left_join(dplyr::select(gene_info, symbol, entrezgene), by = c("symbol"))
  return(fp_cancer_drivers)

}


get_opentarget_associations_v2 <-
  function(basedir = NULL,
           min_overall_score = 0.1,
           min_num_sources = 2,
           release = "2021_04",
           direct_associations_only = F){

    opentarget_targets <- as.data.frame(
      readRDS(
        paste0(basedir,
               "/data-raw/opentargets/opentargets_target_2021.04.rds")
      ) %>%
        dplyr::select(target_symbol,
                      tractability_antibody,
                      tractability_small_molecule) %>%
        dplyr::rename(symbol = target_symbol) %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(
          ot_tractability_compound = paste(unique(tractability_small_molecule), collapse = "&"),
          ot_tractability_antibody = paste(unique(tractability_antibody), collapse = "&"),
          .groups = "drop"))


    opentarget_associations <- as.data.frame(
      readRDS(
        paste0(basedir,
               "/data-raw/opentargets/opentargets_association_direct_HC_2021.04.rds")
      ) %>%
        dplyr::rename(disease_efo_id = efo_id, symbol = target_symbol) %>%
        dplyr::mutate(disease_efo_id = stringr::str_replace(disease_efo_id,":","_")) %>%
        dplyr::mutate(overall_datatype_harmonic_score =
                        round(overall_datatype_harmonic_score, digits = 7)) %>%
        dplyr::filter(overall_datatype_harmonic_score >= min_overall_score) %>%
        dplyr::mutate(association_key =
                        paste(disease_efo_id, "T",
                              overall_datatype_harmonic_score,
                              sep=":")) %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(
          ot_association = paste(association_key, collapse = "&"),
          .groups = "drop")
    )

    ot_associations <- opentarget_targets %>%
      dplyr::left_join(opentarget_associations, by = "symbol")

    return(ot_associations)
  }

get_opentarget_associations <- function(basedir = NULL,
                                        min_overall_score = 0.4,
                                        min_num_sources = 2,
                                        version = "2021_02",
                                        direct_associations_only = F){

  rlogging::message("Retrieving target-disease associations and target tractability evidence from OpenTargets platform (",version,")")
  opentarget <- as.data.frame(
    readr::read_tsv(file=paste0(basedir,"/data-raw/opentargets/opentargets_associations.tsv.gz"), progress = F) %>%
      dplyr::mutate(disease_efo_id = stringr::str_replace(disease_efo_id,"_",":"))
  )

  #efo_map <- readRDS(file=paste0(basedir,"/phenotype_ontology/efo_map.rds")) %>%
  #dplyr::select(efo_id, cui)

  if(direct_associations_only == T){
    opentarget <- opentarget %>% dplyr::filter(association_is_direct == T)
  }

  opentarget_associations <- as.data.frame(
    opentarget %>%
      dplyr::filter(association_overall >= min_overall_score) %>%
      dplyr::mutate(num_sources = 0) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_progeny != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_sysbio != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_expression_atlas != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_europepmc != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_uniprot_literature != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_phenodigm != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_eva != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_gene2phenotype != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_slapenrich != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_genomics_england != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_postgap != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_chembl != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_cancer_gene_census != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_reactome != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_eva_somatic != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_phewas_catalog != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_crispr != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_ot_genetics_portal != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_clingen != 0,num_sources + 1,num_sources)) %>%
      dplyr::mutate(num_sources = dplyr::if_else(association_intogen != 0,num_sources + 1,num_sources)) %>%
      dplyr::filter(num_sources >= min_num_sources) %>%
      dplyr::mutate(disease_efo_id = stringr::str_replace(disease_efo_id,":","_")) %>%
      dplyr::mutate(association_key =
                      paste(disease_efo_id, association_is_direct,
                            association_overall,sep=":")) %>%
      dplyr::rename(symbol = target_symbol) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(
        ot_association = paste(association_key, collapse = "&"),
        ot_tractability_compound = paste(unique(tractability_small_molecule), collapse = "&"),
        ot_tractability_antibody = paste(unique(tractability_antibody), collapse = "&"),
        .groups = "drop"))

  return(opentarget_associations)

}

get_dbnsfp_gene_annotations <- function(basedir = '/Users/sigven/research/DB/var_annotation_tracks'){

  rlogging::message("Retrieving gene damage scores/OMIM annotation from dbNSFP_gene")
  dbnsfp_gene <- read.table(file=gzfile(paste0(basedir,"/data-raw/dbnsfp/dbNSFP_gene.gz")),sep="\t",
                            header = T, stringsAsFactors = F, na.strings = c(".",""), comment.char="",
                            quote = NULL) %>%
    janitor::clean_names() %>%
    dplyr::select(gene_name, entrez_gene_id, function_description) %>%
    dplyr::rename(symbol = gene_name) %>%
    dplyr::mutate(function_description = stringr::str_replace_all(function_description,"HAMAP- Rule","HAMAP-Rule")) %>%
    dplyr::mutate(function_description = stringr::str_replace_all(function_description,"FUNCTION: |\\{ECO:[0-9]{1,}(\\|(PubMed|HAMAP-Rule|UniProtKB):\\S+){0,1}(, ECO:[0-9]{1,}(\\|(PubMed|UniProtKB|HAMAP-Rule):\\S+){0,}){0,}\\}\\.","")) %>%
    dplyr::mutate(function_description = stringr::str_replace_all(function_description,"( )?; $","")) %>%
    dplyr::mutate(function_description = stringr::str_replace_all(function_description,"( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)","")) %>%
    dplyr::mutate(function_description = stringr::str_replace_all(function_description,"( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)","")) %>%
    dplyr::select(symbol, function_description)


  return(dbnsfp_gene)
}

get_gencode_transcripts <- function(basedir = "/Users/sigven/research/DB/var_annotation_tracks",
                                            build = "grch38",
                                            append_regulatory_region = FALSE,
                                            gencode_version = "38",
                                            gene_info = NULL){
  gencode_ftp_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"


  rlogging::message(paste0("Retrieving GENCODE transcripts - version ",
                           gencode_version,", build ",build))
  if(build == 'grch37'){
    if(!file.exists(paste0(basedir,"/data-raw/gencode/",build,"/gencode.",
                           build,".annotation.gtf.gz"))){
      rlogging::message(paste0(gencode_ftp_url,
                               "release_19/gencode.v19.annotation.gtf.gz"))
      download.file(paste0(gencode_ftp_url,"release_19/gencode.v19.annotation.gtf.gz"),
                    destfile = paste0(basedir,"/data-raw/gencode/",
                                      build,"/gencode.",
                                      build,
                                      ".annotation.gtf.gz"),quiet = T)
    }
  }
  if(build == 'grch38'){
    if(!file.exists(paste0(basedir,"/data-raw/gencode/",
                           build,"/gencode.",
                           build,".annotation.gtf.gz"))){
      rlogging::message(
        paste0("Downloading ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
               gencode_version,
               "/gencode.v",
               gencode_version,".annotation.gtf.gz"))
      download.file(paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                           gencode_version,
                           "/gencode.v",
                           gencode_version,".annotation.gtf.gz"),
                    destfile = paste0(basedir,
                                      "/data-raw/gencode/",
                                      build,
                                      "/gencode.",
                                      build,".annotation.gtf.gz"),
                    quiet = T)
    }
  }
  gencode_gtf <- read.table(gzfile(paste0(basedir,"/data-raw/gencode/",build,
                                          "/gencode.",build,".annotation.gtf.gz")),
                            sep="\t", stringsAsFactors = F,na.strings=".", skip=5,
                            comment.char = "#", quote = "") %>%
    dplyr::filter(V3 == 'transcript')

  gencode_transcript_start <- gencode_gtf$V4
  gencode_transcript_end <- gencode_gtf$V5
  gencode_transcript_strand <- gencode_gtf$V7
  gencode_chrom <- gencode_gtf$V1
  gencode_gtf_annotation <- stringr::str_replace_all(gencode_gtf$V9,'"','')
  gencode_ensembl_gene_id <- stringr::str_match(gencode_gtf_annotation,"gene_id ENSG[0-9]{1,}")[,1]
  gencode_ensembl_gene_id <- stringr::str_replace(gencode_ensembl_gene_id,"gene_id ","")
  gencode_ensembl_trans_id <- stringr::str_match(gencode_gtf_annotation,"transcript_id ENST[0-9]{1,}\\.[0-9]{1,}")[,1]
  gencode_ensembl_trans_id <- stringr::str_replace(gencode_ensembl_trans_id,"transcript_id ","")
  gencode_ensembl_prot_id <- stringr::str_match(gencode_gtf_annotation,"protein_id ENSP[0-9]{1,}\\.[0-9]{1,}")[,1]
  gencode_ensembl_prot_id <- stringr::str_replace(gencode_ensembl_prot_id,"protein_id ","")
  gencode_ensembl_trans_id_original <- gencode_ensembl_trans_id
  gencode_ensembl_trans_id <- stringr::str_replace(gencode_ensembl_trans_id,"\\.[0-9]{1,}$","")
  gencode_gene_type <- stringr::str_match(gencode_gtf_annotation,"gene_type \\S+")[,1]
  gencode_gene_type <- stringr::str_replace_all(gencode_gene_type,"gene_type |;$","")
  gencode_trans_type <- stringr::str_match(gencode_gtf_annotation,"transcript_type \\S+")[,1]
  gencode_trans_type <- stringr::str_replace_all(gencode_trans_type,"transcript_type |;$","")
  gencode_tag <- stringr::str_match(gencode_gtf_annotation,"tag \\S+")[,1]
  gencode_tag <- stringr::str_replace_all(gencode_tag,"tag |;$","")
  gencode_genestatus <- stringr::str_match(gencode_gtf_annotation,"gene_status \\S+")[,1]
  gencode_genestatus <- stringr::str_replace_all(gencode_genestatus,"gene_status |;$","")
  gencode_symbol <- stringr::str_match(gencode_gtf_annotation,"gene_name \\S+")[,1]
  gencode_symbol <- stringr::str_replace_all(gencode_symbol,"gene_name |;$","")
  gencode <- data.frame('chrom' = gencode_chrom,
                        'transcript_start' = gencode_transcript_start,
                        'transcript_end' = gencode_transcript_end,
                        'strand' = gencode_transcript_strand,
                        'ensembl_gene_id' = gencode_ensembl_gene_id,
                        'symbol' = gencode_symbol,
                        'gencode_gene_biotype' = gencode_gene_type,
                        'gencode_genestatus' = gencode_genestatus,
                        'gencode_tag' = gencode_tag,
                        'gencode_transcript_type' = gencode_trans_type,
                        'ensembl_transcript_id' = gencode_ensembl_trans_id,
                        'ensembl_transcript_id_full' = gencode_ensembl_trans_id_original,
                        'ensembl_protein_id' = gencode_ensembl_prot_id,
                        stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ensembl_gene_id) &
                    !is.na(ensembl_transcript_id)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(gencode_release = gencode_version)

  rlogging::message(paste0("A total of ",nrow(gencode)," transcripts parsed"))

  ## include regulatory region (for VEP annotation)
  if(append_regulatory_region == T){
    rlogging::message("Parameter 'append_regulatory_region' is TRUE: expanding transcript start/end with 5kb (for VEP consequence compliance)")
    chromosome_lengths <- data.frame('chrom' = head(names(seqlengths(BSgenome.Hsapiens.UCSC.hg19)),25),
                                     'chrom_length' = head(seqlengths(BSgenome.Hsapiens.UCSC.hg19),25),
                                     stringsAsFactors = F, row.names = NULL)
    if(build == 'grch38'){
      chromosome_lengths <- data.frame('chrom' = head(names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)),25),
                                       'chrom_length' = head(seqlengths(BSgenome.Hsapiens.UCSC.hg38),25),
                                       stringsAsFactors = F, row.names = NULL)
    }
    gencode <- as.data.frame(
      gencode %>%
        dplyr::left_join(chromosome_lengths, by=c("chrom")) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(start = dplyr::if_else(chrom != "chrM",
                                             as.integer(
                                               max(1, transcript_start - 5001)
                                             ),
                                             as.integer(transcript_start)),
                      end = dplyr::if_else(chrom != "chrM",
                                           as.integer(
                                             min(chrom_length,transcript_end + 5001)
                                           ),
                                           as.integer(transcript_end))
        ) %>%
        dplyr::select(-chrom_length)
    )
  }

  ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")
  if(build == 'grch37'){
    ensembl_mart <- biomaRt::useMart(biomart = "ensembl")
  }
  ensembl_genes <- biomaRt::useDataset("hsapiens_gene_ensembl",
                                       mart = ensembl_mart)
  queryAttributes <- c('ensembl_transcript_id','refseq_mrna')
  ensembl_genes_xref1 <- as.data.frame(
    biomaRt::getBM(attributes = queryAttributes, mart = ensembl_genes) %>%
      dplyr::mutate(refseq_mrna =
                      dplyr::if_else(refseq_mrna == '',
                                     as.character(NA),
                                     as.character(refseq_mrna))) %>%
      dplyr::filter(!is.na(refseq_mrna)) %>%
      dplyr::group_by(ensembl_transcript_id) %>%
      dplyr::summarise(refseq_mrna = paste(unique(refseq_mrna),
                                           collapse = "&"),
                       .groups = "drop") %>%
      dplyr::ungroup()
  )




  queryAttributes <- c('ensembl_gene_id','hgnc_symbol','hgnc_id')
  ensembl_genes_xref3 <- as.data.frame(
    biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes) %>%
      dplyr::mutate(hgnc_symbol = dplyr::if_else(hgnc_symbol == '',
                                                 as.character(NA),
                                                 as.character(hgnc_symbol))) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(hgnc_symbol = paste(unique(hgnc_symbol),
                                           collapse="&"),
                       hgnc_id = paste(unique(hgnc_id),
                                       collapse="&"),
                       .groups = "drop") %>%
      dplyr::ungroup()
  )

  queryAttributes <- c('hgnc_symbol','entrezgene_id')
  ensembl_genes_xref4 <- as.data.frame(
    biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes) %>%
      dplyr::mutate(entrezgene = as.integer(entrezgene_id)) %>%
      dplyr::filter(hgnc_symbol != "") %>%
      dplyr::select(-entrezgene_id) %>%
      dplyr::filter(!is.na(entrezgene))
  )

  tmp1 <- as.data.frame(ensembl_genes_xref4 %>%
                          dplyr::group_by(hgnc_symbol) %>%
                          dplyr::summarise(n = dplyr::n(), .groups = "drop"))

  b_nonduplicates <- ensembl_genes_xref4 %>%
    dplyr::left_join(tmp1,by=c("hgnc_symbol")) %>%
    dplyr::filter(n == 1)
  b_duplicates <- ensembl_genes_xref4 %>%
    dplyr::left_join(tmp1,by=c("hgnc_symbol")) %>%
    dplyr::filter(n > 1)
  b <- as.data.frame(mygene::queryMany(unique(b_duplicates$hgnc_symbol),
                                       scopes="symbol",
                                       fields="entrezgene",
                                       species="human")) %>%
    dplyr::rename(hgnc_symbol = query) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene)) %>%
    dplyr::select(hgnc_symbol, entrezgene)

  b_duplicates <- b_duplicates %>% dplyr::inner_join(b,by=c("hgnc_symbol","entrezgene"))
  ensembl_genes_xref4 <- dplyr::bind_rows(b_nonduplicates, b_duplicates) %>%
    dplyr::select(-n)

  queryAttributes <- c('entrezgene_id','entrezgene_description')
  ensembl_genes_xref5 <- biomaRt::getBM(attributes = queryAttributes, mart=ensembl_genes) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene_id), description = entrezgene_description) %>%
    dplyr::mutate(description = dplyr::if_else(description == '',as.character(NA),as.character(description))) %>%
    dplyr::select(-c(entrezgene_id, entrezgene_description)) %>%
    dplyr::mutate(description = stringr::str_replace(description," \\[.+$","")) %>%
    dplyr::filter(!is.na(entrezgene))

  queryAttributes <- c('ensembl_transcript_id','refseq_peptide')
  ensembl_genes_xref6 <- as.data.frame(
    biomaRt::getBM(attributes = queryAttributes, mart = ensembl_genes) %>%
      dplyr::mutate(refseq_peptide =
                      dplyr::if_else(refseq_peptide == '',
                                     as.character(NA),
                                     as.character(refseq_peptide))) %>%
      dplyr::filter(!is.na(refseq_peptide)) %>%
      dplyr::group_by(ensembl_transcript_id) %>%
      dplyr::summarise(refseq_peptide = paste(unique(refseq_peptide),
                                              collapse = "&"),
                       .groups = "drop") %>%
      dplyr::ungroup()
  )


  ensembl_genes_xref2 <- as.data.frame(
    ensembl_genes_xref3 %>%
      dplyr::left_join(ensembl_genes_xref4, by=c("hgnc_symbol")) %>%
      dplyr::left_join(ensembl_genes_xref5, by=c("entrezgene")) %>%
      dplyr::rename(symbol = hgnc_symbol, name = description)
  )
  ensembl_genes_xref2$entrezgene <- as.integer(ensembl_genes_xref2$entrezgene)

  gencode <- gencode %>%
    dplyr::left_join(dplyr::select(gene_info, hgnc_id,
                                   entrezgene,name, symbol),
                     by=c("symbol" = "symbol")) %>%
    dplyr::left_join(ensembl_genes_xref1, by=c("ensembl_transcript_id")) %>%
    dplyr::left_join(ensembl_genes_xref6, by=c("ensembl_transcript_id"))

  gencode_not_missing <- gencode %>%
    dplyr::filter(!is.na(symbol) & !is.na(entrezgene) & !is.na(name))
  gencode_missing_name_entrez <- gencode %>%
    dplyr::filter(!is.na(symbol) & is.na(entrezgene)) %>%
    dplyr::select(-c(name,entrezgene)) %>%
    dplyr::left_join(dplyr::select(ensembl_genes_xref2,
                                   entrezgene,
                                   ensembl_gene_id),by=c("ensembl_gene_id")) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::left_join(dplyr::select(gene_info, entrezgene, name),
                     by = c("entrezgene")) %>%
    dplyr::filter(!is.na(name))

  gencode_found <- dplyr::bind_rows(gencode_not_missing, gencode_missing_name_entrez)
  gencode_remain <- dplyr::anti_join(gencode, gencode_found, by=c("ensembl_gene_id"))
  gencode_remain <- gencode_remain %>%
    dplyr::filter(!is.na(symbol) & is.na(entrezgene)) %>%
    dplyr::select(-c(name,entrezgene, hgnc_id)) %>%
    dplyr::left_join(dplyr::select(ensembl_genes_xref2,
                                   entrezgene, name,
                                   ensembl_gene_id,
                                   hgnc_id),by=c("ensembl_gene_id"))


  gencode_remain$hgnc_id <- as.character(gencode_remain$hgnc_id)
  gencode <- dplyr::bind_rows(gencode_found, gencode_remain)
  #gencode <- gencode %>% dplyr::mutate(symbol_entrez = symbol)

  gencode <- gencode %>%
    dplyr::filter(!is.na(ensembl_transcript_id)) %>%
    dplyr::distinct()
  if(nrow(gencode[!is.na(gencode$refseq_mrna) & gencode$refseq_mrna == "NA",]) > 0){
    gencode[!is.na(gencode$refseq_mrna) & gencode$refseq_mrna == "NA",]$refseq_mrna <- NA
  }

  tmp <- gene_info %>%
    dplyr::select(symbol, name) %>%
    dplyr::rename(name_entrez = name)
  gencode <- gencode %>%
    dplyr::left_join(tmp, by=c("symbol" = "symbol")) %>%
    dplyr::mutate(name = dplyr::if_else(
      stringr::str_detect(name,"^(U|u)ncharacterized") &
        name != name_entrez,name_entrez,as.character(name))) %>%
    dplyr::select(-name_entrez) %>%
    dplyr::distinct()


  gencode_n <- gencode %>% dplyr::group_by(ensembl_transcript_id) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  gencode <- gencode %>% dplyr::left_join(gencode_n,by=c("ensembl_transcript_id"))

  gencode_1 <- dplyr::filter(gencode, n == 1) %>%
    dplyr::select(-n)
  gencode_2 <- dplyr::filter(gencode, n > 1) %>%
    dplyr::filter(!is.na(hgnc_id)) %>%
    dplyr::distinct() %>%
    dplyr::select(-n)

  gencode <- dplyr::bind_rows(gencode_1, gencode_2)

  gencode_1 <- gencode %>%
    dplyr::filter(is.na(entrezgene) & is.na(name)) %>%
    dplyr::select(-c(entrezgene,name))
  gene_info_name <- dplyr::select(gene_info, ensembl_gene_id,
                                  entrezgene, name) %>%
    dplyr::filter(!is.na(entrezgene))
  gencode_1 <- gencode_1 %>%
    dplyr::left_join(gene_info_name,by=c("ensembl_gene_id"))

  gencode_2 <- gencode %>%
    dplyr::anti_join(dplyr::select(gencode_1, ensembl_transcript_id),
                     by = c("ensembl_transcript_id"))
  gencode <- dplyr::bind_rows(gencode_1, gencode_2)


  rlogging::message(paste0("A total of ",nrow(gencode)," valid transcripts remaining"))

  return(gencode)

}

resolve_duplicates <- function(gencode_df){

  ## duplicates
  duplicate_symbols <- as.data.frame(
    gencode_df %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n > 1)
  )
  resolved_duplicates <- data.frame()

  i <- 1
  while(i <= nrow(duplicate_symbols)){
    d <- duplicate_symbols[i,]$symbol
    entries <- dplyr::filter(gencode_df, symbol == d)
    entry_1 <- as.integer(stringr::str_replace(
      entries[1,"ensembl_gene_id"],"ENSG(0{1,})",""))
    entry_2 <- as.integer(stringr::str_replace(
      entries[2,"ensembl_gene_id"],"ENSG(0{1,})",""))
    if(entry_1 < entry_2){
      resolved_duplicates <-
        dplyr::bind_rows(resolved_duplicates, entries[1,])
    }else{
      resolved_duplicates <-
        dplyr::bind_rows(resolved_duplicates, entries[2,])
    }
    i <- i + 1
  }

  gencode_df <- gencode_df %>%
    dplyr::anti_join(duplicate_symbols, by = c("symbol")) %>%
    dplyr::bind_rows(resolved_duplicates) %>%
    dplyr::arrange(desc(symbol))

  return(gencode_df)

}

get_uniprot_map <- function(basedir = NULL,
                            uniprot_release = "2020_04"){

  ## corum protein-complex data

  #download.file("http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip",destfile="data-raw/chorum/coreComplexes.txt.zip", quiet = T)
  corum_complexes <-
    read.table(file=paste0(basedir,"/data-raw/corum/coreComplexes.txt"),
               sep="\t",header=T, quote="",comment.char="",stringsAsFactors = F) %>%
    dplyr::select(c(Organism, ComplexID, subunits.UniProt.IDs.,
                    ComplexName,FunCat.description)) %>%
    dplyr::filter(Organism == 'Human') %>%
    dplyr::mutate(n_complexes = nrow(.)) %>%
    dplyr::rename(uniprot_acc = subunits.UniProt.IDs.,
                  corum_id = ComplexID, corum_name = ComplexName,
                  corum_funcat_description = FunCat.description) %>%
    tidyr::separate_rows(uniprot_acc, sep=";") %>%
    dplyr::select(-Organism)

  n_complexes <- unique(corum_complexes$n_complexes)
  n_protein_complex_interactions <- nrow(corum_complexes)
  rlogging::message("Retrieving protein complex annotation from CORUM")
  rlogging::message("A total of ", n_complexes,
                    " protein complexes were parsed, containing ",
                    n_protein_complex_interactions," protein-complex interactions")

  uniprot_acc_to_corum_ids <- as.data.frame(
    dplyr::group_by(corum_complexes, uniprot_acc) %>%
      dplyr::summarise(corum_id = paste(unique(corum_id), collapse = "&"), .groups = "drop"))
  rlogging::message("Retrieving UniProtKB annotation")

  ## read uniprot ID mapping
  idmapping_fname <- paste0(basedir,'/data-raw/uniprot/',uniprot_release,'/HUMAN_9606_idmapping.dat.gz')
  idmapping_up_kb <- read.table(gzfile(idmapping_fname),sep="\t",header = F,quote = "", stringsAsFactors = F) %>%
    dplyr::filter(V2 == 'Ensembl_TRS' | V2 == 'UniProtKB-ID' | V2 == 'Ensembl' | V2 == 'GeneID' | V2 == 'RefSeq_NT') %>%
    magrittr::set_colnames(c('acc','type','name')) %>%
    dplyr::mutate(acc = stringr::str_replace(acc,"-[0-9]{1,}",""))

  reviewed_uniprot_accs <- readRDS(file=paste0(basedir,"/data-raw/uniprot/",uniprot_release,"/reviewed_uniprot_accessions.rds"))

  ## UniProt accession to Ensembl transcript ID
  ensembl_up_acc <- dplyr::filter(idmapping_up_kb, type == 'Ensembl_TRS') %>%
    dplyr::select(acc, name) %>%
    dplyr::rename(uniprot_acc = acc, ensembl_transcript_id = name) %>%
    dplyr::distinct()
  ## remove transcripts mapped to more than one UniProt acc
  unique_ensembl_trans <- as.data.frame(dplyr::group_by(ensembl_up_acc, ensembl_transcript_id) %>%
                                          dplyr::summarise(n_trans = dplyr::n(), .groups = "drop") %>%
                                          dplyr::filter(n_trans == 1))
  ensembl_up_acc <- dplyr::semi_join(ensembl_up_acc,
                                     dplyr::select(unique_ensembl_trans, ensembl_transcript_id),
                                     by=c("ensembl_transcript_id"))

  ## UniProt accession to UniProt ID
  up_id_acc <- as.data.frame(dplyr::filter(idmapping_up_kb, type == 'UniProtKB-ID') %>%
                               dplyr::group_by(acc) %>%
                               dplyr::summarise(uniprot_id = paste(unique(name), collapse="&"),
                                                .groups = "drop") %>%
                               dplyr::rename(uniprot_acc = acc))
  rlogging::message("A total of ",nrow(up_id_acc)," UniProt protein accessions were parsed")


  ## UniProt accession to UniProt ID
  refseq_id_acc <- as.data.frame(dplyr::filter(idmapping_up_kb, type == 'RefSeq_NT') %>%
                                   dplyr::mutate(name = stringr::str_replace(name,"\\.[0-9]{1,}$","")) %>%
                                   dplyr::group_by(acc) %>%
                                   dplyr::summarise(refseq_mrna = paste(unique(name), collapse="&"),
                                                    .groups = "drop") %>%
                                   dplyr::rename(uniprot_acc = acc))


  uniprot_map <- dplyr::full_join(dplyr::left_join(refseq_id_acc, up_id_acc,by=c("uniprot_acc")),
                                  dplyr::left_join(ensembl_up_acc, up_id_acc,by=c("uniprot_acc")),by=c("uniprot_acc","uniprot_id"))

  ## UniProt accession and ID for Ensembl transcript IDs
  uniprot_map <- uniprot_map %>%
    dplyr::left_join(uniprot_acc_to_corum_ids,by=c("uniprot_acc")) %>%
    dplyr::left_join(reviewed_uniprot_accs, by=c("uniprot_acc","uniprot_id")) %>%
    dplyr::rename(uniprot_reviewed = reviewed)

  return(list('uniprot_map' = uniprot_map, 'corum_complexes' = corum_complexes))

}

map_uniprot_accession <- function(gencode, uniprot_map){

  gencode <- as.data.frame(
    gencode %>%
      dplyr::left_join(dplyr::select(dplyr::filter(uniprot_map$uniprot_map,!is.na(ensembl_transcript_id)), -refseq_mrna),
                       by=c("ensembl_transcript_id"))
  )

  tmp0 <- dplyr::select(uniprot_map$uniprot_map, refseq_mrna, uniprot_id, uniprot_acc, uniprot_reviewed) %>%
    tidyr::separate_rows(refseq_mrna, sep="&") %>%
    dplyr::filter(!is.na(refseq_mrna)) %>%
    dplyr::rename(uniprot_id2 = uniprot_id, uniprot_acc2 = uniprot_acc, uniprot_reviewed2 = uniprot_reviewed) %>%
    dplyr::distinct()

  tmp1 <- as.data.frame(
    dplyr::select(gencode, ensembl_transcript_id, refseq_mrna,
                  uniprot_id, uniprot_reviewed, uniprot_acc) %>%
      tidyr::separate_rows(refseq_mrna, sep="&") %>%
      dplyr::filter(!is.na(refseq_mrna)) %>%
      dplyr::left_join(tmp0,by=c("refseq_mrna")) %>%
      dplyr::mutate(uniprot_id = dplyr::if_else((is.na(uniprot_id) | uniprot_id != uniprot_id2 ) & !is.na(uniprot_id2),
                                                uniprot_id2,
                                                as.character(uniprot_id))) %>%
      dplyr::mutate(uniprot_acc = dplyr::if_else((is.na(uniprot_acc) | uniprot_acc != uniprot_acc2 ) & !is.na(uniprot_acc2),
                                                 uniprot_acc2,
                                                 as.character(uniprot_acc))) %>%
      dplyr::mutate(uniprot_reviewed = dplyr::if_else((is.na(uniprot_reviewed) | uniprot_reviewed != uniprot_reviewed2 ) & !is.na(uniprot_reviewed2),as.logical(uniprot_reviewed2),as.logical(uniprot_reviewed))) %>%
      dplyr::filter(uniprot_reviewed == T) %>%
      dplyr::group_by(ensembl_transcript_id) %>%
      dplyr::summarise(refseq_mrna = paste(unique(refseq_mrna),collapse="&"),
                       uniprot_id = paste(unique(uniprot_id), collapse="&"),
                       uniprot_reviewed = paste(unique(uniprot_reviewed),collapse="&"),
                       uniprot_acc = paste(unique(uniprot_acc),collapse="&"), .groups = "drop") %>%
      dplyr::mutate(uniprot_reviewed = as.logical(uniprot_reviewed))
  )

  tmp2 <- dplyr::select(gencode, ensembl_transcript_id, refseq_mrna, uniprot_id, uniprot_acc, uniprot_reviewed) %>%
    dplyr::filter(is.na(refseq_mrna))

  uniprot_mapping_final <- dplyr::bind_rows(tmp1, tmp2)

  gencode <- dplyr::select(gencode, -c(uniprot_id,uniprot_acc,uniprot_reviewed,refseq_mrna)) %>%
    dplyr::left_join(uniprot_mapping_final,by=c("ensembl_transcript_id"))

  rm(tmp0)
  rm(tmp1)
  rm(tmp2)

  return(gencode)
}

get_cancer_hallmarks <- function(basedir = NULL,
                                 gene_info = NULL){

  hallmark_data_long <- as.data.frame(
    readRDS(file=paste0(basedir,"/data-raw/opentargets/opentargets_target_2021.04.rds")) %>%
      dplyr::select(target_symbol, cancer_hallmark) %>%
      dplyr::filter(!is.na(cancer_hallmark)) %>%
      tidyr::separate_rows(cancer_hallmark, sep="@@@") %>%
      tidyr::separate(cancer_hallmark,
                      into = c('hallmark','promote','suppress','literature_support'),
                      remove = F, sep = '\\|') %>%
      dplyr::rename(promotes = promote, suppresses = suppress) %>%
      dplyr::mutate(hallmark = stringr::str_to_title(hallmark)) %>%
      tidyr::separate_rows(literature_support, sep="&") %>%
      dplyr::mutate(promotes = as.logical(stringr::str_replace(promotes,"PROMOTE:",""))) %>%
      dplyr::mutate(suppresses = as.logical(stringr::str_replace(suppresses,"SUPPRESS:",""))) %>%
      tidyr::separate(literature_support, into = c('pmid','pmid_summary'), sep="%%%") %>%
      dplyr::mutate(pmid_summary = R.utils::capitalize(pmid_summary))
  )

  pmid_data <- get_citations_pubmed(unique(hallmark_data_long$pmid)) %>%
    dplyr::mutate(pmid = as.character(pmid))


  hallmark_data_long <- hallmark_data_long %>%
    dplyr::left_join(pmid_data, by = "pmid")

  hallmark_data_short <- as.data.frame(
    hallmark_data_long %>%
      dplyr::select(target_symbol, hallmark,
                    promotes, suppresses, pmid_summary,
                    link) %>%
      dplyr::mutate(literature_support = paste0(pmid_summary," (",link,")")) %>%
      dplyr::group_by(target_symbol, hallmark,
                      promotes, suppresses) %>%
      dplyr::summarise(literature_support = paste(
        literature_support, collapse=". "),
        .groups = "drop") %>%
      dplyr::left_join(dplyr::select(gene_info,
                                     symbol, entrezgene),
                       by = c("target_symbol" = "symbol")) %>%
      dplyr::mutate(symbol =
        paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                 entrezgene,"' target='_blank'>",target_symbol,"</a>")) %>%
      dplyr::select(-entrezgene) %>%
      dplyr::select(target_symbol, symbol, hallmark, promotes, suppresses,
                    literature_support)
  )

  return(list('long' = hallmark_data_long,
              'short' = hallmark_data_short))

}


