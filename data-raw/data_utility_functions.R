
#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed <- function(pmid){

  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmid, ceiling(length(pmid)/100))
  j <- 0
  all_citations <- data.frame()
  cat('Retrieving PubMed citations for PMID list, total length', length(pmid))
  cat('\n')
  while(j < length(pmid_chunks)){
    #cat(unlist(pmid_chunks),"\n")
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    cat(unlist(pmid_chunk),"\n")
    cat('Processing chunk ',j,' with ',length(pmid_chunk),'PMIDS')
    cat('\n')
    pmid_string <- paste(pmid_chunk,collapse = " ")
    res <- RISmed::EUtilsGet(
      RISmed::EUtilsSummary(pmid_string, type="esearch", db="pubmed", retmax = 5000)
    )


    year <- RISmed::YearPubmed(res)
    authorlist <- RISmed::Author(res)
    pmid_list <- RISmed::PMID(res)
    i <- 1
    first_author <- c()
    while(i <= length(authorlist)){

      first_author <- c(first_author,paste(authorlist[[i]][1,]$LastName," et al.",sep=""))
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(res)
    citations <- data.frame('pmid' = as.integer(pmid_list),
                            'citation' = paste(first_author,year,journal,sep=", "),
                            stringsAsFactors = F)
    citations$link <- paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',
                             citations$pmid,'\' target=\'_blank\'>',
                             citations$citation,'</a>')
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
    # dplyr::mutate(ensembl_gene_id = stringr::str_replace(
    #   stringr::str_match(V6,"Ensembl:ENSG[0-9]{1,}"), "Ensembl:", "")) %>%
    dplyr::mutate(hgnc_id = stringr::str_replace(
      stringr::str_match(V6,"HGNC:HGNC:[0-9]{1,}"), "HGNC:HGNC:", "")) %>%
    dplyr::select(-V6) %>%
    dplyr::mutate(symbol_entrez = dplyr::if_else(
      is.na(symbol_entrez) | symbol_entrez == "-",
      symbol, as.character(symbol_entrez))) %>%
    dplyr::mutate(gene_biotype = dplyr::if_else(
      gene_biotype ==  "protin-coding",
      "protein_coding", as.character(gene_biotype))) %>%
    dplyr::filter(nchar(symbol_entrez) > 0) %>%
    ## TEC
    dplyr::filter(entrezgene != 100124696) %>%
    ## HBD
    dplyr::filter(entrezgene != 100187828) %>%
    ## MMD2
    dplyr::filter(entrezgene != 100505381) %>%
    ## MEMO1
    dplyr::filter(entrezgene != 7795) %>%
    dplyr::distinct()
    #dplyr::select(-symbol) %>%
    #dplyr::rename(symbol = symbol_entrez)

  return(gene_info)

}

get_gene_aliases <- function(gene_info = NULL){

  invisible(assertthat::assert_that(!is.null(gene_info)))
  assertable::assert_colnames(
    gene_info, colnames = c("synonyms","symbol"),
    only_colnames = F, quiet = T)

  unique_synonyms <- as.data.frame(
    gene_info %>%
      dplyr::select(synonyms, entrezgene) %>%
      tidyr::separate_rows(synonyms, sep="\\|") %>%
      dplyr::group_by(synonyms) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n == 1)) %>%
    dplyr::select(-n)

  alias2primary <- as.data.frame(
    gene_info %>%
      dplyr::select(synonyms, entrezgene) %>%
      tidyr::separate_rows(synonyms,sep="\\|") %>%
      dplyr::inner_join(unique_synonyms, by = "synonyms") %>%
      dplyr::rename(value = synonyms) %>%
      dplyr::mutate(property = "alias") %>%
      dplyr::select(entrezgene, property, value) %>%
      dplyr::distinct()
  )

  return(alias2primary)
}

get_msigdb_signatures <- function(basedir = NULL,
                                  db_version = 'v7.5.1 (January 2022)'){

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
                                        version = "39",
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
      ncg_tumor_suppressor = dplyr::if_else(
        stringr::str_detect(cgc_annotation,"TSG") | ncg_tsg == 1,
        as.logical(TRUE),as.logical(FALSE))) %>%
    dplyr::mutate(
      ncg_oncogene = dplyr::if_else(
        stringr::str_detect(cgc_annotation,"oncogene") | ncg_oncogene == 1,
        as.logical(TRUE),as.logical(FALSE))) %>%
    dplyr::select(entrez, ncg_tumor_suppressor, ncg_oncogene) %>%
    dplyr::rename(entrezgene = entrez) %>%
    dplyr::filter(!(ncg_tumor_suppressor == F & ncg_oncogene == F)) %>%
    dplyr::anti_join(dplyr::select(fp_cancer_drivers, entrezgene), by = "entrezgene") %>%
    dplyr::distinct()

  pmids <- as.data.frame(
    read.table(gzfile(paste0(basedir,'/data-raw/cancermine/cancermine_sentences.tsv.gz')),
               header = T, comment.char = "", quote = "", sep="\t", stringsAsFactors = F) %>%
      dplyr::filter(predictprob >= 0.9) %>%
      dplyr::group_by(role, gene_entrez_id, pmid) %>%
      dplyr::summarise(doid = paste(unique(cancer_id), collapse=","), .groups = "drop") %>%
      dplyr::rename(entrezgene = gene_entrez_id) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, entrezgene), by = "entrezgene") %>%
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
      dplyr::arrange(entrezgene, desc(pmid)) %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(
        pmids_oncogene = paste(pmid,collapse=", "),
        citations_oncogene = paste(citation,collapse="; "),
        citation_links_oncogene = paste(head(citation_link,50),collapse=", "),
        .groups = "drop") %>%
      dplyr::mutate(n_citations_oncogene = as.integer(stringr::str_count(pmids_oncogene,", ")) + 1) %>%
      dplyr::filter(n_citations_oncogene > 1) %>%
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene > 4,
        "MC",
        as.character("LC"))) %>%
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene >= 10,
        "HC",
        as.character(oncogene_cancermine)))
  )

  pmids_tsgene <- as.data.frame(
    pmids %>%
      dplyr::filter(role == "Tumor_Suppressor") %>%
      dplyr::arrange(entrezgene, desc(pmid)) %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(
        pmids_tsgene = paste(unique(pmid),collapse=", "),
        citations_tsgene = paste(citation,collapse="; "),
        citation_links_tsgene = paste(head(citation_link,50),collapse=", "),
        .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_citations_tsgene = as.integer(stringr::str_count(pmids_tsgene,", ")) + 1) %>%
      dplyr::filter(n_citations_tsgene > 1) %>%
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene > 4,
        "MC",
        as.character("LC"))) %>%
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene >= 10,
        "HC",
        as.character(tumor_suppressor_cancermine)))
  )

  pmids_cdriver <- as.data.frame(
    pmids %>%
      dplyr::filter(role == "Driver") %>%
      dplyr::arrange(entrezgene, desc(pmid)) %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(
        pmids_cdriver = paste(pmid,collapse=", "),
        citations_cdriver = paste(citation,collapse="; "),
        citation_links_cdriver = paste(head(citation_link,50),collapse=", "),
        .groups = "drop") %>%
      dplyr::mutate(n_citations_cdriver = as.integer(stringr::str_count(pmids_cdriver,", ")) + 1) %>%
      dplyr::filter(n_citations_cdriver > 1) %>%
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver > 4,
        "MC",
        as.character("LC"))) %>%
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver >= 10,
        "HC",
        as.character(cancer_driver_cancermine)))
  )


  rlogging::message("Retrieving known proto-oncogenes/tumor suppressor genes from CancerMine")
  oncogene <- as.data.frame(
    readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                    col_names = T,na ="-", comment = "#", quote = "",
                    show_col_types = F) %>%
      dplyr::filter(role == "Oncogene") %>%
      dplyr::rename(entrezgene = gene_entrez_id) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, entrezgene), by = "entrezgene") %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(doid_oncogene = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
      dplyr::inner_join(pmids_oncogene, by = "entrezgene") %>%
      dplyr::distinct())

  n_hc_oncogene <- oncogene %>%
    dplyr::filter(oncogene_cancermine == "HC") %>%
    nrow()


  tsgene <- as.data.frame(
    readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                    col_names = T,na ="-", comment = "#", quote = "",
                    show_col_types = F) %>%
      dplyr::filter(role == "Tumor_Suppressor") %>%
      dplyr::rename(entrezgene = gene_entrez_id) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, entrezgene), by = "entrezgene") %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(doid_tsgene = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
      dplyr::inner_join(pmids_tsgene, by = "entrezgene") %>%
      dplyr::distinct()
  )
  n_hc_tsgene <- tsgene %>%
    dplyr::filter(tumor_suppressor_cancermine == "HC") %>%
    nrow()


  cdriver <- as.data.frame(
    readr::read_tsv(paste0(basedir,'/data-raw/cancermine/cancermine_collated.tsv.gz'),
                    col_names = T,na ="-", comment = "#", quote = "",
                    show_col_types = F) %>%
      dplyr::filter(role == "Driver") %>%
      dplyr::rename(entrezgene = gene_entrez_id) %>%
      dplyr::anti_join(dplyr::select(fp_cancer_drivers, entrezgene), by = "entrezgene") %>%
      dplyr::group_by(entrezgene) %>%
      dplyr::summarise(doid_cdriver = paste(unique(cancer_id),collapse=","), .groups = "drop") %>%
      dplyr::inner_join(pmids_cdriver, by = "entrezgene") %>%
      dplyr::distinct()
  )
  n_hc_cdriver <- cdriver %>%
    dplyr::filter(cancer_driver_cancermine == "HC") %>%
    nrow()

  tsgene_full <- dplyr::full_join(tsgene, oncogene, by = c("entrezgene")) %>%
    dplyr::full_join(cdriver, by="entrezgene") %>%
    dplyr::full_join(ncg, by = "entrezgene") %>%
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
                                                    paste0("NCG: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",entrezgene,"\" target=\"_blank\">Tumor Suppressor</a>"),
                                                    as.character("NCG: Not defined"))) %>%
    dplyr::mutate(ncg_links_tsgene = dplyr::if_else(is.na(ncg_links_tsgene),
                                                    as.character("NCG: Not defined"),
                                                    as.character(ncg_links_tsgene))) %>%
    dplyr::mutate(ncg_links_oncogene = dplyr::if_else(ncg_oncogene == T,
                                                      paste0("NCG: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",entrezgene,"\" target=\"_blank\">Oncogene</a>"),
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
                                            as.logical(oncogene))) %>%
    dplyr::mutate(
      citation_links_oncogene =
        stringr::str_replace(
          citation_links_oncogene,"(NA, ){1,}","")) %>%
    dplyr::mutate(
      citation_links_tsgene =
        stringr::str_replace(
          citation_links_tsgene,"(NA, ){1,}",""))


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
    dplyr::inner_join(dplyr::select(gene_info, symbol, entrezgene), by = c("symbol")) %>%
    dplyr::select(entrezgene, fp_driver_gene)
  return(fp_cancer_drivers)

}

get_opentarget_associations <-
  function(basedir = '/Users/sigven/research/DB/var_annotation_tracks',
           min_overall_score = 0.1,
           min_num_sources = 2,
           release = "2022.02",
           direct_associations_only = F){

    opentarget_targets <- as.data.frame(
      readRDS(
        paste0(basedir,
               "/data-raw/opentargets/opentargets_target_",
               release,
               ".rds")
      ) %>%
        dplyr::select(target_ensembl_gene_id,
                      SM_tractability_category,
                      SM_tractability_support,
                      AB_tractability_category,
                      AB_tractability_support) %>%
        dplyr::rename(ensembl_gene_id = target_ensembl_gene_id)
        # dplyr::mutate(target_tractability_symbol = paste0(
        #   "<a href='https://platform.opentargets.org/target/", target_ensembl_gene_id,"' target='_blank'>",symbol,"</a>")
        # )
      )



    opentarget_associations <- as.data.frame(
      readRDS(
        paste0(basedir,
               "/data-raw/opentargets/opentargets_association_direct_HC_",
               release,
               ".rds")
      ) %>%
        dplyr::rename(disease_efo_id = disease_id, ensembl_gene_id = target_ensembl_gene_id) %>%
        dplyr::mutate(disease_efo_id = stringr::str_replace(disease_efo_id,":","_")) %>%
        dplyr::filter(score >= min_overall_score) %>%
        dplyr::filter(stringr::str_count(datatype_items,",") >= min_num_sources - 1) %>%
        dplyr::mutate(association_key =
                        paste(disease_efo_id, "T",
                              score,
                              sep=":")) %>%
        dplyr::group_by(ensembl_gene_id) %>%
        dplyr::summarise(
          ot_association = paste(association_key, collapse = "&"),
          .groups = "drop")
    )

    ot_associations <- opentarget_targets %>%
      dplyr::left_join(opentarget_associations, by = "ensembl_gene_id")


    return(ot_associations)

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
    dplyr::select(symbol, function_description) %>%
    dplyr::mutate(
      gene_summary_uniprot =
        dplyr::if_else(
          !is.na(function_description),
          paste0("<b>UniProt:</b> ",
                 function_description),
          as.character(function_description)
        )) %>%
    dplyr::select(symbol, gene_summary_uniprot)


  return(dbnsfp_gene)
}

get_gencode_ensembl_transcripts <-
  function(basedir = NULL,
           build = "grch38",
           gencode_version = "39",
           gene_info = NULL,
           update = F){

    gencode_rds <- file.path(
      basedir, "data-raw","gencode", build,
      paste0("gencode_transcripts_", gencode_version, ".rds"))


    if(update == F & file.exists(gencode_rds)){
      rlogging::message(
        paste0("Loading pre-processed GENCODE transcripts - version ",
               gencode_version,", build ",build))
      gencode_transcripts <-
        readRDS(file = gencode_rds)
      return(gencode_transcripts)
    }

    gencode_ftp_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    rlogging::message(paste0("Retrieving GENCODE transcripts - version ",
                             gencode_version,", build ",build))

    destfile_gtf <-
      file.path(basedir, "data-raw","gencode", build,
                paste0("gencode.", build, ".annotation.v",
                       gencode_version, ".gtf.gz"))

    remote_gtf <-
      paste0(gencode_ftp_url,"release_", gencode_version, "/",
             "gencode.v", gencode_version,".annotation.gtf.gz")

    if(!file.exists(destfile_gtf)){
      rlogging::message(paste0("Downloading ",remote_gtf))
      download.file(remote_gtf,
                    destfile = destfile_gtf,
                    quiet = T)
    }

    gencode_gtf <- valr::read_gtf(destfile_gtf) %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::rename(transcript_start = start,
                    transcript_end = end,
                    ensembl_gene_id = gene_id,
                    symbol = gene_name,
                    gencode_gene_biotype = gene_type,
                    gencode_transcript_type = transcript_type,
                    ensembl_transcript_id = transcript_id,
                    ensembl_protein_id = protein_id) %>%
      dplyr::select(chrom,
                    transcript_start,
                    transcript_end,
                    strand,
                    symbol,
                    ensembl_gene_id,
                    ensembl_transcript_id,
                    gencode_gene_biotype,
                    gencode_transcript_type,
                    ensembl_protein_id) %>%
      dplyr::mutate(ensembl_transcript_id =
                      stringr::str_replace(ensembl_transcript_id,
                                           "\\.[0-9]{1,}$","")) %>%
      dplyr::mutate(ensembl_gene_id =
                      stringr::str_replace(ensembl_gene_id,
                                           "\\.[0-9]{1,}$","")) %>%
      dplyr::mutate(ensembl_protein_id =
                      stringr::str_replace(ensembl_protein_id,
                                           "\\.[0-9]{1,}$","")) %>%
      dplyr::distinct()

    ## Temporary fix: Need additional processing to get multiple tag-value
    ## pairs in a proper fashion
    gencode_gtf_fix <- as.data.frame(
      data.table::fread(destfile_gtf, skip = 5, verbose = F) %>%
        dplyr::filter(V3 == "transcript")
    )
    tag_cols <- stringr::str_match_all(
      gencode_gtf_fix$V9, "tag \\\"(\\S+)\\\";")  %>%
      purrr::map_chr(~stringr::str_c(.x[, ncol(.x)], collapse = "&"))

    gencode_gtf$gencode_tag <- tag_cols
    gencode <- gencode_gtf %>%
      dplyr::mutate(gencode_tag = dplyr::if_else(
        nchar(gencode_tag) == 0,
        as.character(NA),
        as.character(gencode_tag)
      )) %>%
      dplyr::filter(!is.na(ensembl_gene_id) &
                      !is.na(ensembl_transcript_id)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(gencode_release = gencode_version)

    #
    rlogging::message(paste0("A total of ",nrow(gencode)," transcripts parsed"))

    transcript_df <- gencode %>%
      dplyr::select(ensembl_transcript_id,
                    ensembl_gene_id,
                    symbol) %>%
      dplyr::distinct()

    ## get gene and transcript cross-references (biomart + gene_info)
    gencode_xrefs <- resolve_ensembl_transcript_xrefs(
      transcript_df = transcript_df,
      build = build,
      gene_info = gene_info
    )

    gencode <- gencode %>% dplyr::left_join(
      gencode_xrefs, by = c("ensembl_gene_id",
                            "ensembl_transcript_id",
                            "symbol"))

    gencode <- annotate_other_gencode_basic_members(gencode)
    attr(gencode, "groups") <- NULL

    if(update == T){
      saveRDS(gencode, file=gencode_rds)
    }

    rlogging::message(paste0("A total of ",nrow(gencode)," valid transcripts remaining"))

    return(gencode)

  }

annotate_other_gencode_basic_members <- function(gencode){

  ## get all genes that has an annotated basic transcript

  ## perform anti_join against the set defined below

  gencode_basicset <- gencode %>%
    dplyr::select(ensembl_gene_id, gencode_tag) %>%
    dplyr::filter(!is.na(gencode_tag) &
                    stringr::str_detect(gencode_tag,"basic"))


  single_transcripts_per_gene <- as.data.frame(
    gencode %>%
      dplyr::filter(
        stringr::str_detect(
          gencode_gene_biotype,"^(IG_|TR_)|scaRNA|snRNA|snoRNA|sRNA|scRNA|pseudogene")) %>%
      dplyr::select(ensembl_gene_id, symbol, gencode_tag,
                    ensembl_transcript_id) %>%
      dplyr::anti_join(gencode_basicset, by = "ensembl_gene_id") %>%
      dplyr::group_by(symbol, ensembl_gene_id) %>%
      dplyr::summarise(n = dplyr::n(),
                       gencode_tag = paste(gencode_tag, collapse=","),
                       .groups = "drop") %>%
      dplyr::filter(n == 1) %>%
      dplyr::filter(!stringr::str_detect(gencode_tag,"basic")) %>%
      dplyr::select(ensembl_gene_id, symbol, gencode_tag) %>%
      dplyr::mutate(basic_expanded = T)
  )

  gencode <- gencode %>%
    dplyr::left_join(
      dplyr::select(single_transcripts_per_gene,
                    ensembl_gene_id, basic_expanded),
      by = "ensembl_gene_id") %>%
    dplyr::mutate(gencode_tag = dplyr::case_when(
      is.na(gencode_tag) & basic_expanded == T ~ "basic",
      !is.na(gencode_tag) &
        !stringr::str_detect(gencode_tag,"basic") &
        basic_expanded == T ~ paste0("basic&", gencode_tag),
      TRUE ~ as.character(gencode_tag)
    )) %>%
    dplyr::select(-basic_expanded)

  return(gencode)

}



resolve_ensembl_transcript_xrefs <- function(
  transcript_df = NULL,
  build = "grch38",
  gene_info = NULL){

  invisible(assertable::assert_colnames(
    transcript_df, colnames = c('ensembl_transcript_id',
                                'ensembl_gene_id',
                                'symbol'),
    quiet = T))

  ensembl_mart <- NULL
  if(build == 'grch38'){
    ensembl_mart <- biomaRt::useEnsembl(biomart = 'genes',
                            dataset = 'hsapiens_gene_ensembl',
                            version = 105)
  }
  if(build == 'grch37'){
    ensembl_mart <- biomaRt::useEnsembl(biomart = 'genes',
                                        dataset = 'hsapiens_gene_ensembl',
                                        version = "GRCh37")
  }


  queryAttributes1 <- c('ensembl_transcript_id',
                        'refseq_mrna',
                        'ensembl_gene_id',
                        'hgnc_id',
                        'entrezgene_id')
  queryAttributes2 <- c('ensembl_transcript_id',
                        'uniprotswissprot',
                        'refseq_peptide')

  xref_biomart_1 <- biomaRt::getBM(attributes = queryAttributes1,
                                   mart = ensembl_mart) %>%
    dplyr::rename(entrezgene = entrezgene_id)
  xref_biomart_2 <- biomaRt::getBM(attributes = queryAttributes2,
                                   mart = ensembl_mart) %>%
    dplyr::rename(uniprot_acc = uniprotswissprot)
  xref_biomart <- xref_biomart_1 %>% dplyr::left_join(
    xref_biomart_2, by = "ensembl_transcript_id"
  ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(hgnc_id = as.integer(
      stringr::str_replace(hgnc_id, "HGNC:","")
    ))

  for(n in c('refseq_peptide',
             'uniprot_acc',
             'refseq_mrna')){
    xref_biomart[!is.na(xref_biomart[,n]) &
                   (xref_biomart[,n] == "" |
                      xref_biomart[,n] == "NA"),][,n] <- NA

  }

  for(xref in c('refseq_peptide',
                'uniprot_acc',
                'refseq_mrna')){

    ensXref <- as.data.frame(
      xref_biomart %>%
        dplyr::select(ensembl_gene_id,
                      ensembl_transcript_id,
                      !!rlang::sym(xref)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(ensembl_transcript_id,
                        ensembl_gene_id) %>%
        dplyr::summarise(
          !!rlang::sym(xref) := paste(unique(!!rlang::sym(xref)),
                                      collapse = "&"),
          .groups = "drop") %>%
        dplyr::mutate(!!rlang::sym(xref) := dplyr::if_else(
          !!rlang::sym(xref) == "NA", as.character(NA),
          as.character(!!rlang::sym(xref))
        ))
    )
    transcript_df <- transcript_df %>%
      dplyr::left_join(ensXref, by = c("ensembl_gene_id",
                                       "ensembl_transcript_id"))
  }

  for(xref in c('entrezgene',
                'hgnc_id')){

    ensXref <- as.data.frame(
      xref_biomart %>%
        dplyr::select(ensembl_gene_id,

                      !!rlang::sym(xref)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(
          ensembl_gene_id) %>%
        dplyr::summarise(
          !!rlang::sym(xref) := paste(unique(!!rlang::sym(xref)),
                                      collapse = "&"),
          .groups = "drop") %>%
        dplyr::mutate(!!rlang::sym(xref) := dplyr::if_else(
          !!rlang::sym(xref) == "NA", as.character(NA),
          as.character(!!rlang::sym(xref))
        ))
    )
    transcript_df <- transcript_df %>%
      dplyr::left_join(ensXref, by = c("ensembl_gene_id"))
  }


  gene_xrefs_maps <- list()

  ## map gene cross-references (name) by hgnc identifier
  gene_xrefs_maps[['by_hgnc']] <- as.data.frame(
    transcript_df %>%
      dplyr::filter(!is.na(hgnc_id)) %>%
      dplyr::select(-entrezgene) %>%
      dplyr::left_join(
        dplyr::select(gene_info,
                      hgnc_id,
                      entrezgene,
                      symbol_entrez,
                      name),
        by = "hgnc_id")
  )

  ## map gene cross-references by Entrez gene identifier
  ## given that hgnc id is not provided
  gene_xrefs_maps[['by_entrezgene']] <- as.data.frame(
    transcript_df %>%
      dplyr::filter(is.na(hgnc_id) &
                      !is.na(entrezgene) &
                      !stringr::str_detect(entrezgene, "&")) %>%
      dplyr::select(-hgnc_id) %>%
      dplyr::mutate(entrezgene = as.integer(entrezgene)) %>%
      dplyr::left_join(
        dplyr::select(gene_info,
                      hgnc_id,
                      entrezgene,
                      symbol_entrez,
                      name),
        by = "entrezgene")
  )

  ## map gene cross-references by symbol identifier
  ## given that hgnc id and entrezgene is not provided
  gene_xrefs_maps[['by_symbol']] <- as.data.frame(
    transcript_df %>%
      dplyr::filter(is.na(hgnc_id) &
                      is.na(entrezgene)) %>%
      dplyr::select(-c(entrezgene, hgnc_id)) %>%
      dplyr::left_join(
        dplyr::select(gene_info,
                      hgnc_id,
                      entrezgene,
                      symbol,
                      symbol_entrez,
                      name),
        by = "symbol")
  )

  gene_xrefs_maps[['remain']] <- as.data.frame(
    transcript_df %>%
      dplyr::filter(is.na(hgnc_id) &
                      !is.na(entrezgene) &
                      stringr::str_detect(entrezgene, "&")) %>%
      dplyr::mutate(name = NA, symbol_entrez = NA)
  )


  gencode_transcripts_xref <-
    do.call(rbind, gene_xrefs_maps)

  ## Map remaining gene cross-refs against unambiguous gene alias
  alias2gene <- as.data.frame(
    gene_info %>%
      dplyr::select(name, synonyms, entrezgene,
                    hgnc_id, symbol_entrez) %>%
      #dplyr::rename(symbol_entrez = symbol) %>%
      dplyr::mutate(entrezgene = as.character(entrezgene)) %>%
      dplyr::filter(!is.na(synonyms)) %>%
      tidyr::separate_rows(synonyms, sep = "\\|")
  )

  UnAmbigAliases <- as.data.frame(
    alias2gene %>%
      dplyr::group_by(synonyms) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::filter(n == 1) %>%
      dplyr::inner_join(alias2gene, by = "synonyms")) %>%
    dplyr::select(-n) %>%
    dplyr::rename(symbol = synonyms)

  gene_xrefs_maps[['by_alias']] <- gencode_transcripts_xref %>%
    dplyr::filter(is.na(entrezgene) & is.na(name) & !is.na(symbol))  %>%
    dplyr::select(-c(entrezgene, name, hgnc_id, symbol_entrez)) %>%
    dplyr::left_join(UnAmbigAliases, by = "symbol")

  gencode_transcripts_xref_final <- as.data.frame(
    dplyr::anti_join(gencode_transcripts_xref, gene_xrefs_maps[['by_alias']],
                     by = c("ensembl_gene_id","ensembl_transcript_id","symbol")) %>%
      dplyr::bind_rows(gene_xrefs_maps[['by_alias']])
  )

  rownames(gencode_transcripts_xref_final) <- NULL
  attr(gencode_transcripts_xref_final, "groups") <- NULL

  return(gencode_transcripts_xref_final)

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


get_cancer_hallmarks <- function(basedir = NULL,
                                 gene_info = NULL,
                                 opentargets_version = "2021.11",
                                 update = T){

  rds_fname <-
    file.path(basedir, "data-raw", "opentargets",
              paste0("opentargets_hallmarkdata_",
              opentargets_version, ".rds"))

  if(update == F & file.exists(rds_fname)){
    hallmark_data <- readRDS(file = rds_fname)
    return(hallmark_data)
  }

  opentargets_target_rds <- file.path(
    basedir, "data-raw", "opentargets",
    paste0("opentargets_target_", opentargets_version, ".rds")
  )

  hallmark_data_long <- as.data.frame(
    readRDS(file = opentargets_target_rds) %>%
      dplyr::select(target_ensembl_gene_id, hgnc_id, cancer_hallmark) %>%
      dplyr::rename(ensembl_gene_id = target_ensembl_gene_id) %>%
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
      dplyr::select(hgnc_id, ensembl_gene_id, hallmark,
                    promotes, suppresses, pmid_summary,
                    link) %>%
      dplyr::mutate(literature_support = paste0(pmid_summary," (",link,")")) %>%
      dplyr::group_by(hgnc_id, ensembl_gene_id, hallmark,
                      promotes, suppresses) %>%
      dplyr::summarise(literature_support = paste(
        literature_support, collapse=". "),
        .groups = "drop") %>%
      dplyr::left_join(dplyr::select(gene_info, hgnc_id,
                                     symbol, entrezgene),
                       by = "hgnc_id") %>%
      dplyr::mutate(symbol =
                      paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                             entrezgene,"' target='_blank'>",symbol,"</a>")) %>%
      dplyr::select(-c(entrezgene, hgnc_id)) %>%
      dplyr::select(symbol, ensembl_gene_id, hallmark, promotes, suppresses,
                    literature_support)
  )

  hallmark_data <- list('long' = hallmark_data_long,
                        'short' = hallmark_data_short)
  saveRDS(hallmark_data, file = rds_fname)
  return(hallmark_data)

}


get_gencode_data <- function(
  basedir = NULL,
  gencode_version = "39",
  build = "grch38",
  gene_info = NULL,
  update = T){

  rds_fname <- file.path(
    basedir, "data-raw",
    "gencode", paste0("gencode_OE_", gencode_version, ".rds")
  )

  if(update == F & file.exists(rds_fname)){
    gencode <- readRDS(file = rds_fname)
    return(gencode)
  }

  gencode <- as.data.frame(
    get_gencode_ensembl_transcripts(
      basedir = basedir,
      build = build,
      gene_info = gene_info,
      gencode_version = gencode_version,
      update = update) %>%
      dplyr::select(ensembl_gene_id,
                    symbol,
                    symbol_entrez,
                    refseq_mrna,
                    refseq_peptide,
                    uniprot_acc,
                    name,
                    ensembl_transcript_id,
                    ensembl_protein_id,
                    entrezgene,
                    gencode_gene_biotype) %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(entrezgene) & !is.na(symbol_entrez))
  )

  saveRDS(gencode, file = rds_fname)


  return(gencode)
}

get_tf_target_interactions <- function(
  basedir = NULL,
  update = F){


  rds_fname = file.path(
    basedir, "data-raw",
    "omnipathdb", "omnipath_tf_interactions.rds")

  if(update == F & file.exists(rds_fname)){
    tf_target_interactions <- readRDS(file = rds_fname)
    return(tf_target_interactions)
  }


  tf_target_interactions_dorothea_all <-
    OmnipathR::import_transcriptional_interactions(
      organism = "9606",
      dorothea_levels = c("A","B","C","D"),
      references_by_resource = F) %>%
    dplyr::filter(!is.na(dorothea_level)) %>%
    dplyr::group_by(source_genesymbol, target_genesymbol) %>%
    dplyr::summarise(
      references = paste(unique(references), collapse=";"),
      interaction_sources = paste(unique(sources), collapse="; "),
      dorothea_level = paste(unique(dorothea_level), collapse=";"),
      .groups = "drop") %>%
    dplyr::filter(!stringr::str_detect(source_genesymbol,"_")) %>%
    dplyr::mutate(interaction_sources = stringr::str_squish(
      stringr::str_replace_all(
        interaction_sources, ";","; "))) %>%

    ## Not available for non-academic use (entries provided solely with TRED)
    dplyr::filter(!stringr::str_detect(interaction_sources,"^(RegNetwork_DoRothEA; TRED_DoRothEA)$"))

  tf_target_pmids <- tf_target_interactions_dorothea_all %>%
    dplyr::filter(!is.na(references) & nchar(references) > 0) %>%
    dplyr::select(source_genesymbol, target_genesymbol, references) %>%
    tidyr::separate_rows(references, sep=";") %>%
    dplyr::filter(nchar(references) > 0) %>%
    dplyr::distinct()

  tf_target_citations <- get_citations_pubmed(
    pmid = unique(tf_target_pmids$references))

  tf_target_pmids <- as.data.frame(
    tf_target_pmids %>%
      dplyr::mutate(references = as.integer(references)) %>%
      dplyr::left_join(tf_target_citations, by = c("references" = "pmid")) %>%
      dplyr::group_by(source_genesymbol, target_genesymbol) %>%
      dplyr::summarise(
        tf_target_literature_support = paste(
          link, collapse= ","),
        tf_target_literature = paste(
          citation, collapse="|"
        )
      )
  )

  tf_target_interactions_dorothea_all <- tf_target_interactions_dorothea_all %>%
    dplyr::left_join(tf_target_pmids,
                     by = c("source_genesymbol",
                            "target_genesymbol")) %>%
    dplyr::rename(regulator = source_genesymbol,
                  target = target_genesymbol) %>%
    dplyr::select(-references) %>%
    dplyr::mutate(confidence_level = dplyr::case_when(
      stringr::str_detect(dorothea_level,"A|A;B|A;B;D|A;D") ~ "A",
      stringr::str_detect(dorothea_level,"B|B;D|B;C") ~ "B",
      stringr::str_detect(dorothea_level,"C;D|C") ~ "C",
      stringr::str_detect(dorothea_level,"D") ~ "D",
      TRUE ~ as.character(dorothea_level),
    )) %>%
    dplyr::select(-dorothea_level) %>%
    dplyr::mutate(mode_of_regulation = NA)


  tf_target_interactions_dorothea_pancancer <-
    dorothea::dorothea_hs_pancancer %>%
    dplyr::mutate(interaction_sources = "DoRothEA_hs_pancancer") %>%
    dplyr::rename(regulator = tf,
                  confidence_level = confidence,
                  mode_of_regulation = mor) %>%
    dplyr::mutate(mode_of_regulation = dplyr::if_else(
      mode_of_regulation == 1,"Stimulation",
      "Repression")) %>%
    dplyr::mutate(tf_target_literature_support = NA,
                  tf_target_literature = NA) %>%
    dplyr::select(regulator, target,
                  interaction_sources,
                  confidence_level,
                  mode_of_regulation,
                  tf_target_literature_support,
                  tf_target_literature)

  tf_target_interactions <- list()
  tf_target_interactions[['global']] <-
    tf_target_interactions_dorothea_all
  tf_target_interactions[['pancancer']] <-
    tf_target_interactions_dorothea_pancancer

  saveRDS(tf_target_interactions, file = rds_fname)

  return(tf_target_interactions)
}

get_function_summary_ncbi <- function(
  basedir = NULL,
  gene_df = NULL,
  update = F){

  rds_fname <- file.path(
    basedir, "data-raw",
    "ncbi_gene_summary.rds")

  if(update == F & file.exists(rds_fname)){
    ncbi_gene_summary <- readRDS(rds_fname)
    return(ncbi_gene_summary)
  }
  assertable::assert_colnames(
    gene_df, c("entrezgene","gene_biotype"), only_colnames = F, quiet = T)

  pcg <- gene_df %>%
    dplyr::filter(gene_biotype == "protein-coding") %>%
    dplyr::select(entrezgene) %>%
    dplyr::distinct()

  i <- 1
  ncbi_gene_summary <- data.frame()
  while(i < NROW(pcg) - 300){

    queryset_stop <- i + 299
    queryset <- pcg[i:queryset_stop,"entrezgene"]

    summary_results <- suppressMessages(as.data.frame(
      mygene::queryMany(queryset,
                        scopes = "entrezgene",
                        fields= "summary",
                        species = "human")
    )) %>%
      dplyr::select(query, summary) %>%
      dplyr::rename(entrezgene = query,
                    gene_summary_ncbi = summary) %>%
      dplyr::mutate(entrezgene = as.integer(entrezgene))

    ncbi_gene_summary <-
      dplyr::bind_rows(ncbi_gene_summary,
                       summary_results)

    i <- i + 300

  }

  queryset <- pcg[i:nrow(pcg),"entrezgene"]

  summary_results <- suppressMessages(
    as.data.frame(
      mygene::queryMany(queryset,
                        scopes = "entrezgene",
                        fields="summary",
                        species="human")
    )) %>%
    dplyr::select(query, summary) %>%
    dplyr::rename(entrezgene = query,
                  gene_summary_ncbi = summary) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene))


  ncbi_gene_summary <-
    dplyr::bind_rows(ncbi_gene_summary,
                     summary_results)

  ncbi_gene_summary <- ncbi_gene_summary %>%
    dplyr::mutate(
      gene_summary_ncbi = stringr::str_replace(
        gene_summary_ncbi,
        " \\[provided by RefSeq, [A-Za-z]{3} [0-9]{4}\\]\\.",
        "")
    ) %>%
    dplyr::mutate(
      gene_summary_ncbi = stringr::str_replace(
        gene_summary_ncbi,
        "\\[supplied by OMIM, [A-Za-z]{3} [0-9]{4}\\]\\.",
        "")
    ) %>%
    dplyr::mutate(
      gene_summary_ncbi = stringr::str_replace(
        gene_summary_ncbi,
        " \\[PubMed [0-9]{1,}\\]",
        "")
    ) %>%
    dplyr::mutate(
      gene_summary_ncbi =
        dplyr::if_else(
          !is.na(gene_summary_ncbi),
          paste0("<b>NCBI/RefSeq/OMIM:</b> ",
                 gene_summary_ncbi),
          as.character(gene_summary_ncbi)
        )
    )

  saveRDS(ncbi_gene_summary, file = rds_fname)
  return(ncbi_gene_summary)

}

get_cancer_drugs <- function(){

  cancer_drugs <- list()
  target_drug_edges <- list()
  target_drug_nodes <- list()
  drugs_per_gene <- list()
  for(p in c('early_phase','late_phase')){
    cancer_drugs[[p]] <- data.frame()
    target_drug_edges[[p]] <- data.frame()
    target_drug_nodes[[p]] <- data.frame()
    drugs_per_gene[[p]] <- data.frame()
  }


  cancer_drugs[['late_phase']] <-
    oncoPharmaDB::get_onco_drugs(
      drug_is_targeted = T,
      source_opentargets_only = T,
      drug_minimum_phase_any_indication = 3) %>%
    dplyr::filter(!is.na(molecule_chembl_id)) %>%
    dplyr::select(target_entrezgene,
                  drug_name,
                  primary_site,
                  molecule_chembl_id,
                  drug_max_ct_phase) %>%
    dplyr::distinct() %>%
    dplyr::rename(entrezgene = target_entrezgene) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene)) %>%
    dplyr::distinct()

  cancer_drugs[['early_phase']] <-
    oncoPharmaDB::get_onco_drugs(
      drug_is_targeted = T,
      source_opentargets_only = T,
      drug_minimum_phase_any_indication = 0) %>%
    dplyr::filter(!is.na(molecule_chembl_id)) %>%
    dplyr::anti_join(
      cancer_drugs[['late_phase']], by = "molecule_chembl_id") %>%
    dplyr::select(target_entrezgene,
                  drug_name,
                  primary_site,
                  molecule_chembl_id,
                  drug_max_ct_phase) %>%
    dplyr::rename(entrezgene = target_entrezgene) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene)) %>%
    dplyr::distinct()

  for(p in c('early_phase','late_phase')){

    cancer_drugs[[p]] <- cancer_drugs[[p]] %>%
      dplyr::mutate(to = paste0("s", molecule_chembl_id)) %>%
      dplyr::mutate(from = paste0("s", entrezgene)) %>%
      dplyr::distinct()

    drugs_per_gene[[p]] <- as.data.frame(
      cancer_drugs[[p]] %>%
        dplyr::select(entrezgene,
                      drug_name,
                      drug_max_ct_phase,
                      molecule_chembl_id) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          drug_link =
            paste0("<a href = 'https://platform.opentargets.org/drug/",
                   molecule_chembl_id,"' target='_blank'>",drug_name,"</a>")) %>%
        dplyr::arrange(entrezgene, desc(drug_max_ct_phase)) %>%
        dplyr::group_by(entrezgene) %>%
        dplyr::summarise(
          targeted_cancer_drugs =
            paste(unique(drug_link),
                  collapse = ", "), .groups = "drop")
    )
    if(p == 'early_phase'){
      drugs_per_gene[[p]] <- drugs_per_gene[[p]] %>%
        dplyr::rename(
          targeted_cancer_drugs_ep = targeted_cancer_drugs)
    }else{
      drugs_per_gene[[p]] <- drugs_per_gene[[p]] %>%
        dplyr::rename(
          targeted_cancer_drugs_lp = targeted_cancer_drugs)
    }

    target_drug_edges[[p]] <- cancer_drugs[[p]] %>%
      dplyr::select(drug_name,
                    from,
                    to,
                    #target_symbol,
                    entrezgene,
                    molecule_chembl_id) %>%
      dplyr::distinct() %>%
      dplyr::mutate(width = 1)

    target_drug_nodes[[p]] <- as.data.frame(
      cancer_drugs[[p]] %>%
        dplyr::select(to,
                      drug_name,
                      primary_site) %>%
        dplyr::distinct() %>%
        dplyr::rename(id = to,
                      label = drug_name) %>%
        dplyr::group_by(id, label) %>%
        dplyr::summarise(
          primary_site = paste(sort(unique(primary_site)),
                               collapse = ", "),
          .groups = "drop") %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          primary_site =
            dplyr::if_else(is.na(primary_site) |
                             nchar(primary_site) == 0,
                           "Unspecific indication",
                           as.character(primary_site))) %>%
        dplyr::mutate(
          title = paste0(label, " (", primary_site,")")) %>%
        dplyr::mutate(
          shadow = T, shape = "diamond",
          color.background = "orange",
          font.color = "black")
    )
    if(p == 'early_phase'){
      target_drug_nodes[[p]] <- target_drug_nodes[[p]] %>%
        dplyr::mutate(shadow = T,
                      shape = "diamond",
                      color.background = "purple",
                      font.color = "black")
    }

  }
  cancerdrugdb <- list()
  cancerdrugdb[['drug_per_target']] <- list()
  cancerdrugdb[['drug_per_target']][['early_phase']] <-
    drugs_per_gene[['early_phase']]
  cancerdrugdb[['drug_per_target']][['late_phase']] <-
    drugs_per_gene[['late_phase']]

  cancerdrugdb[['ppi']] <- list()
  cancerdrugdb[['ppi']][['edges']] <-
    dplyr::bind_rows(target_drug_edges[['early_phase']],
                     target_drug_edges[['late_phase']]) %>%
    dplyr::select(from, to, width) %>%
    dplyr::distinct()

  cancerdrugdb[['ppi']][['nodes']] <-
    dplyr::bind_rows(target_drug_nodes[['early_phase']],
                     target_drug_nodes[['late_phase']]) %>%
    dplyr::distinct()

  return(cancerdrugdb)


}

get_pathway_annotations <- function(
  basedir = NULL,
  gene_info = NULL,
  wikipathways = T,
  kegg = T,
  netpath = T,
  msigdb = T,
  omnipathdb = NULL,
  msigdb_version = NULL,
  wikipathways_version = NULL,
  kegg_version = NULL,
  netpath_version = NULL){


  wikipathways_gmt <- file.path(
    basedir, "data-raw", "wikipathways",
    paste0("wikipathways-", wikipathways_version,
           "-gmt-Homo_sapiens.gmt")
  )

  netpath_mapping_fname <- file.path(
    basedir, "data-raw", "netpath",
    "id_mapping.tsv")

  keggdb_fname <- file.path(
    basedir, "data-raw",
    "kegg", "kegg.pathway.gene.tsv"
  )

  assertable::assert_colnames(
    gene_info,
    c("symbol","entrezgene"),
    only_colnames = F,
    quiet = T)

  pathwaydb <- list()

  if(wikipathways == T){
    if(!file.exists(wikipathways_gmt)){
      rWikiPathways::downloadPathwayArchive(
        organism = "Homo sapiens",
        date = wikipathways_version,
        destpath = "data-raw/wikipathways/",
        format = "gmt")
    }

    wikipathways <- clusterProfiler::read.gmt(
      wikipathways_gmt)
    wp2gene <- wikipathways %>%
      tidyr::separate(term, c("name","version","wpid","org"), "%")
    wikipathwaydb <- list()
    wikipathwaydb[['VERSION']] <- wikipathways_version
    wikipathwaydb[['TERM2GENE']] <- wp2gene %>%
      dplyr::select(wpid, gene) %>%
      dplyr::rename(standard_name = wpid, entrez_gene = gene)
    wikipathwaydb[['TERM2NAME']] <- wp2gene %>%
      dplyr::select(wpid, name) %>%
      dplyr::rename(standard_name = wpid) %>%
      dplyr::distinct()

    pathwaydb[['wikipathways']] <- wikipathwaydb
  }


  if(netpath == T & !is.null(omnipathdb) & !is.null(gene_info)){
      netpath_idmapping <-
      read.table(netpath_mapping_fname,
                 stringsAsFactors = F, header = F, sep = "\t",
                 col.names = c("name", "standard_name")) %>%
      dplyr::mutate(standard_name = paste0("NetPath_",standard_name))

    netpath_pathway_data <- omnipathdb %>%
      dplyr::filter(source == "NetPath") %>%
      dplyr::select(genesymbol, value) %>%
      dplyr::rename(name = value, symbol = genesymbol) %>%
      dplyr::left_join(netpath_idmapping, by = "name") %>%
      dplyr::mutate(name = paste0(name," signaling pathway")) %>%
      dplyr::arrange(name, standard_name, symbol) %>%
      dplyr::left_join(dplyr::select(
        gene_info, entrezgene, symbol_entrez),
        by = c("symbol" = "symbol_entrez")) %>%
      dplyr::rename(entrez_gene = entrezgene) %>%
      dplyr::mutate(entrez_gene = as.character(entrez_gene)) %>%
      dplyr::select(standard_name, name, entrez_gene)
    netpathdb <- list()
    netpathdb[['VERSION']] <- netpath_version
    netpathdb[['TERM2GENE']] <- netpath_pathway_data %>%
      dplyr::select(standard_name, entrez_gene)
    netpathdb[['TERM2NAME']] <- netpath_pathway_data %>%
      dplyr::select(standard_name, name) %>%
      dplyr::distinct()

    pathwaydb[['netpath']] <- netpathdb
  }


  if(kegg == T){
    kegg_pathway_genes <-
      read.table(file = keggdb_fname,
                 header=T,sep="\t",stringsAsFactors = F,quote="")
    keggdb <- list()
    keggdb[['VERSION']] <- kegg_version
    keggdb[['TERM2NAME']] <- kegg_pathway_genes %>%
      dplyr::select(name,pathway_id) %>%
      dplyr::rename(standard_name = pathway_id) %>%
      dplyr::select(standard_name, name) %>%
      dplyr::distinct()

    keggdb[['TERM2GENE']] <- kegg_pathway_genes %>%
      dplyr::select(pathway_id, gene_id) %>%
      dplyr::rename(standard_name = pathway_id,
                    entrez_gene = gene_id) %>%
      dplyr::select(standard_name,entrez_gene) %>%
      dplyr::distinct()

    pathwaydb[['kegg']] <- keggdb
  }

  if(msigdb == T){

    msigdb <-
      get_msigdb_signatures(basedir = basedir,
                            db_version = msigdb_version)

    pathwaydb[['msigdb']] <- msigdb$db

  }

  return(pathwaydb)

}

get_protein_complexes <- function(
    basedir = NULL,
    update = F){

    rds_fname <- file.path(
      basedir, "data-raw",
      "omnipathdb", "omnipath_complexdb.rds")

    humap2_complexes_tsv <-file.path(
      basedir, "data-raw",
      "humap2", "humap2_complexes.tsv")

    corum_complexes_tsv <- file.path(
      basedir, "data-raw",
      "corum", "coreComplexes.txt")

    humap_url <-
      "http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt"


    if(update == F & file.exists(rds_fname)){
      proteincomplexdb <- readRDS(rds_fname)
      return(proteincomplexdb)
    }

    complexes_curated <- OmnipathR::import_omnipath_complexes() %>%
      dplyr::filter(!is.na(references) & !is.na(name)) %>%
      dplyr::rename(pmid = references,
                    uniprot_acc = components,
                    complex_name = name) %>%
      dplyr::select(complex_name, uniprot_acc,
                    pmid, sources) %>%
      dplyr::filter(stringr::str_detect(uniprot_acc,"_")) %>%
      dplyr::distinct() %>%
      dplyr::arrange(complex_name) %>%
      dplyr::group_by(complex_name) %>%
      dplyr::mutate(id2 = dplyr::row_number()) %>%
      dplyr::group_by(complex_name, uniprot_acc) %>%
      dplyr::summarise(pmid = paste(sort(unique(pmid)), collapse=","),
                       sources = paste(sort(unique(sources)), collapse=";"),
                       id2 = paste(id2, collapse=";"),
                       .groups = "drop") %>%
      dplyr::mutate(num_sources = stringr::str_count(sources,";") + 1) %>%
      dplyr::arrange(complex_name, desc(nchar(uniprot_acc)), desc(num_sources)) %>%
      dplyr::mutate(id2 = dplyr::row_number()) %>%
      #dplyr::filter(id2 == 1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(complex_id = dplyr::row_number()) %>%
      dplyr::mutate(complex_name = stringr::str_replace_all(
        complex_name,"\\\\u2013","-"
      )) %>%
      dplyr::mutate(sources = stringr::str_replace_all(
        sources, ";", "; ")) %>%
      dplyr::select(complex_id, complex_name, uniprot_acc,
                    sources, pmid)

    pmid2complex <- complexes_curated %>%
      dplyr::select(complex_id,
                    pmid) %>%
      tidyr::separate_rows(pmid, sep=";") %>%
      dplyr::filter(stringr::str_detect(pmid,"^[0-9]{4,}$"))

    protein_complex_citations <- get_citations_pubmed(
      pmid = unique(pmid2complex$pmid))

    pmid2complex <- pmid2complex %>%
      dplyr::mutate(pmid = as.integer(pmid)) %>%
      dplyr::left_join(protein_complex_citations, by = "pmid") %>%
      dplyr::group_by(complex_id) %>%
      dplyr::summarise(
        complex_literature_support = paste(
          link, collapse= ","),
        complex_literature = paste(
          citation, collapse="|"
        )
      )

    complex_omnipath <- complexes_curated %>%
      dplyr::left_join(pmid2complex,
                       by = "complex_id")

    ## CORUM annotations (complex comments, disease comments, purification methods)
    corum_complex_annotations <- as.data.frame(
      readr::read_tsv(
        corum_complexes_tsv,
        show_col_types = F) %>%
        janitor::clean_names() %>%
        dplyr::filter(organism == "Human") %>%
        dplyr::group_by(complex_name) %>%
        dplyr::summarise(
          purification_method =
            paste(unique(protein_complex_purification_method),
                  collapse = "; "),
          complex_comment = paste(unique(complex_comment),
                                  collapse="; "),
          disease_comment = paste(unique(disease_comment),
                                  collapse="; ")
        ) %>%
        dplyr::distinct()
    )

    complex_omnipath <- complex_omnipath %>%
      dplyr::left_join(corum_complex_annotations, by = "complex_name") %>%
      dplyr::select(-pmid) %>%
      dplyr::mutate(confidence = NA) %>%
      dplyr::mutate(complex_id = as.character(complex_id))

    if(!file.exists(humap2_complexes_tsv)){
      download.file(url = humap_url,
                    destfile = humap2_complexes_tsv)
    }

    complex_humap <- readr::read_csv(
      file = humap2_complexes_tsv,
      show_col_types = F
    ) %>%
      janitor::clean_names() %>%
      dplyr::rename(complex_id = hu_map2_id,
                    uniprot_acc = uniprot_ac_cs) %>%
      dplyr::select(-genenames) %>%
      dplyr::mutate(uniprot_acc = stringr::str_replace_all(
        uniprot_acc, " ","_"
      )) %>%
      dplyr::mutate(complex_name = complex_id) %>%
      dplyr::mutate(sources = "hu.MAP 2.0")

    complex_all <- complex_omnipath %>%
      dplyr::bind_rows(complex_humap)

    proteincomplexdb <- list()
    proteincomplexdb[["db"]] <- complex_all %>%
      dplyr::select(-uniprot_acc) %>%
      dplyr::distinct()

    proteincomplexdb[["up_xref"]] <- complex_all %>%
      dplyr::select(complex_id, uniprot_acc) %>%
      tidyr::separate_rows(uniprot_acc, sep="_") %>%
      dplyr::distinct()


    saveRDS(proteincomplexdb,
            file = rds_fname)
    return(proteincomplexdb)
}

get_crispr_scores <- function(
    basedir = NULL,
    gene_info = NULL){

    cell_lines <- read.table(
      file = file.path(
        basedir, "data-raw",
        "project_score",
        "model_list.csv"),
      quote="",
      comment.char="",sep=",", header = F, stringsAsFactors = F,
      fill = T,na.strings = c("")) %>%
      dplyr::select(V1, V2, V3, V4, V18, V22, V24,
                    V25, V32, V40, V39, V41) %>%
      dplyr::filter(V4 == "Cell Line" & V39 == "Homo Sapiens") %>%
      dplyr::filter(!stringr::str_detect(V1,"model_id")) %>%
      dplyr::mutate(V2 = stringr::str_replace_all(V2,"-",".")) %>%
      dplyr::filter(V18 == "True") %>%
      magrittr::set_colnames(
        c('model_id','model_name','synonyms','model_type',
          'crispr_ko_data','tissue','cancer_type',
          'cancer_type_detail','sample_site',
          'gender','species','ethnicity')) %>%
      dplyr::select(
        -c(gender,species,ethnicity,cancer_type_detail,cancer_type,
           crispr_ko_data,synonyms)) %>%
      dplyr::filter(tissue != "Eye" & tissue != "Biliary Tract" &
                      tissue != "Prostate" & tissue != "Soft Tissue") %>%
      dplyr::mutate(tissue = dplyr::if_else(tissue == "Central Nervous System",
                                            "CNS/Brain",
                                            as.character(tissue))) %>%
      dplyr::mutate(tissue = dplyr::if_else(tissue == "Large Intestine",
                                            "Colon/Rectum",
                                            as.character(tissue)))

    gene_identifiers <-
      read.csv(file = file.path(
        basedir, "data-raw",
        "project_score",
        "gene_identifiers.csv"),
        stringsAsFactors = F) %>%
      dplyr::rename(gene_id_project_score = gene_id,
                    entrezgene = entrez_id) %>%
      dplyr::filter(!is.na(entrezgene)) %>%
      dplyr::left_join(
        dplyr::select(gene_info, symbol,
                      entrezgene), by = "entrezgene") %>%
      dplyr::select(gene_id_project_score, entrezgene, symbol)

    bayes_factors <-
      read.table(file = file.path(
        basedir, "data-raw",
        "project_score", "03_scaledBayesianFactors.tsv"),
        header=T, quote="",sep="\t", stringsAsFactors = F)

    rownames(bayes_factors) <- bayes_factors$Gene
    bayes_factors$Gene <- NULL
    bayes_factors <- as.matrix(bayes_factors)
    bayes_factors <-
      as.data.frame(setNames(reshape2::melt(bayes_factors, na.rm = T),
                             c('symbol', 'model_name', 'scaled_BF'))) %>%
      dplyr::mutate(model_name = as.character(model_name),
                    symbol = as.character(symbol)) %>%
      dplyr::filter(scaled_BF > 0) %>%
      dplyr::left_join(cell_lines, by = c("model_name")) %>%
      dplyr::filter(!is.na(model_id)) %>%
      dplyr::left_join(gene_identifiers,by = c("symbol")) %>%
      dplyr::filter(!is.na(gene_id_project_score))


    ## Fitness scores
    cell_fitness_scores <-
      read.table(
        file = file.path(
          basedir, "data-raw",
          "project_score",
          "binaryDepScores.tsv"),header=T,
        quote="",sep="\t", stringsAsFactors = F)

    rownames(cell_fitness_scores) <- cell_fitness_scores$Gene
    cell_fitness_scores$Gene <- NULL
    cell_fitness_scores <- as.matrix(cell_fitness_scores)
    projectscoredb <- list()
    projectscoredb[['fitness_scores']] <-
      as.data.frame(setNames(reshape2::melt(cell_fitness_scores, na.rm = T),
                             c('symbol', 'model_name', 'loss_of_fitness'))) %>%
      dplyr::mutate(model_name = as.character(model_name), symbol = as.character(symbol)) %>%
      dplyr::filter(loss_of_fitness == 1) %>%
      dplyr::left_join(cell_lines, by = c("model_name")) %>%
      dplyr::filter(!is.na(model_id)) %>%
      dplyr::left_join(gene_identifiers,by=c("symbol")) %>%
      dplyr::filter(!is.na(gene_id_project_score)) %>%
      dplyr::left_join(bayes_factors, by =
                         c("symbol", "model_name","model_id",
                           "model_type","tissue","sample_site",
                           "gene_id_project_score","entrezgene"))

    ## Target priority scores
    projectscoredb[['target_priority_scores']] <-
      read.csv(file = file.path(
        basedir, "data-raw",
        "project_score", "depmap-priority-scores.csv"),
        header = T, quote="", stringsAsFactors = F) %>%
      magrittr::set_colnames(c("tractability_bucket","gene_id","priority_score",
                               "symbol","analysis_id","tumor_type")) %>%
      dplyr::mutate(priority_score = as.numeric(
        stringr::str_replace_all(priority_score,"\\\"",""))) %>%
      dplyr::mutate(symbol = as.character(
        stringr::str_replace_all(symbol,"\\\"",""))) %>%
      dplyr::mutate(gene_id = as.character(
        stringr::str_replace_all(gene_id,"\\\"",""))) %>%
      dplyr::mutate(tumor_type = as.character(
        stringr::str_replace_all(tumor_type,"\\\"",""))) %>%
      dplyr::select(symbol, gene_id, tumor_type,
                    priority_score) %>%
      dplyr::mutate(priority_score =
                      round(priority_score,
                            digits = 3)) %>%
      dplyr::mutate(tumor_type = dplyr::case_when(
        tumor_type == "Pan-Cancer" ~ "Pancancer",
        tumor_type == "Haematopoietic and Lyphoid" ~
          "Haematopoietic and Lymphoid",
        TRUE ~ as.character(tumor_type)
      )) %>%
      dplyr::distinct()

    ## order symbols by pancancer priority rank
    priority_order <-
      projectscoredb$target_priority_scores %>%
      dplyr::filter(tumor_type == "Pancancer") %>%
      dplyr::arrange(desc(priority_score))

    projectscoredb[['target_priority_scores']] <-
      projectscoredb[['target_priority_scores']] %>%
      dplyr::mutate(symbol = factor(symbol,
                                    levels = priority_order$symbol))

    return(projectscoredb)
    #saveRDS(projectscoredb, file="db/projectscoredb.rds")

}


get_survival_associations <- function(
  gene_info = NULL,
  basedir = NULL){

  for(feature_type in
      c('CNAs','Mutations','Gene expression','Methylation',
        'miRNA expression','Protein expression')){
    destfile_fname <-
      file.path(basedir, "data-raw", "km_survival_cshl",
                paste0(
                  tolower(
                    stringr::str_replace(
                      stringr::str_replace(feature_type," ","_"),
                      "s$","")),
                  ".xlsx"))
    if(!file.exists(destfile_fname)){
      download.file(url = paste0(
        "https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/full-downloads/",
        stringr::str_replace(feature_type," ","%20"),".xlsx"),
        destfile = destfile_fname,
        quiet = T)
    }

  }

  project_survival <- list()
  for(feature_type in c('cna','mutation',
                        'gene_expression',
                        'methylation')){

    fname <- file.path(basedir, "data-raw",
                       "km_survival_cshl",
                       paste0(feature_type,".xlsx"))
    data <- openxlsx::read.xlsx(fname)

    if(feature_type == 'protein_expression'){
      rownames(data) <- data$Protein
      data$Protein <- NULL
    }else{
      rownames(data) <- data$Gene
      data$Gene <- NULL
    }

    feattype_brief <- feature_type
    if(feattype_brief == "mutation"){
      feattype_brief <- "mut"
    }
    if(feattype_brief == "gene_expression"){
      feattype_brief <- "exp"
    }
    if(feattype_brief == "methylation"){
      feattype_brief <- "meth"
    }

    data <- as.matrix(data)
    project_survival_df <-
      as.data.frame(setNames(reshape2::melt(data, na.rm = T),
                             c('symbol', 'tcga_cohort', 'z_score'))) %>%
      #dplyr::mutate(feature_type = feattype_brief) %>%
      dplyr::mutate(z_score = as.numeric(stringr::str_trim(z_score))) %>%
      dplyr::mutate(symbol = stringr::str_replace(
        symbol,"^'","")) %>%
      dplyr::filter(!stringr::str_detect(tcga_cohort,"^Stouffer")) %>%
      dplyr::left_join(dplyr::select(gene_info,
                                     symbol, gene_biotype),
                       by = "symbol") %>%
      ## limit associations to protein-coding genes
      dplyr::filter(gene_biotype == "protein-coding") %>%
      dplyr::select(-gene_biotype)

    pancancer_trend <- as.data.frame(
      project_survival_df %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(tot_z_score = sum(z_score)) %>%
        dplyr::arrange(desc(tot_z_score)))

    project_survival[[feattype_brief]] <- project_survival_df %>%
      dplyr::mutate(
        symbol = factor(symbol, levels = pancancer_trend$symbol))
  }
  projectsurvivaldb <- project_survival

  return(projectsurvivaldb)

}

get_unique_transcript_xrefs <- function(
  basedir = NULL,
  gencode = NULL,
  gene_info = NULL,
  update = F){

  rds_fname <- file.path(
    basedir, "data-raw",
    "transcript_xref", "transcript_xref_db.rds"
  )

  ## gene synonyms
  alias2entrez <-
    get_gene_aliases(gene_info = gene_info) %>%
    dplyr::mutate(entrezgene = as.character(entrezgene))


  if(update == F & file.exists(rds_fname)){
    transcript_xref_db <- readRDS(file = rds_fname)
    return(transcript_xref_db)
  }

  gencode_merged <-
    gencode$grch37 %>%
    dplyr::bind_rows(gencode$grch38) %>%
    dplyr::filter(!stringr::str_detect(entrezgene,"&"))

  transcript_xref_db <- data.frame()

  for(xref in c('ensembl_transcript_id',
                'uniprot_acc',
                'ensembl_gene_id',
                'ensembl_protein_id',
                'refseq_mrna',
                'symbol',
                'refseq_peptide')){

    tmp <- gencode_merged %>%
      dplyr::select(entrezgene, !!rlang::sym(xref)) %>%
      tidyr::separate_rows(!!rlang::sym(xref), sep ="&") %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(!!rlang::sym(xref))) %>%
      dplyr::group_by(!!rlang::sym(xref)) %>%
      dplyr::summarise(n = dplyr::n(),
                       entrezgene = paste(unique(entrezgene),collapse=",")) %>%
      dplyr::filter(n == 1) %>%
      dplyr::mutate(property = xref,
                    value = !!rlang::sym(xref)) %>%
      dplyr::select(entrezgene, property, value)

    transcript_xref_db <- dplyr::bind_rows(
      transcript_xref_db, tmp
    )

  }

  transcript_xref_db <- dplyr::bind_rows(
    transcript_xref_db, alias2entrez
  ) %>%
    dplyr::mutate(entrezgene = as.integer(entrezgene))

  saveRDS(transcript_xref_db, file = rds_fname)
  return(transcript_xref_db)

}

get_omnipath_gene_annotations <- function(
  basedir = NULL,
  gene_info = NULL,
  update = F){

  rds_fname <- file.path(
    basedir, "data-raw",
    "omnipathdb", "omnipathdb.rds")

  if(update == F & file.exists(rds_fname)){
    omnipathdb <- readRDS(file = rds_fname)
    return(omnipathdb)
  }

  protein_coding_genes <- gene_info %>%
    dplyr::filter(gene_biotype == "protein-coding") %>%
    dplyr::select(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol,"-")) %>%
    dplyr::filter(symbol != "XYLB") %>%
    dplyr::distinct()

  i <- 1
  omnipathdb <- data.frame()
  while(i <= nrow(protein_coding_genes) - 200){
    genes <-
      protein_coding_genes$symbol[i:min(nrow(protein_coding_genes),i + 199)]
    annotations <-
      OmnipathR::import_omnipath_annotations(proteins = genes) %>%
      dplyr::filter(
        !stringr::str_detect(
          source, "^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$")) %>%
      dplyr::rename(uniprot_acc = uniprot)
    omnipathdb <-
      dplyr::bind_rows(
        omnipathdb, annotations)
    cat(i, min(nrow(protein_coding_genes), i + 199), '\n')
    i <- i + 200
  }
  genes <-
    protein_coding_genes$symbol[i:nrow(protein_coding_genes)]
  annotations <-
    OmnipathR::import_omnipath_annotations(
      select_genes = genes) %>%
    dplyr::filter(!stringr::str_detect(
      source, "^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$")) %>%
    dplyr::rename(uniprot_acc = uniprot)
  if(nrow(annotations) > 0){
    omnipathdb <-
      dplyr::bind_rows(omnipathdb, annotations)
  }

  saveRDS(omnipathdb, file = rds_fname)
  return(omnipathdb)
}

get_gene_go_terms <- function(
  basedir = NULL){


  goa_human_fname <-
    file.path(basedir, "data-raw",
              "gene_ontology", "goa_human.gaf.gz")

  go_annotations_qc <-
    read.table(gzfile(goa_human_fname),
               sep = "\t", comment.char = "!", header = F,
               stringsAsFactors = F, quote = "") %>%
    dplyr::select(V3, V5, V7, V9, V14) %>%
    dplyr::rename(symbol = V3, go_id = V5,
                  go_evidence_code = V7,
                  go_ontology = V9, go_annotation_date = V14) %>%
    dplyr::distinct() %>%
    dplyr::filter(go_evidence_code != "IEA" &
                    go_evidence_code != "ND") %>%
    dplyr::filter(go_ontology != "C")

  go_function_terms_prgene <- go_annotations_qc %>%
    dplyr::filter(go_ontology == "F") %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(
      num_ids_function = length(unique(go_id)),
      .groups = "drop"
    )

  go_process_terms_prgene <- go_annotations_qc %>%
    dplyr::filter(go_ontology == "P") %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(
      num_ids_process = length(unique(go_id)),
      .groups = "drop"
    )

  go_terms <- as.data.frame(
    dplyr::full_join(go_function_terms_prgene,
                     go_process_terms_prgene,
                     by = "symbol") %>%
      dplyr::filter(symbol != "") %>%
      dplyr::mutate(num_ids_process =
                      dplyr::if_else(
                        is.na(num_ids_process),as.integer(0),
                        as.integer(num_ids_process)
                      )
      ) %>%
      dplyr::mutate(num_ids_function =
                      dplyr::if_else(
                        is.na(num_ids_function),as.integer(0),
                        as.integer(num_ids_function)
                      )
      ) %>%
      dplyr::mutate(num_go_terms = num_ids_function +
                      num_ids_process) %>%

      dplyr::select(symbol, num_go_terms)
  )

  return(go_terms)
}
assign_unknown_function_rank <- function(
  gene_xref = NULL){

  gene_xref <- gene_xref %>%
    dplyr::mutate(
      unknown_function_rank =
        dplyr::case_when(
          stringr::str_detect(tolower(name),
                              "uncharacterized|open reading frame") &
            is.na(gene_summary) &
            num_go_terms == 0 ~ as.integer(1),
          stringr::str_detect(tolower(name),
                              "uncharacterized|open reading frame") &
            is.na(gene_summary) &
            num_go_terms > 0 &
            num_go_terms <= 3 ~ as.integer(2),
          is.na(gene_summary) &
            num_go_terms == 0 &
            !stringr::str_detect(tolower(name),
                                 "uncharacterized|open reading frame") ~ as.integer(3),
          is.na(gene_summary) &
            num_go_terms == 1 &
            !stringr::str_detect(tolower(name),
                                 "uncharacterized|open reading frame") ~ as.integer(4),
          !is.na(gene_summary) &
            num_go_terms == 0 &
            !stringr::str_detect(tolower(name),
                                 "uncharacterized|open reading frame") ~ as.integer(5),
          is.na(gene_summary) &
            num_go_terms == 2 &
            !stringr::str_detect(tolower(name),
                                 "uncharacterized|open reading frame") ~ as.integer(6),
          TRUE ~ as.integer(NA)
        )
    )

  return(gene_xref)

}

remove_duplicate_ensembl_genes <- function(
  gene_xref = NULL){

  ## Trick to remove duplicates
  ensembl2entrez <- gene_xref %>%
    dplyr::select(entrezgene, ensembl_gene_id) %>%
    dplyr::filter(!is.na(ensembl_gene_id))

  duplicates <- ensembl2entrez %>%
    dplyr::group_by(entrezgene) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::inner_join(
      dplyr::filter(gene_xref, !is.na(ensembl_gene_id)),
      by = "entrezgene")

  duplicate_entrez <- unique(duplicates$entrezgene)
  all_ranked_dups <- data.frame()
  chosen_candidates <- data.frame()

  for(i in 1:length(duplicate_entrez)){
    dset <- duplicates %>%
      dplyr::filter(entrezgene == duplicate_entrez[i])

    dset_ranked <- data.frame()
    symbol <- ""

    for(j in 1:nrow(dset)){
      e <- dset[j,]
      symbol <- e$symbol
      e$n <- NULL
      e$duplicate_rank <- 0

      for(m in c('SM_tractability_category',
                 'AB_tractability_category',
                 'gene_summary')){
        if(!is.na(e[,m])){
          e$duplicate_rank <- e$duplicate_rank + 1
          if(e[,m] != "Unknown"){
            e$duplicate_rank <- e$duplicate_rank + 1
          }

        }
      }
      if(e$cancer_max_rank > 0){
        e$duplicate_rank <- e$duplicate_rank + 1
      }

      dset_ranked <- dset_ranked %>%
        dplyr::bind_rows(e)

    }

    dset_ranked <- dset_ranked %>%
      dplyr::arrange(desc(duplicate_rank))

    tied_ranking <-
      length(unique(dset_ranked$duplicate_rank)) == 1

    if(tied_ranking){
      #cat('Not possible to rank - ', symbol, "\n")
      dset_ranked <- dset_ranked %>%
        dplyr::arrange(ensembl_gene_id)
    }

    dset_ranked$duplicate_rank <- NULL

    chosen_candidates <- chosen_candidates %>%
      dplyr::bind_rows(head(dset_ranked, 1))

  }

  to_remove <- duplicates %>%
    dplyr::anti_join(chosen_candidates, by = c("ensembl_gene_id"))

  gene_xref <- gene_xref %>%
    dplyr::anti_join(to_remove, by = c("ensembl_gene_id"))

  return(gene_xref)

}

generate_gene_xref_df <- function(
  basedir = NULL,
  gene_info = NULL,
  transcript_xref_db = NULL,
  ts_oncogene_annotations = NULL,
  opentarget_associations = NULL,
  uniprot_gene_summary = NULL,
  ncbi_gene_summary = NULL,
  cancerdrugdb = NULL,
  otdb = NULL,
  go_terms_pr_gene = NULL,
  update = T){

  rds_fname <- file.path(
    basedir, "data-raw", "gene_xref", "gene_xref.rds"
  )

  if(update == F & file.exists(rds_fname)){
    gene_xref <- readRDS(file = rds_fname)
  }


  invisible(assertthat::assert_that(!is.null(gene_info)))
  invisible(assertthat::assert_that(!is.null(transcript_xref_db)))
  invisible(assertthat::assert_that(!is.null(ts_oncogene_annotations)))
  invisible(assertthat::assert_that(!is.null(ncbi_gene_summary)))
  invisible(assertthat::assert_that(!is.null(uniprot_gene_summary)))
  invisible(assertthat::assert_that(!is.null(cancerdrugdb)))
  invisible(assertthat::assert_that(!is.null(otdb)))
  invisible(assertthat::assert_that(!is.null(go_terms_pr_gene)))
  invisible(assertthat::assert_that(!is.null(opentarget_associations)))

  invisible(assertthat::assert_that(is.data.frame(gene_info)))
  invisible(assertthat::assert_that(is.data.frame(transcript_xref_db)))
  invisible(assertthat::assert_that(is.data.frame(ts_oncogene_annotations)))
  invisible(assertthat::assert_that(is.data.frame(ncbi_gene_summary)))
  invisible(assertthat::assert_that(is.data.frame(uniprot_gene_summary)))
  invisible(assertthat::assert_that(typeof(cancerdrugdb) == "list"))
  invisible(assertthat::assert_that(typeof(otdb) == "list"))
  invisible(assertthat::assert_that(is.data.frame(go_terms_pr_gene)))
  invisible(assertthat::assert_that(is.data.frame(opentarget_associations)))

  invisible(assertthat::assert_that("drug_per_target" %in% names(cancerdrugdb)))
  invisible(assertthat::assert_that("early_phase" %in% names(cancerdrugdb[['drug_per_target']])))
  invisible(assertthat::assert_that("late_phase" %in% names(cancerdrugdb[['drug_per_target']])))
  invisible(assertthat::assert_that("max_site_rank" %in% names(otdb)))


  assertable::assert_colnames(
    otdb[['max_site_rank']],
    c('ensembl_gene_id','cancer_max_rank'), only_colnames = F,
    quiet = T)
  assertable::assert_colnames(
    gene_info, c('entrezgene', 'name', 'symbol_entrez',
                  'hgnc_id', 'gene_biotype'), only_colnames = F,
    quiet = T)
  assertable::assert_colnames(
    transcript_xref_db, c('entrezgene', 'property',
                 'value'), quiet = T)
  assertable::assert_colnames(
    uniprot_gene_summary,
      c('symbol','gene_summary_uniprot'), quiet = T)
  assertable::assert_colnames(
    ncbi_gene_summary,
    c('entrezgene','gene_summary_ncbi'), quiet = T)

  assertable::assert_colnames(
    go_terms_pr_gene,
    c('symbol','num_go_terms'), quiet = T)

  assertable::assert_colnames(
    ts_oncogene_annotations,
    c('entrezgene','tumor_suppressor','oncogene', 'cancer_driver'),
    quiet = T, only_colnames = F)

  assertable::assert_colnames(
    opentarget_associations,
    c('ensembl_gene_id','SM_tractability_category','SM_tractability_support',
      'AB_tractability_category','AB_tractability_support'),
    quiet = T, only_colnames = F)

  go2entrez <- transcript_xref_db %>%
    dplyr::filter(property == "symbol") %>%
    dplyr::rename(symbol = value) %>%
    dplyr::left_join(go_terms_pr_gene, by = "symbol") %>%
    dplyr::filter(!is.na(num_go_terms)) %>%
    dplyr::select(-c(property, symbol)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(entrezgene) %>%
    dplyr::summarise(num_go_terms = max(num_go_terms))

  upsummary2entrez <- transcript_xref_db %>%
    dplyr::filter(property == "symbol") %>%
    dplyr::rename(symbol = value) %>%
    dplyr::left_join(uniprot_gene_summary, by = "symbol") %>%
    dplyr::filter(!is.na(gene_summary_uniprot)) %>%
    dplyr::distinct() %>%
    dplyr::select(entrezgene, gene_summary_uniprot)

  ensembl2entrez <- transcript_xref_db %>%
    dplyr::filter(property == "ensembl_gene_id") %>%
    dplyr::rename(ensembl_gene_id = value) %>%
    dplyr::select(entrezgene, ensembl_gene_id) %>%
    dplyr::distinct()

  gene_xref <- gene_info %>%
    dplyr::select(symbol_entrez, hgnc_id, entrezgene,
                  name, gene_biotype) %>%
    dplyr::rename(symbol = symbol_entrez) %>%
    dplyr::left_join(ensembl2entrez, by = "entrezgene") %>%
    dplyr::left_join(go2entrez, by = "entrezgene") %>%
    dplyr::filter(!is.na(entrezgene) & !is.na(symbol)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      genename =
        paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
               entrezgene,
               "' target='_blank'>",name,"</a>")) %>%
    dplyr::left_join(opentarget_associations,
                     by = "ensembl_gene_id") %>%
    dplyr::left_join(ts_oncogene_annotations,
                     by = "entrezgene") %>%
    dplyr::left_join(upsummary2entrez,
                     by = "entrezgene") %>%
    dplyr::left_join(ncbi_gene_summary,
                     by = "entrezgene") %>%
    dplyr::left_join(cancerdrugdb[['drug_per_target']][['early_phase']],
                     by = "entrezgene") %>%
    dplyr::left_join(cancerdrugdb[['drug_per_target']][['late_phase']],
                     by = "entrezgene") %>%
    dplyr::mutate(tumor_suppressor = dplyr::if_else(
      is.na(tumor_suppressor),
      FALSE,
      as.logical(tumor_suppressor))) %>%
    dplyr::mutate(oncogene = dplyr::if_else(
      is.na(oncogene),
      FALSE,
      as.logical(oncogene))) %>%
    dplyr::mutate(cancer_driver = dplyr::if_else(
      is.na(cancer_driver),
      FALSE,
      as.logical(cancer_driver))) %>%
    dplyr::mutate(gene_summary = dplyr::case_when(
      !is.na(gene_summary_ncbi) &
        !is.na(gene_summary_uniprot) ~
        paste0(gene_summary_ncbi," ",
               gene_summary_uniprot),
      !is.na(gene_summary_ncbi) &
        is.na(gene_summary_uniprot) ~
        gene_summary_ncbi,
      is.na(gene_summary_ncbi) &
        !is.na(gene_summary_uniprot) ~
        gene_summary_uniprot,
      TRUE ~ as.character(NA)
    )) %>%
    dplyr::select(-c(gene_summary_ncbi, gene_summary_uniprot)) %>%

    dplyr::mutate(has_gene_summary = dplyr::if_else(
      !is.na(gene_summary),TRUE,FALSE
    )) %>%
    dplyr::left_join(otdb$max_site_rank,
                     by = "ensembl_gene_id") %>%
    dplyr::mutate(cancer_max_rank = dplyr::if_else(
      is.na(cancer_max_rank), 0, as.numeric(cancer_max_rank)
    )) %>%
    assign_unknown_function_rank() %>%
    remove_duplicate_ensembl_genes()

  saveRDS(gene_xref, file = rds_fname)
  return(gene_xref)

}

get_tissue_celltype_specificity <- function(
  basedir = NULL){

    ## https://www.proteinatlas.org/about/assays+annotation#gtex_rna
    ##
    ## Transcript expression levels summarized per gene in 36 tissues
    ## based on RNA-seq. The tab-separated file includes Ensembl
    ## gene identifier ("Gene"), analysed sample ("Tissue"),
    ## transcripts per million ("TPM"), protein-transcripts
    ## per million ("pTPM") and normalized expression ("NX").

    tissue_cell_expr <- list()
    tissue_cell_expr[['tissue']] <- list()
    tissue_cell_expr[['tissue']][['expr_df']] <- data.frame()
    tissue_cell_expr[['tissue']][['unit']] <- NULL
    tissue_cell_expr[['tissue']][['te_SE']] <- NULL
    tissue_cell_expr[['tissue']][['te_df']] <- data.frame()

    ## https://www.proteinatlas.org/about/assays+annotation#singlecell_rna
    ##
    ## Transcript expression levels summarized per gene in 51 cell
    ## types from 13 tissues. The tab-separated file includes Ensembl
    ## gene identifier ("Gene"), gene name ("Gene name"), analysed
    ## sample ("Cell type") and normalized expresion ("NX").

    tissue_cell_expr[['single_cell']] <- list()
    tissue_cell_expr[['single_cell']][['expr_df']] <- data.frame()
    tissue_cell_expr[['single_cell']][['unit']] <- NULL
    tissue_cell_expr[['single_cell']][['te_SE']] <- NULL
    tissue_cell_expr[['single_cell']][['te_df']] <- data.frame()

    for(t in c("rna_tissue_gtex",
               "rna_single_cell_type")){

      local_hpa_file <- file.path(
        basedir, "data-raw", "hpa",
        paste0(t, ".tsv.zip")
      )

      if(!file.exists(local_hpa_file)){
        download.file(
          url = paste0("https://www.proteinatlas.org/download/",
                       t, ".tsv.zip"),
          destfile = local_hpa_file
        )
      }

      ##set max number of tissues/cell_types for
      ##determination of groups in group-enriched genes
      max_types_in_group <- 10
      if(t == "rna_tissue_gtex"){
        max_types_in_group <- 5
      }

      data <- readr::read_tsv(
        local_hpa_file,
        show_col_types = F,
        col_names = T)

      num_tissues_per_gene <- data %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(n = dplyr::n()) %>%
        dplyr::filter(n == 1)

      data <- data %>%
        dplyr::anti_join(num_tissues_per_gene,
                         by = "Gene")

      if(t == "rna_tissue_gtex"){
        data <- as.data.frame(
          data[, c(1,3,6)] %>%
            magrittr::set_colnames(c("ensembl_gene_id",
                                     "category","exp1")) %>%
            dplyr::mutate(
              category =
                stringr::str_replace_all(category," |, ","_")) %>%
            dplyr::group_by(ensembl_gene_id, category) %>%
            dplyr::summarise(exp = mean(exp1), .groups = "drop") %>%
            tidyr::pivot_wider(names_from = category,
                               values_from = exp)
        )
      }else{
        data <- as.data.frame(
          data[, c(1,3,4)] %>%
            magrittr::set_colnames(c("ensembl_gene_id",
                                     "category","exp1")) %>%
            dplyr::mutate(
              category =
                stringr::str_replace_all(category," |, ","_")) %>%
            dplyr::group_by(ensembl_gene_id, category) %>%
            dplyr::summarise(exp = mean(exp1), .groups = "drop") %>%
            tidyr::pivot_wider(names_from = category,
                               values_from = exp)
        )
      }
      rownames(data) <- data$ensembl_gene_id
      data$ensembl_gene_id <- NULL
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(
          as.matrix(data)),
        rowData = row.names(data),
        colData = colnames(data)
      )
      te_gene_retrieval_se <-
        TissueEnrich::teGeneRetrieval(
          expressionData = se,
          foldChangeThreshold = 4,
          maxNumberOfTissues = max_types_in_group)

      if(t == "rna_tissue_gtex"){
        tissue_cell_expr[['tissue']][['expr_df']] <- data
        tissue_cell_expr[['tissue']][['te_SE']] <-
          te_gene_retrieval_se
        tissue_cell_expr[['tissue']][['unit']] <- 'nTPM'
        tissue_cell_expr[['tissue']][['te_df']] <-
          as.data.frame(
            as.data.frame(
              SummarizedExperiment::assay(
                tissue_cell_expr[['tissue']][['te_SE']])
            ) %>%
              dplyr::rename(ensembl_gene_id = Gene,
                            category = Group) %>%
              dplyr::group_by(ensembl_gene_id,
                              category) %>%
              dplyr::summarise(
                tissue = paste(unique(Tissue), collapse=", "),
                .groups = "drop")
          ) %>%
          dplyr::mutate(
            category = dplyr::case_when(
              category == "Expressed-In-All" ~ "Low tissue specificity",
              category == "Group-Enriched" ~ "Group enriched",
              category == "Tissue-Enhanced" ~ "Tissue enhanced",
              category == "Tissue-Enriched" ~ "Tissue enriched",
              category == "Not-Expressed" ~ "Not detected",
              TRUE ~ as.character("Mixed"))
          ) %>%
          dplyr::mutate(
            category = factor(category,
                              levels = c("Tissue enriched",
                                         "Group enriched",
                                         "Tissue enhanced",
                                         "Mixed",
                                         "Low tissue specificity",
                                         "Not detected"))
          )


      }else{
        tissue_cell_expr[['single_cell']][['expr_df']] <- data
        tissue_cell_expr[['single_cell']][['unit']] <- 'nTPM'
        tissue_cell_expr[['single_cell']][['te_SE']] <-
          te_gene_retrieval_se
        tissue_cell_expr[['single_cell']][['te_df']] <-
          as.data.frame(
            as.data.frame(
              SummarizedExperiment::assay(
                tissue_cell_expr[['single_cell']][['te_SE']])
            ) %>%
              dplyr::rename(ensembl_gene_id = Gene,
                            category = Group) %>%
              dplyr::group_by(ensembl_gene_id,
                              category) %>%
              dplyr::summarise(
                cell_type = paste(unique(Tissue), collapse=", "),
                .groups = "drop")
          ) %>%
          dplyr::mutate(
            category = dplyr::case_when(
              category == "Expressed-In-All" ~ "Low cell type specificity",
              category == "Group-Enriched" ~ "Group enriched",
              category == "Tissue-Enhanced" ~ "Cell type enhanced",
              category == "Tissue-Enriched" ~ "Cell type enriched",
              category == "Not-Expressed" ~ "Not detected",
              TRUE ~ as.character("Mixed"))
          ) %>%
          dplyr::mutate(
            category =
              factor(category,
                     levels = c("Cell type enriched",
                                "Group enriched",
                                "Cell type enhanced",
                                "Mixed",
                                "Low cell type specificity",
                                "Not detected"))
          )
      }

    }
    return(tissue_cell_expr)

  }

quantify_gene_cancer_relevance <- function(
  ot_associations = NULL){


  phenotype_cancer_efo <-
    oncoPhenoMap::auxiliary_maps$umls$concept %>%
    dplyr::filter(main_term == T) %>%
    dplyr::select(cui, cui_name) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(
      dplyr::select(oncoPhenoMap::oncotree_expanded_full,
                    efo_id, cui,
                    cui_name, primary_site),
      by = c("cui", "cui_name")) %>%
    dplyr::rename(disease_efo_id = efo_id) %>%
    dplyr::filter(!is.na(disease_efo_id)) %>%
    dplyr::mutate(cancer_phenotype = TRUE) %>%
    dplyr::distinct()


  otdb_tmp <- as.data.frame(
    dplyr::select(ot_associations,
                  ot_association,
                  #symbol,
                  ensembl_gene_id) %>%
      dplyr::filter(!is.na(ot_association)) %>%
      tidyr::separate_rows(ot_association, sep="&") %>%
      tidyr::separate(
        ot_association,
        sep=":",
        c('disease_efo_id',
          'direct_ot_association',
          'ot_association_score')) %>%
      dplyr::mutate(ot_association_score =
                      as.numeric(ot_association_score)) %>%
      dplyr::filter(!is.na(disease_efo_id)) %>%

      dplyr::mutate(
        disease_efo_id = stringr::str_replace(disease_efo_id,"_",":")) %>%
      dplyr::left_join(
        dplyr::select(oncoPhenoMap::auxiliary_maps$efo$efo2name,
                      efo_id, efo_name),
        by = c("disease_efo_id" = "efo_id")) %>%
      ## exclude non-specific cancer associations
      dplyr::filter(
        !stringr::str_detect(
          tolower(efo_name),
          "^(((urogenital|mixed|papillary|hereditary|familial|susceptibility|syndrome|small cell|clear cell|squamous cell|digestive system|endocrine|metastatic|invasive|pharynx|vascular|intestinal|cystic|mucinous|epithelial|benign) )?(neoplasm|carcinoma|adenocarcinoma|cancer))$")
      ) %>%
      dplyr::filter(
        !stringr::str_detect(
          tolower(efo_name)," neoplasm$"
        )
      ) %>%
      dplyr::left_join(
        dplyr::select(phenotype_cancer_efo,
                      disease_efo_id,
                      primary_site,
                      cancer_phenotype),
        by=c("disease_efo_id")) %>%
      dplyr::mutate(
        cancer_phenotype =
          dplyr::if_else(
            is.na(cancer_phenotype) &
              stringr::str_detect(
                tolower(efo_name),
                "carcinoma|tumor|cancer|neoplasm|melanoma|myeloma|hemangioma|astrocytoma|leiomyoma|leukemia|lymphoma|barrett|glioma|sarcoma|blastoma|teratoma|seminoma"),
            TRUE, as.logical(cancer_phenotype))) %>%
      dplyr::mutate(
        cancer_phenotype =
          dplyr::if_else(
            is.na(cancer_phenotype) &
              stringr::str_detect(
                tolower(efo_name),
                "hereditary|familial|susceptibility|syndrome") &
              stringr::str_detect(
                tolower(efo_name),
                "tumor|cancer|lynch|carcinoma|sarcoma|melanoma|blastoma"),
            TRUE, as.logical(cancer_phenotype))) %>%
      dplyr::mutate(
        ot_link = paste0(
          "<a href='https://platform.opentargets.org/evidence/",
          ensembl_gene_id,"/",
          stringr::str_replace(disease_efo_id,":","_"),
          "' target=\"_blank\">", stringr::str_to_title(efo_name),"</a>")) %>%
      dplyr::distinct() %>%
      dplyr::distinct()
  )

  set1 <- otdb_tmp %>%
    dplyr::filter(!is.na(primary_site)) %>%
    dplyr::distinct()
  set2 <- otdb_tmp %>%
    dplyr::filter(is.na(primary_site) & cancer_phenotype == T) %>%
    dplyr::select(-c(efo_name, primary_site)) %>%
    dplyr::anti_join(
      dplyr::select(set1, disease_efo_id, ensembl_gene_id),
      by = c("disease_efo_id","ensembl_gene_id")) %>%
    dplyr::inner_join(
      dplyr::select(oncoPhenoMap::auxiliary_maps$efo$efo2name,
                    efo_id, efo_name, primary_site),
      by = c("disease_efo_id" = "efo_id")) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(primary_site))

  set12 <- dplyr::bind_rows(set1, set2)

  set3 <- otdb_tmp %>%
    dplyr::anti_join(
      dplyr::select(set12, disease_efo_id, ensembl_gene_id),
      by = c("disease_efo_id","ensembl_gene_id")) %>%
    dplyr::distinct()

  otdb_all <- set12 %>%
    dplyr::bind_rows(set3) %>%
    dplyr::arrange(primary_site, ensembl_gene_id,
                   desc(ot_association_score))

  otdb_tissue_scores <- as.data.frame(
    otdb_all %>%
      dplyr::filter(!is.na(primary_site)) %>%
      dplyr::group_by(primary_site, ensembl_gene_id) %>%
      dplyr::summarise(
        tissue_assoc_score =
          round(sum(ot_association_score), digits = 3),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        primary_site, desc(tissue_assoc_score))
  )


  otdb_site_rank <- data.frame()
  for(s in unique(otdb_tissue_scores$primary_site)){
    tmp <- otdb_tissue_scores %>%
      dplyr::filter(primary_site == s) %>%
      dplyr::mutate(
        tissue_assoc_rank =
          round(dplyr::percent_rank(tissue_assoc_score),
                digits = 3)
      )
    otdb_site_rank <- otdb_site_rank %>%
      dplyr::bind_rows(tmp)
  }

  otdb_max_site_rank <- as.data.frame(
    otdb_site_rank %>%
      dplyr::arrange(desc(tissue_assoc_rank)) %>%
      dplyr::select(ensembl_gene_id, tissue_assoc_rank) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(
        cancer_max_rank = max(tissue_assoc_rank),
        .groups = "drop")
  )

  otdb <- list()
  otdb$all <- otdb_all
  otdb$site_rank <- otdb_site_rank
  otdb$max_site_rank <- otdb_max_site_rank

  return(otdb)

}

get_ligand_receptors <- function(
  basedir = NULL,
  keggdb = NULL,
  update = F){

  rds_fname <- file.path(
    basedir, "data-raw",
    "cellchatdb", "cellchatdb_interactions.rds")

  if(update == F & file.exists(rds_fname)){
    ligandreceptordb <- readRDS(file = rds_fname)
    return(ligandreceptordb)
  }

  ligand_receptor_db <- CellChat::CellChatDB.human$interaction %>%
    magrittr::set_rownames(NULL) %>%
    dplyr::rename(interaction_members = interaction_name) %>%
    dplyr::rename(interaction_name = interaction_name_2) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC: 4393358", "PMID: 25772309"
    )) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC4571854", "PMID: 26124272"
    )) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC2194005", "PMID: 12208882"
    )) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC2194217", "PMID: 14530377"
    )) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC1237098", "PMID: 16093349"
    )) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,"PMC2431087", "PMID: 18250165"
    )) %>%
    dplyr::mutate(evidence = stringr::str_squish(
      stringr::str_replace_all(
        evidence,
        ",","; "))
    ) %>%
    dplyr::mutate(evidence = stringr::str_replace_all(
      evidence,
      " ","")) %>%
    dplyr::mutate(interaction_id = dplyr::row_number()) %>%
    dplyr::select(interaction_id, interaction_name,
                  annotation, pathway_name, dplyr::everything())


  ligand_receptor_kegg_support <- ligand_receptor_db %>%
    dplyr::select(interaction_id, evidence) %>%
    tidyr::separate_rows(evidence, sep=";") %>%
    dplyr::filter(stringr::str_detect(evidence,"KEGG")) %>%
    dplyr::mutate(evidence = stringr::str_replace(
      evidence, "KEGG:",""
    )) %>%
    dplyr::left_join(keggdb$TERM2NAME, by = c("evidence" = "standard_name")) %>%
    dplyr::mutate(ligand_receptor_kegg_support = paste0(
      "<a href='https://www.genome.jp/pathway/",
      stringr::str_replace(evidence,"hsa","map"),
      "' target='_blank'>",name,"</a>")) %>%
    dplyr::select(-evidence)

  ligand_receptor_literature <- as.data.frame(
    ligand_receptor_db %>%
      dplyr::select(interaction_id, evidence) %>%
      tidyr::separate_rows(evidence, sep=";") %>%
      dplyr::filter(stringr::str_detect(evidence,"PMID")) %>%
      dplyr::mutate(evidence = stringr::str_replace(
        evidence, "PMID:",""
      )) %>%
      dplyr::rename(pmid = evidence)
  )

  ligand_receptor_citations <- get_citations_pubmed(
    pmid = unique(ligand_receptor_literature$pmid))

  ligand_receptor_literature <- as.data.frame(
    ligand_receptor_literature %>%
      dplyr::mutate(pmid = as.integer(pmid)) %>%
      dplyr::left_join(ligand_receptor_citations, by = c("pmid" = "pmid")) %>%
      dplyr::group_by(interaction_id) %>%
      dplyr::summarise(
        ligand_receptor_literature_support = paste(
          link, collapse= ","),
        ligand_receptor_literature = paste(
          citation, collapse="|"
        )
      )
  )

  ligand_receptor_db_final <- ligand_receptor_db %>%
    dplyr::left_join(dplyr::select(
      ligand_receptor_kegg_support, interaction_id,
      ligand_receptor_kegg_support), by = "interaction_id") %>%
    dplyr::left_join(dplyr::select(
      ligand_receptor_literature, interaction_id,
      ligand_receptor_literature_support)) %>%
    dplyr::mutate(literature_support = dplyr::if_else(
      is.na(ligand_receptor_kegg_support),
      as.character(ligand_receptor_literature_support),
      paste(ligand_receptor_kegg_support,
            ligand_receptor_literature_support,
            sep = ", ")
    )) %>%
    dplyr::mutate(literature_support = stringr::str_replace_all(
      literature_support,",?NA$",""
    )) %>%
    dplyr::select(-c(ligand_receptor_literature_support,
                     ligand_receptor_kegg_support,
                     evidence))

  interaction2ligands <- as.data.frame(
    ligand_receptor_db_final %>%
      dplyr::select(interaction_id, ligand) %>%
      dplyr::rename(symbol = ligand) %>%
      dplyr::mutate(class = "ligand")
  )

  interaction2receptor <- as.data.frame(
    ligand_receptor_db_final %>%
      dplyr::select(interaction_id, receptor) %>%
      dplyr::rename(symbol = receptor) %>%
      tidyr::separate_rows(symbol, sep = "_") %>%
      dplyr::mutate(symbol = dplyr::if_else(
        stringr::str_detect(symbol,"^(R2|TGFbR2)$"),
        "TGFBR2",
        as.character(symbol)
      )) %>%
      dplyr::mutate(class = "receptor") %>%
      dplyr::mutate(symbol = stringr::str_replace(
        symbol,"b","B"
      ))
  )

  ligandreceptordb <- list()
  ligandreceptordb[['db']] <- ligand_receptor_db_final
  ligandreceptordb[['xref']] <- interaction2ligands %>%
    dplyr::bind_rows(interaction2receptor)


  saveRDS(ligandreceptordb,
          file = rds_fname)
  return(ligandreceptordb)
}

get_hpa_associations <- function(
  basedir = NULL,
  gene_xref = NULL,
  update = F){

  assertable::assert_colnames(
    gene_xref,
    c("symbol", "ensembl_gene_id"),
    only_colnames = F,
    quiet = T
  )

  rds_fname <- file.path(basedir, "data-raw", "hpa", "hpa.rds")

  if(update == F & file.exists(rds_fname)){
    hpa <- readRDS(file = rds_fname)
    return(hpa)
  }

  variables_json <- c('RNA tissue specificity',
                      'RNA tissue distribution',
                      'RNA tissue specific NX',
                      'RNA single cell type specificity',
                      'RNA single cell type distribution',
                      'RNA single cell type specific NX',
                      'RNA cancer specificity',
                      'RNA cancer distribution',
                      'RNA cancer specific FPKM',
                      'RNA cell line specificity',
                      'RNA cell line distribution',
                      'RNA cell line specific NX',
                      'Pathology prognostics - Breast cancer',
                      'Pathology prognostics - Cervical cancer',
                      'Pathology prognostics - Colorectal cancer',
                      'Pathology prognostics - Endometrial cancer',
                      'Pathology prognostics - Glioma',
                      'Pathology prognostics - Head and neck cancer',
                      'Pathology prognostics - Liver cancer',
                      'Pathology prognostics - Lung cancer',
                      'Pathology prognostics - Melanoma',
                      'Pathology prognostics - Ovarian cancer',
                      'Pathology prognostics - Pancreatic cancer',
                      'Pathology prognostics - Prostate cancer',
                      'Pathology prognostics - Renal cancer',
                      'Pathology prognostics - Stomach cancer',
                      'Pathology prognostics - Testis cancer',
                      'Pathology prognostics - Thyroid cancer',
                      'Pathology prognostics - Urothelial cancer',
                      'Antibody RRID')

  variables_df <- c('rna_tissue_specificity',
                    'rna_tissue_distribution',
                    'rna_tissue_specific_NX',
                    'rna_single_cell_type_specificity',
                    'rna_single_cell_type_distribution',
                    'rna_single_cell_type_specificity NX',
                    'rna_cancer_specificity',
                    'rna_cancer_distribution',
                    'rna_cancer_specific_FPKM',
                    'rna_cell_line_specificity',
                    'rna_cell_line_distribution',
                    'rna_cell_line_specific_NX',
                    'pathology_prognostics_breast_cancer',
                    'pathology_prognostics_cervical_cancer',
                    'pathology_prognostics_colorectal_cancer',
                    'pathology_prognostics_endometrial_cancer',
                    'pathology_prognostics_glioma',
                    'pathology_prognostics_head_and_neck_cancer',
                    'pathology_prognostics_liver_cancer',
                    'pathology_prognostics_lung_cancer',
                    'pathology_prognostics_melanoma',
                    'pathology_prognostics_ovarian_cancer',
                    'pathology_prognostics_pancreatic_cancer',
                    'pathology_prognostics_prostate_cancer',
                    'pathology_prognostics_renal_cancer',
                    'pathology_prognostics_stomach_cancer',
                    'pathology_prognostics_testis_cancer',
                    'pathology_prognostics_thyroid_cancer',
                    'pathology_prognostics_urothelial_cancer',
                    'antibody_rrid')

  for (i in 1:nrow(gene_xref)) {
    g <- gene_xref[i,]

    local_json <-
      file.path("data-raw", "hpa", "json",
                paste0(g$ensembl_gene_id,".json"))

    if(!file.exists(local_json)){
      protein_atlas_url <-
        paste0("https://www.proteinatlas.org/",
               g$ensembl_gene_id,".json")
      if(RCurl::url.exists(protein_atlas_url)){
        download.file(
          protein_atlas_url,
          destfile = local_json,
          quiet = T)
        Sys.sleep(1)
      }
    }

    if(file.exists(local_json)){
      pa_data <- jsonlite::fromJSON(txt = local_json)

      for(j in 1:length(variables_json)){
        var_name_json <- variables_json[j]
        var_name_df <- variables_df[j]
        if(!is.null(pa_data[[var_name_json]])){
          if(!startsWith(var_name_df, "pathology") &
             !stringr::str_detect(var_name_df,"(NX|FPKM|rrid)$")){
            df <- data.frame(
              'ensembl_gene_id' = g$ensembl_gene_id,
              'property' = var_name_df,
              'value' = pa_data[[var_name_json]],
              stringsAsFactors = F)

          }else{
            if(startsWith(var_name_df, "pathology")
               | endsWith(var_name_df,"rrid")){
              value <- paste(unlist(pa_data[[var_name_json]]),
                             collapse="|")
              if(stringr::str_detect(value,"TRUE") |
                 var_name_df == "antibody_rrid"){
                df <- data.frame(
                  'ensembl_gene_id' = g$ensembl_gene_id,
                  'property' = var_name_df,
                  'value' = value,
                  stringsAsFactors = F)
              }
            }else{
              df <- data.frame(
                'ensembl_gene_id' = g$ensembl_gene_id,
                'property' = var_name_df,
                'value' = paste(sort(names(pa_data[[var_name_json]])),
                                collapse=", "),
                stringsAsFactors = F)
            }
          }
          if(NROW(df) > 0){
            hpa <- dplyr::bind_rows(hpa, df)
            df <- data.frame()
          }
        }
      }
    }

  }

  saveRDS(hpa, file = rds_fname)
  return(hpa)


}


get_tcga_db <- function(
  basedir = NULL,
  gene_xref = NULL,
  update = F){

  rds_fname <- file.path(
    basedir, "data-raw", "tcga", "tcgadb.rds")

  co_expression_tsv <- file.path(
    basedir, "data-raw", "tcga",
    "co_expression_strong_moderate.release31_20211029.tsv.gz")

  tcga_clinical_rds <- file.path(
    basedir, "data-raw", "tcga",
    "tcga_clinical.rds"
  )

  tcga_aberration_rds <- file.path(
    basedir, "data-raw", "tcga",
    "tcga_gene_aberration_stats.rds"
  )

  maf_codes_tsv <- file.path(
    basedir, "data-raw", "maf_codes.tsv"
  )
  maf_path <- file.path(basedir, "data-raw", "tcga")


  recurrent_variants_rds <- file.path(
    basedir, "data-raw", "tcga",
    "tcga_recurrent_coding_gvanno_grch38.tsv.gz"
  )

  pfam_domains_tsv <- file.path(
    basedir, "data-raw", "pfam",
    "pfam.domains.tsv.gz"
  )

  if(update == F & file.exists(rds_fname)){
    tcgadb <- readRDS(file = rds_fname)
    return(tcgadb)
  }


  tcga_clinical <-
    readRDS(file = tcga_clinical_rds)
  tcga_aberration_stats <-
    readRDS(file = tcga_aberration_rds) %>%
    dplyr::filter(
      clinical_strata == "site" |
        (clinical_strata == "site_diagnosis" &
           percent_mutated >= 1))

  tcga_aberration_stats$genomic_strata <- NULL
  tcga_aberration_stats$entrezgene <- NULL
  tcga_aberration_stats$consensus_calls <- NULL
  tcga_aberration_stats$fp_driver_gene <- NULL

  site_code <- tcga_aberration_stats %>%
    dplyr::select(primary_site) %>%
    dplyr::distinct() %>%
    dplyr::arrange(primary_site) %>%
    dplyr::mutate(site_code = dplyr::row_number())

  diagnosis_code <- tcga_aberration_stats %>%
    dplyr::select(primary_diagnosis_very_simplified) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      primary_diagnosis = primary_diagnosis_very_simplified) %>%
    dplyr::arrange(primary_diagnosis) %>%
    dplyr::mutate(
      diagnosis_code = dplyr::row_number())

  clinical_strata_code <- tcga_aberration_stats %>%
    dplyr::select(clinical_strata) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      clinical_strata_code = dplyr::row_number())


  tcga_aberration_stats <- tcga_aberration_stats %>%
    dplyr::rename(primary_diagnosis = primary_diagnosis_very_simplified) %>%
    dplyr::left_join(clinical_strata_code, by = "clinical_strata") %>%
    dplyr::left_join(diagnosis_code, by = "primary_diagnosis") %>%
    dplyr::left_join(site_code, by = "primary_site") %>%
    dplyr::select(-c(primary_site, primary_diagnosis,
                     clinical_strata))

  maf_codes <- read.table(file = maf_codes_tsv,
                          header = T, sep = "\t", quote ="")
  i <- 1
  while(i <= nrow(maf_codes)){
    primary_site <- maf_codes[i,]$primary_site
    maf_code <- maf_codes[i,]$code
    maf_file <- paste0(maf_path,"/tcga_mutation_grch38_release31_20211029.",maf_code,"_0.maf.gz")
    if(file.exists(maf_file)){
      tmp <- read.table(gzfile(maf_file), quote="", header = T, stringsAsFactors = F, sep="\t", comment.char="#")
      tmp$primary_site <- NULL
      tmp$site_diagnosis_code <- NULL
      tmp$Tumor_Sample_Barcode <- stringr::str_replace(tmp$Tumor_Sample_Barcode,"-[0-9][0-9][A-Z]$","")

      clinical <- tcga_clinical %>% dplyr::filter(primary_site == primary_site) %>%
        dplyr::select(bcr_patient_barcode, primary_diagnosis_very_simplified,
                      MSI_status, Gleason_score, ER_status,
                      PR_status, HER2_status,
                      pancan_subtype_selected) %>%
        dplyr::rename(Diagnosis = primary_diagnosis_very_simplified,
                      Tumor_Sample_Barcode = bcr_patient_barcode,
                      PanCancer_subtype = pancan_subtype_selected) %>%
        dplyr::semi_join(tmp, by = "Tumor_Sample_Barcode")

      write.table(tmp, file="tmp.maf", quote=F, row.names = F,
                  col.names = T, sep ="\t")
      system('gzip -f tmp.maf', intern = T)
      maf <- maftools::read.maf("tmp.maf.gz", verbose = F, clinicalData = clinical)
      saveRDS(maf, file = file.path(basedir, "data-raw", "maf", paste0(maf_code,".maf.rds")))
    }
    i <- i + 1
    cat(primary_site,'\n')

  }
  system('rm -f tmp.maf.gz')


  pfam_domains <- as.data.frame(
    readr::read_tsv(
      pfam_domains_tsv,

      show_col_types = F
    )) %>%
    dplyr::rename(PFAM_ID = pfam_id, PROTEIN_DOMAIN = url) %>%
    dplyr::rename(PFAM_DOMAIN_NAME = name) %>%
    dplyr::select(PFAM_ID, PFAM_DOMAIN_NAME) %>%
    dplyr::mutate(PFAM_ID = stringr::str_replace(
      PFAM_ID,"\\.[0-9]{1,}$","")
    )

  ##TCGA recurrent SNVs/InDels
  recurrent_tcga_variants <- as.data.frame(readr::read_tsv(
    file = recurrent_variants_rds,
    skip = 1, na = c("."), show_col_types = F) %>%
      dplyr::select(TCGA_SITE_RECURRENCE,
                    CHROM, POS, REF, ALT,
                    TCGA_TOTAL_RECURRENCE,
                    PFAM_DOMAIN, LoF,
                    MUTATION_HOTSPOT,
                    HGVSp_short,
                    ENSEMBL_TRANSCRIPT_ID,
                    SYMBOL,
                    COSMIC_MUTATION_ID,
                    Consequence) %>%
      dplyr::mutate(VAR_ID = paste(
        CHROM, POS, REF, ALT, sep = "_")
      ) %>%
      dplyr::rename(PFAM_ID = PFAM_DOMAIN) %>%
      dplyr::mutate(LOSS_OF_FUNCTION = FALSE) %>%
      dplyr::mutate(LOSS_OF_FUNCTION = dplyr::if_else(
        !is.na(LoF) & LoF == "HC",
        as.logical(TRUE),
        as.logical(LOSS_OF_FUNCTION),
      )) %>%
      dplyr::rename(CONSEQUENCE = Consequence,
                    TOTAL_RECURRENCE = TCGA_TOTAL_RECURRENCE,
                    PROTEIN_CHANGE = HGVSp_short) %>%
      dplyr::select(SYMBOL,
                    VAR_ID,
                    CONSEQUENCE,
                    PROTEIN_CHANGE,
                    PFAM_ID,
                    MUTATION_HOTSPOT,
                    LOSS_OF_FUNCTION,
                    ENSEMBL_TRANSCRIPT_ID,
                    COSMIC_MUTATION_ID,
                    TCGA_SITE_RECURRENCE,
                    TOTAL_RECURRENCE) %>%
      tidyr::separate_rows(TCGA_SITE_RECURRENCE, sep=",") %>%
      tidyr::separate(TCGA_SITE_RECURRENCE, into =
                        c("PRIMARY_SITE","SITE_RECURRENCE", "TCGA_SAMPLES"),
                      sep = ":",
                      remove = T) %>%
      dplyr::select(-TCGA_SAMPLES) %>%
      dplyr::distinct() %>%
      dplyr::mutate(PRIMARY_SITE = dplyr::case_when(
        PRIMARY_SITE == "CNS_Brain" ~ "CNS/Brain",
        PRIMARY_SITE == "Colon_Rectum" ~ "Colon/Rectum",
        PRIMARY_SITE == "Head_and_Neck" ~ "Head and Neck",
        PRIMARY_SITE == "Esophagus_Stomach" ~ "Esophagus/Stomach",
        PRIMARY_SITE == "Bladder_Urinary_Tract" ~ "Bladder/Urinary Tract",
        PRIMARY_SITE == "Adrenal_Gland" ~ "Adrenal Gland",
        PRIMARY_SITE == "Biliary_Tract" ~ "Biliary Tract",
        PRIMARY_SITE == "Ovary_Fallopian_Tube" ~ "Ovary/Fallopian Tube",
        PRIMARY_SITE == "Soft_Tissue" ~ "Soft Tissue",
        TRUE ~ as.character(PRIMARY_SITE))
      ) %>%
      dplyr::filter(!stringr::str_detect(
        CONSEQUENCE,"^(intron|intergenic|mature|non_coding|synonymous|upstream|downstream|3_prime|5_prime)"))
  )

  ##TCGA co-expression data
  raw_coexpression <-
    readr::read_tsv(co_expression_tsv,
                    col_names = c("symbol_A","symbol_B",
                                  "r","p_value","tumor"),
                    show_col_types = F)
  co_expression_genes1 <- raw_coexpression %>%
    dplyr::rename(symbol = symbol_A, symbol_partner = symbol_B) %>%
    dplyr::left_join(
      dplyr::select(gene_xref, symbol, tumor_suppressor,
                    oncogene, cancer_driver), by = "symbol") %>%
    dplyr::filter(tumor_suppressor == T |
                    oncogene == T |
                    cancer_driver == T)

  co_expression_genes2 <- raw_coexpression %>%
    dplyr::rename(symbol = symbol_B, symbol_partner = symbol_A) %>%
    dplyr::left_join(
      dplyr::select(gene_xref, symbol, tumor_suppressor,
                    oncogene, cancer_driver), by = "symbol") %>%
    dplyr::filter(tumor_suppressor == T |
                    oncogene == T |
                    cancer_driver == T)

  tcga_coexp_db <- dplyr::bind_rows(co_expression_genes1,
                                    co_expression_genes2) %>%
    dplyr::arrange(desc(r)) %>%
    dplyr::mutate(p_value = signif(p_value, digits = 4)) %>%
    dplyr::filter(p_value < 1e-6) %>%
    dplyr::mutate(r = signif(r, digits = 3)) %>%
    dplyr::filter((r <= -0.7 & r < 0) | (r >= 0.7))


  tcgadb <- list()
  tcgadb[['co_expression']] <- tcga_coexp_db
  tcgadb[['aberration']] <- tcga_aberration_stats
  tcgadb[['recurrent_variants']] <- recurrent_tcga_variants
  tcgadb[['pfam']] <- pfam_domains
  tcgadb[['maf_codes']] <- maf_codes
  tcgadb[['site_code']] <- site_code
  tcgadb[['diagnosis_code']] <- diagnosis_code
  tcgadb[['clinical_strata_code']] <- clinical_strata_code

  saveRDS(tcgadb, file = rds_fname)

  return(tcgadb)

}

get_subcellular_annotations <- function(
  basedir = NULL,
  transcript_xref_db = NULL){


  uniprot_xref <- transcript_xref_db %>%
    dplyr::filter(property == "uniprot_acc") %>%
    dplyr::select(entrezgene, value) %>%
    dplyr::rename(uniprot_acc = value)

  go_terms <- data.frame()
  go_structure <- as.list(GO.db::GOTERM)
  for(n in names(go_structure)){
    term_df <-
      data.frame(
        'go_id' = as.character(n),
        'go_term' =
          as.character(AnnotationDbi::Term(go_structure[[n]])),
        'go_ontology' =
          AnnotationDbi::Ontology(go_structure[[n]]),
        stringsAsFactors = F)
    go_terms <- dplyr::bind_rows(go_terms, term_df)
  }

  comppi_fname <- file.path(
    basedir, "data-raw", "compppi", "comppi_proteins.txt")

  gganatogram_map_xlsx <- file.path(
    basedir, "data-raw", "compartment_mapping_jw.xlsx")

  go_gganatogram_map <-
    openxlsx::read.xlsx(gganatogram_map_xlsx,
                        sheet = 2) %>%
    janitor::clean_names() %>%
    dplyr::select(organ, go_id) %>%
    dplyr::mutate(organ = stringr::str_trim(organ)) %>%
    dplyr::left_join(gganatogram::cell_key$cell, by = "organ") %>%
    dplyr::select(-c(type,value,colour)) %>%
    dplyr::rename(ggcompartment = organ)

  comppidb <- as.data.frame(
    read.table(
      file = comppi_fname,
      sep = "\t", header = T, quote = "",
      stringsAsFactors = F,comment.char = "") %>%
      janitor::clean_names() %>%
      dplyr::filter(naming_convention == "UniProtKB/Swiss-Prot/P") %>%
      dplyr::select(
        major_loc_with_loc_score,
        protein_name,
        minor_loc,
        experimental_system_type,
        localization_source_database,
        pubmed_id) %>%
      tidyr::separate_rows(
        minor_loc, pubmed_id,
        localization_source_database,
        experimental_system_type, sep="\\|") %>%
      dplyr::rename(
        go_id = minor_loc,
        uniprot_acc = protein_name,
        annotation_source = localization_source_database) %>%
      dplyr::left_join(go_terms,by=c("go_id")) %>%
      dplyr::filter(!is.na(go_term)) %>%
      dplyr::mutate(
        annotation_type = dplyr::if_else(
          stringr::str_detect(experimental_system_type,
                              "Experimental"),
          "Experimental","Predicted/Unknown")) %>%
      dplyr::select(-experimental_system_type) %>%
      dplyr::distinct() %>%
      dplyr::left_join(uniprot_xref,by=c("uniprot_acc")) %>%
      dplyr::filter(!is.na(uniprot_acc)) %>%
      dplyr::group_by(uniprot_acc,
                      go_id,
                      go_term) %>%
      dplyr::summarise(
        annotation_source =
          paste(sort(unique(annotation_source)), collapse="|"),
        annotation_type =
          paste(sort(unique(annotation_type)), collapse="|"),
        .groups = "drop") %>%
      dplyr::mutate(
        confidence =
          stringr::str_count(annotation_source,
                             pattern = "\\|") + 1)
  )


  subcell_figure_legend <- list()
  cell_key[['cell']]$colour <- "#ffd700"
  for (i in 1:nrow(cell_key[['cell']])) {
    subcell_figure_legend[[i]] <-
      gganatogram::gganatogram(
        data=cell_key[['cell']][i,],
        outline = T,
        fillOutline='#a6bddb',
        organism="cell",
        fill="colour")  +
      ggplot2::theme_void() +
      ggplot2::ggtitle(cell_key[['cell']][i,]$organ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust=0.5, size=10)
      ) +
      ggplot2::coord_fixed()
  }

  subcelldb <- list()
  subcelldb[['comppidb']] <- comppidb
  subcelldb[['go_gganatogram_map']] <- go_gganatogram_map
  subcelldb[['gganatogram_legend']] <- subcell_figure_legend

  return(subcelldb)

}
