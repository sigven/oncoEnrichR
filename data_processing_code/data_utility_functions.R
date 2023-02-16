
#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @param chunk_size Size of PMID chunks
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed <- function(
    pmid,
    chunk_size = 100){

  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmid, ceiling(length(pmid)/chunk_size))
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

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @param raw_db_dir base directory
#' @param chunk_size Size of PMID chunks
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed2 <- function(
    pmid,
    raw_db_dir = NULL,
    chunk_size = 300){

  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(pmid, ceiling(length(pmid)/chunk_size))
  j <- 0
  all_citations <- data.frame()
  cat('Retrieving PubMed citations for PMID list, total length', length(pmid))
  cat('\n')

  eutils_medline_fname <-
    file.path(
      raw_db_dir,
      paste0("tmp_edirect_results_",
             stringi::stri_rand_strings(
               1, 20, pattern = "[A-Za-z0-9]"),
             ".txt"
      )
    )

  eutils_medline_sh <-
    file.path(
      raw_db_dir,
      paste0("tmp_edirect_cmd_",
             stringi::stri_rand_strings(
               1, 20, pattern = "[A-Za-z0-9]"),
             ".sh"
      )
    )

  write("#!/bin/sh\n\n", file = eutils_medline_sh)

  all_cmds <- c()
  while(j < length(pmid_chunks)){
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    cat('Processing chunk (edirect utils) ',j,' with ',length(pmid_chunk),'PMIDs')
    cat('\n')
    pmid_string <- paste(pmid_chunk,collapse = ",")

    eutils_cmd <- paste0(
      "/Users/sigven/edirect/esearch -db pubmed -query '",
      pmid_string,
      "' | /Users/sigven/edirect/efetch -format medline >> ",
      eutils_medline_fname)

    all_cmds <- c(all_cmds, eutils_cmd)
    j <- j + 1
  }

  write(all_cmds, file = eutils_medline_sh, append = T)
  system(paste0('chmod a+rx ',eutils_medline_sh))
  system(eutils_medline_sh, ignore.stdout = T, ignore.stderr = T)

  all_citations <- data.frame()

  if(file.exists(eutils_medline_fname)){

    con <- file(eutils_medline_fname, open = "r")
    lines <- readLines(con)
    fau <- NA
    pmid <- NA
    journal <- NA
    year <- NA

    for(l in lines){
      if(stringr::str_detect(l,"^PMID-")){

        if(!is.na(fau) & !is.na(pmid) & !is.na(journal) & !is.na(year)){
          citation <- data.frame(
            'pmid' = as.integer(pmid),
            'citation' = paste(fau,year,journal,sep=", "),
            'link' =  paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',
                             pmid,'\' target=\'_blank\'>',
                             paste(fau,year,journal,sep=", ")
                             ,'</a>'),
            stringsAsFactors = F)

          fau <- NA
          pmid <- NA
          journal <- NA
          year <- NA

          all_citations <- all_citations |>
            dplyr::bind_rows(citation)
        }

        pmid <- stringr::str_replace(l, "PMID- ","")
        cat(pmid,'\n')
      }
      if(stringr::str_detect(l,"^FAU - ") & is.na(fau)){
        fau <- paste0(
          stringr::str_split_fixed(
            stringr::str_replace(l, "FAU - ",""),",",2)[1],
          " et al.")
        cat(fau,'\n')
      }
      if(stringr::str_detect(l,"^DP(\\s+)-")){
        year <- stringr::str_split(
          stringr::str_replace(l, "DP(\\s+)- ","")," ")[[1]][1]
        cat(year,'\n')
      }
      if(stringr::str_detect(l,"^TA(\\s+)-")){
        journal <- stringr::str_replace(l, "TA(\\s+)- ","")
        cat(journal,'\n')
      }

    }
    close(con)

  }

  system('rm -f ', eutils_medline_fname)

  return(all_citations)

}

get_msigdb_signatures <- function(
    raw_db_dir = NULL,
    db_version = 'v7.5.1 (January 2022)'){

  ## get full dataset: Broad Institute's Molecular Signatures Database (v7.1)
  invisible(assertthat::assert_that(
    dir.exists(raw_db_dir),
    msg = paste0("Directory '",
                 raw_db_dir,"' does not exist")))
  msigdb_xml_fname <- file.path(
    raw_db_dir, "msigdb", "msigdb.xml")
  invisible(assertthat::assert_that(
    file.exists(msigdb_xml_fname),
    msg = paste0("File '",
                 msigdb_xml_fname,
                 "' does not exist")))

  msig_data_xml <- xml2::read_xml(msigdb_xml_fname)

  ## make data frame with signatures, one record pr. gene-signature association
  all_genesets <- msig_data_xml |> xml2::xml_find_all("//GENESET")
  category_code <- all_genesets |> xml2::xml_attr("CATEGORY_CODE")
  all_msigdb <- data.frame('category_code' = category_code, stringsAsFactors = F)
  all_msigdb$description <- all_genesets |> xml2::xml_attr("DESCRIPTION_BRIEF")
  all_msigdb$standard_name <- all_genesets |> xml2::xml_attr("STANDARD_NAME")
  all_msigdb$organism <- all_genesets |> xml2::xml_attr("ORGANISM")
  all_msigdb$pmid <- all_genesets |> xml2::xml_attr("PMID")
  all_msigdb$systematic_name <- all_genesets |> xml2::xml_attr("SYSTEMATIC_NAME")
  all_msigdb$subcategory_code <- all_genesets |> xml2::xml_attr("SUB_CATEGORY_CODE")
  all_msigdb$entrezgene <- all_genesets |> xml2::xml_attr("MEMBERS_EZID")
  all_msigdb$contributor <- all_genesets |> xml2::xml_attr("CONTRIBUTOR")
  all_msigdb$exact_source <- all_genesets |> xml2::xml_attr("EXACT_SOURCE")
  all_msigdb$external_url <- all_genesets |> xml2::xml_attr("EXTERNAL_DETAILS_URL")
  all_msigdb <- all_msigdb |>
    tidyr::separate_rows(entrezgene,sep=",") |>
    dplyr::filter(organism == "Homo sapiens") |>
    dplyr::filter(subcategory_code != 'CP:KEGG') |>
    dplyr::arrange(category_code) |>
    dplyr::mutate(pmid = dplyr::if_else(
      nchar(pmid) == 0,
      as.character(NA),
      as.character(pmid))) |>
    dplyr::mutate(subcategory_code = dplyr::if_else(
      nchar(subcategory_code) == 0,
      as.character("ALL"),
      as.character(subcategory_code))) |>
    dplyr::filter(
      category_code != "ARCHIVED" &
        category_code != "C1") |>
    dplyr::mutate(external_url = dplyr::if_else(
      nchar(external_url) == 0,
      paste0("http://software.broadinstitute.org/gsea/msigdb/cards/",standard_name),
      as.character(external_url))) |>
    dplyr::mutate(description = stringr::str_replace_all(
      description,
      "( \\[ICI [0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[GeneID=[0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[PubChem=[0-9]{1,}(;[0-9]{1,})*\\]( )?)",""))

  msigdb_category_description <- read.table(
    file = file.path(
      raw_db_dir, "msigdb",
      "msigdb_collection_description.tsv"),
    sep = "\t", header = T, stringsAsFactors = F)

  msigdb_complete <- as.data.frame(all_msigdb |>
    dplyr::left_join(msigdb_category_description,
                     by = c("category_code", "subcategory_code"),
                     multiple = "all") |>
    dplyr::mutate(db = "MSigDB", db_version = db_version) |>
    dplyr::select(db, db_version, category_code, category_description,
                  subcategory_code, subcategory_description,
                  standard_name, description, organism,
                  entrezgene) |>
    dplyr::rename(signature_description = description,
                  signature_name = standard_name)
  )


  msigdb <- list()
  msigdb[['VERSION']] <- version
  msigdb[['TERM2SOURCE']] <- all_msigdb |>
    dplyr::filter(subcategory_code != "MIR:MIR_Legacy") |>
    dplyr::filter(subcategory_code != "TFT:TFT_Legacy") |>
    dplyr::filter(subcategory_code != "CP:WIKIPATHWAYS") |>
    dplyr::filter(subcategory_code != "VAX") |>
    dplyr::select(standard_name, exact_source, external_url,
                  category_code, subcategory_code) |>
    dplyr::distinct() |>
    dplyr::mutate(external_url = stringr::str_replace(
      external_url,"\\|https://reactome.org/PathwayBrowser/","")) |>
    dplyr::mutate(external_url = dplyr::if_else(
      nchar(external_url) == 0,as.character(NA),as.character(external_url))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^GO_"),
      subcategory_code,as.character(NA))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^HP_"),
      subcategory_code,as.character(NA))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^REACTOME_"),
      "REACTOME",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^BIOCARTA_"),
      "BIOCARTA",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^PID_"),
      "PATHWAY_INTERACTION_DB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "H",
      "HALLMARK",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CGP",
      "CHEM_GEN_PERTURB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CP",
      "CANONICAL_PATHWAY",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "MIR:MIRDB",
      "MICRORNA_TARGET_MIRDB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "TFT:GTRD",
      "TF_TARGET_GTRD",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CGN",
      "CANCER_NEIGHBOURHOOD",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CM",
      "CANCER_MODULE",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C6","ONCOGENIC",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C7","IMMUNESIGDB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C8","CELLTYPE_SIGNATURES",as.character(db))) |>
    dplyr::select(-c(category_code,subcategory_code)) |>
    dplyr::filter(!is.na(db))


  msigdb[['COLLECTION']] <- list()

  for(c in c('H','C2','C3','C4','C5','C6','C7','C8')){
    msigdb[['COLLECTION']][[c]] <- list()
    subcategories <-
      unique(all_msigdb[all_msigdb$category_code == c,]$subcategory_code)
    for(scat in subcategories){
      subcat <- stringr::str_replace(scat,"GO:","")
      if(subcat == "CP:WIKIPATHWAYS" |
         subcat == "VAX" |
         subcat == "MIR:MIR_Legacy" |
         subcat == "TFT:TFT_Legacy"){
        next
      }

      msigdb[['COLLECTION']][[c]][[subcat]] <- list()
      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2GENE']] <- data.frame()
      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2GENE']] <-
        dplyr::filter(all_msigdb,
                      category_code == c &
                        subcategory_code == scat) |>
        dplyr::select(standard_name, entrezgene) |>
        dplyr::rename(entrez_gene = entrezgene)

      msigdb[['COLLECTION']][[c]][[subcat]][['TERM2NAME']] <-
        dplyr::filter(all_msigdb,
                      category_code == c &
                        subcategory_code == scat) |>
        dplyr::select(standard_name, description) |>
        dplyr::rename(name = description) |>
        dplyr::distinct()
    }
  }
  rm(all_msigdb)
  return(list('db' = msigdb, 'df' = msigdb_complete))
}


get_opentarget_associations <-
  function(raw_db_dir = NULL,
           min_overall_score = 0.1,
           min_num_sources = 2,
           release = "2022.11",
           direct_associations_only = F){

    opentarget_targets <- as.data.frame(
      readRDS(
        file.path(raw_db_dir,
               "opentargets",
               paste0("opentargets_target_",
               release,
               ".rds"))) |>
        dplyr::select(target_ensembl_gene_id,
                      SM_tractability_category,
                      SM_tractability_support,
                      AB_tractability_category,
                      AB_tractability_support) |>
        dplyr::rename(ensembl_gene_id = target_ensembl_gene_id)
      )

    opentarget_associations_raw <- as.data.frame(
      readRDS(
        file.path(raw_db_dir,
                  "opentargets",
                  paste0(
                    "opentargets_association_direct_HC_",
                    release,
                    ".rds"))
      )
    )

    opentarget_datatype_support <- as.data.frame(
      opentarget_associations_raw |>
        dplyr::select(disease_id, target_ensembl_gene_id, datatype_items) |>
        tidyr::separate_rows(datatype_items, sep=",") |>
        tidyr::separate(datatype_items,
                        into = c("datatype","n_evidence","score"),
                        sep="\\|") |>
        dplyr::group_by(disease_id, target_ensembl_gene_id) |>
        dplyr::summarise(datatype_support = paste(
          unique(sort(datatype)), collapse=";"),
          .groups = "drop")
      )

    opentarget_associations_raw <- opentarget_associations_raw |>
      dplyr::left_join(
        opentarget_datatype_support,
        by = c("disease_id", "target_ensembl_gene_id"),
        multiple = "all")


    opentarget_associations <- as.data.frame(
      opentarget_associations_raw |>
        dplyr::rename(disease_efo_id = disease_id,
                      ensembl_gene_id = target_ensembl_gene_id) |>
        dplyr::mutate(disease_efo_id = stringr::str_replace(
          disease_efo_id,":","_")) |>
        dplyr::filter(score >= min_overall_score) |>
        dplyr::filter(stringr::str_count(datatype_items,",") >= min_num_sources - 1) |>
        dplyr::mutate(association_key =
                        paste(disease_efo_id,
                              "T",
                              datatype_support,
                              score,
                              sep=":")) |>
        dplyr::group_by(ensembl_gene_id) |>
        dplyr::summarise(
          ot_association = paste(association_key, collapse = "&"),
          .groups = "drop")
    )

    ot_associations <- opentarget_targets |>
      dplyr::left_join(opentarget_associations,
                       by = "ensembl_gene_id",
                       multiple = "all")


    return(ot_associations)

  }


get_cancer_hallmarks <- function(raw_db_dir = NULL,
                                 gene_info = NULL,
                                 opentargets_version = "2022.11",
                                 update = T){

  rds_fname <-
    file.path(raw_db_dir, "opentargets",
              paste0("opentargets_hallmarkdata_",
              opentargets_version, ".rds"))

  if(update == F & file.exists(rds_fname)){
    hallmark_data <- readRDS(file = rds_fname)
    return(hallmark_data)
  }

  opentargets_target_rds <- file.path(
    raw_db_dir, "opentargets",
    paste0("opentargets_target_", opentargets_version, ".rds")
  )

  hallmark_data_long <- as.data.frame(
    readRDS(file = opentargets_target_rds) |>
      dplyr::select(target_ensembl_gene_id, hgnc_id, cancer_hallmark) |>
      dplyr::rename(ensembl_gene_id = target_ensembl_gene_id) |>
      dplyr::filter(!is.na(cancer_hallmark)) |>
      tidyr::separate_rows(cancer_hallmark, sep="@@@") |>
      tidyr::separate(cancer_hallmark,
                      into = c('hallmark','promote','suppress','literature_support'),
                      remove = F, sep = '\\|') |>
      dplyr::rename(promotes = promote, suppresses = suppress) |>
      dplyr::mutate(hallmark = stringr::str_to_title(hallmark)) |>
      tidyr::separate_rows(literature_support, sep="&") |>
      dplyr::mutate(promotes = as.logical(stringr::str_replace(promotes,"PROMOTE:",""))) |>
      dplyr::mutate(suppresses = as.logical(stringr::str_replace(suppresses,"SUPPRESS:",""))) |>
      tidyr::separate(literature_support, into = c('pmid','pmid_summary'), sep="%%%") |>
      dplyr::filter(pmid != "37937225") |>
      dplyr::mutate(pmid_summary = R.utils::capitalize(pmid_summary))
  )

  pmid_data <- get_citations_pubmed(
    unique(hallmark_data_long$pmid),
    chunk_size = 50) |>
    dplyr::mutate(pmid = as.character(pmid))


  hallmark_data_long <- hallmark_data_long |>
    dplyr::left_join(pmid_data, by = "pmid",  multiple = "all")

  hallmark_data_short <- as.data.frame(
    hallmark_data_long |>
      dplyr::select(hgnc_id, ensembl_gene_id, hallmark,
                    promotes, suppresses, pmid_summary,
                    link) |>
      dplyr::mutate(literature_support = paste0(pmid_summary," (",link,")")) |>
      dplyr::group_by(hgnc_id, ensembl_gene_id, hallmark,
                      promotes, suppresses) |>
      dplyr::summarise(literature_support = paste(
        literature_support, collapse=". "),
        .groups = "drop") |>
      dplyr::left_join(dplyr::select(gene_info, hgnc_id,
                                     symbol, entrezgene),
                       by = "hgnc_id",  multiple = "all") |>
      dplyr::mutate(symbol =
                      paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
                             entrezgene,"' target='_blank'>",symbol,"</a>")) |>
      dplyr::select(-c("hgnc_id")) |>
      dplyr::select(symbol, entrezgene, ensembl_gene_id,
                    hallmark, promotes, suppresses,
                    literature_support)
  )

  hallmark_data <- list('long' = hallmark_data_long,
                        'short' = hallmark_data_short)
  saveRDS(hallmark_data, file = rds_fname)
  return(hallmark_data)

}

get_tf_target_interactions <- function(
  raw_db_dir = NULL,
  update = F){


  rds_fname = file.path(
    raw_db_dir,
    "omnipathdb",
    "omnipath_tf_interactions.rds")

  if(update == F & file.exists(rds_fname)){
    tf_target_interactions <- readRDS(file = rds_fname)
    return(tf_target_interactions)
  }


  tf_target_interactions_dorothea_all <-
    OmnipathR::import_transcriptional_interactions(
      organism = "9606",
      dorothea_levels = c("A","B","C","D"),
      references_by_resource = F) |>
    dplyr::filter(!is.na(dorothea_level)) |>
    dplyr::group_by(source_genesymbol, target_genesymbol) |>
    dplyr::summarise(
      references = paste(unique(references), collapse=";"),
      interaction_sources = paste(unique(sources), collapse="; "),
      dorothea_level = paste(unique(dorothea_level), collapse=";"),
      .groups = "drop") |>
    dplyr::filter(!stringr::str_detect(source_genesymbol,"_")) |>
    dplyr::mutate(interaction_sources = stringr::str_squish(
      stringr::str_replace_all(
        interaction_sources, ";","; "))) |>

    ## Not available for non-academic use (entries provided solely with TRED)
    dplyr::filter(!stringr::str_detect(
      interaction_sources,"^(RegNetwork_DoRothEA; TRED_DoRothEA)$")) |>
    dplyr::mutate(interaction_source = stringr::str_replace_all(
      interaction_sources, "; TRED_DoRothEA$", ""
    )) |>
    dplyr::mutate(interaction_source = stringr::str_replace_all(
      interaction_sources, "; TRED_DoRothEA; ", "; "
    ))

  tf_target_pmids <- tf_target_interactions_dorothea_all |>
    dplyr::filter(!is.na(references) & nchar(references) > 0) |>
    dplyr::select(source_genesymbol, target_genesymbol, references) |>
    tidyr::separate_rows(references, sep=";") |>
    dplyr::filter(nchar(references) > 0) |>
    dplyr::distinct()

  tf_target_citations <- get_citations_pubmed(
    pmid = unique(tf_target_pmids$references),
    chunk_size = 100)

  tf_target_pmids <- as.data.frame(
    tf_target_pmids |>
      dplyr::mutate(references = as.integer(references)) |>
      dplyr::left_join(
        tf_target_citations, by = c("references" = "pmid"),
        multiple = "all") |>
      dplyr::group_by(source_genesymbol, target_genesymbol) |>
      dplyr::summarise(
        tf_target_literature_support = paste(
          link, collapse= ","),
        tf_target_literature = paste(
          citation, collapse="|"
        ),
        .groups = "drop"
      )
  )

  tf_target_interactions_dorothea_all <- tf_target_interactions_dorothea_all |>
    dplyr::left_join(tf_target_pmids,
                     by = c("source_genesymbol",
                            "target_genesymbol"),
                     multiple = "all") |>
    dplyr::rename(regulator = source_genesymbol,
                  target = target_genesymbol) |>
    dplyr::select(-references) |>
    dplyr::mutate(confidence_level = dplyr::case_when(
      stringr::str_detect(dorothea_level,"A|A;B|A;B;D|A;D") ~ "A",
      stringr::str_detect(dorothea_level,"B|B;D|B;C") ~ "B",
      stringr::str_detect(dorothea_level,"C;D|C") ~ "C",
      stringr::str_detect(dorothea_level,"D") ~ "D",
      TRUE ~ as.character(dorothea_level),
    )) |>
    dplyr::select(-dorothea_level) |>
    dplyr::mutate(mode_of_regulation = NA)


  tf_target_interactions_dorothea_pancancer <-
    dorothea::dorothea_hs_pancancer |>
    dplyr::mutate(interaction_sources = "DoRothEA_hs_pancancer") |>
    dplyr::rename(regulator = tf,
                  confidence_level = confidence,
                  mode_of_regulation = mor) |>
    dplyr::mutate(mode_of_regulation = dplyr::if_else(
      mode_of_regulation == 1,"Stimulation",
      "Repression")) |>
    dplyr::mutate(tf_target_literature_support = NA,
                  tf_target_literature = NA) |>
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

get_cancer_drugs <- function(raw_db_dir = NULL){

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


  approved_inhibitors <- as.data.frame(
    pharmOncoX::get_drugs(
      cache_dir = raw_db_dir,
      drug_is_targeted = T,
      inhibitor_only = T,
      drug_is_approved = T)$records |>
      dplyr::filter(!is.na(primary_site)) |>
      dplyr::filter(
        stringr::str_detect(drug_clinical_source,"FDA|DailyMed|ATC")) |>
      dplyr::select(
        disease_efo_label, target_entrezgene,
        drug_name, molecule_chembl_id) |>
      dplyr::rename(entrezgene = target_entrezgene) |>
      dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
      dplyr::mutate(
        drug_link =
          paste0("<a href = 'https://platform.opentargets.org/drug/",
                 molecule_chembl_id,
                 "' target='_blank'>",drug_name,"</a>")) |>
      dplyr::group_by(entrezgene, drug_link) |>
      dplyr::summarise(indications = paste(
        sort(unique(disease_efo_label)), collapse="; "
      ), .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(indications = stringr::str_replace(
        indications,
        "breast cancer; breast carcinoma",
        "breast cancer"
      )) |>
      dplyr::mutate(drug_indications = paste0(
        drug_link, " (",indications,")")) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(approved_drugs = paste(
        drug_indications, collapse=", "),
        .groups = "drop")
  )

  cancer_drugs[['late_phase']] <-
    pharmOncoX::get_drugs(
      cache_dir = raw_db_dir,
      drug_is_targeted = T,
      inhibitor_only = T,
      drug_minimum_phase_any_indication = 3)$records |>
    dplyr::filter(!is.na(molecule_chembl_id)) |>
    dplyr::select(target_entrezgene,
                  drug_name,
                  primary_site,
                  molecule_chembl_id,
                  drug_max_ct_phase) |>
    dplyr::distinct() |>
    dplyr::rename(entrezgene = target_entrezgene) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::distinct()

  cancer_drugs[['early_phase']] <-
    pharmOncoX::get_drugs(
      cache_dir = raw_db_dir,
      drug_is_targeted = T,
      inhibitor_only = T,
      drug_minimum_phase_any_indication = 0)$records |>
    dplyr::filter(!is.na(molecule_chembl_id)) |>
    dplyr::anti_join(
      cancer_drugs[['late_phase']], by = "molecule_chembl_id") |>
    dplyr::select(target_entrezgene,
                  drug_name,
                  primary_site,
                  molecule_chembl_id,
                  drug_max_ct_phase) |>
    dplyr::rename(entrezgene = target_entrezgene) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::distinct()

  for(p in c('early_phase','late_phase')){

    cancer_drugs[[p]] <- cancer_drugs[[p]] |>
      dplyr::mutate(to = paste0("s", molecule_chembl_id)) |>
      dplyr::mutate(from = paste0("s", entrezgene)) |>
      dplyr::distinct()

    drugs_per_gene[[p]] <- as.data.frame(
      cancer_drugs[[p]] |>
        dplyr::select(entrezgene,
                      drug_name,
                      drug_max_ct_phase,
                      molecule_chembl_id) |>
        dplyr::distinct() |>
        dplyr::mutate(
          drug_link =
            paste0("<a href = 'https://platform.opentargets.org/drug/",
                   molecule_chembl_id,"' target='_blank'>",drug_name,"</a>")) |>
        dplyr::arrange(entrezgene, desc(drug_max_ct_phase)) |>
        dplyr::group_by(entrezgene) |>
        dplyr::summarise(
          targeted_cancer_drugs =
            paste(unique(drug_link),
                  collapse = ", "), .groups = "drop")
    )
    if(p == 'early_phase'){
      drugs_per_gene[[p]] <- drugs_per_gene[[p]] |>
        dplyr::rename(
          targeted_cancer_drugs_ep = targeted_cancer_drugs)
    }else{
      drugs_per_gene[[p]] <- drugs_per_gene[[p]] |>
        dplyr::rename(
          targeted_cancer_drugs_lp = targeted_cancer_drugs)
    }

    target_drug_edges[[p]] <- cancer_drugs[[p]] |>
      dplyr::select(drug_name,
                    from,
                    to,
                    #target_symbol,
                    entrezgene,
                    molecule_chembl_id) |>
      dplyr::distinct() |>
      dplyr::mutate(width = 1)

    target_drug_nodes[[p]] <- as.data.frame(
      cancer_drugs[[p]] |>
        dplyr::select(to,
                      drug_name,
                      primary_site) |>
        dplyr::distinct() |>
        dplyr::rename(id = to,
                      label = drug_name) |>
        dplyr::group_by(id, label) |>
        dplyr::summarise(
          primary_site = paste(sort(unique(primary_site)),
                               collapse = ", "),
          .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::mutate(
          primary_site =
            dplyr::if_else(is.na(primary_site) |
                             nchar(primary_site) == 0,
                           "Unspecific indication",
                           as.character(primary_site))) |>
        dplyr::mutate(
          title = paste0(label, " (", primary_site,")")) |>
        dplyr::mutate(
          shadow = F, shape = "diamond",
          color.background = "orange",
          font.color = "black")
    )
    if(p == 'early_phase'){
      target_drug_nodes[[p]] <- target_drug_nodes[[p]] |>
        dplyr::mutate(shadow = F,
                      shape = "diamond",
                      color.background = "purple",
                      font.color = "black")
    }

  }
  cancerdrugdb <- list()
  cancerdrugdb[['approved_per_target']] <- approved_inhibitors
  cancerdrugdb[['drug_per_target']] <- list()
  cancerdrugdb[['drug_per_target']][['early_phase']] <-
    drugs_per_gene[['early_phase']]
  cancerdrugdb[['drug_per_target']][['late_phase']] <-
    drugs_per_gene[['late_phase']]

  cancerdrugdb[['network']] <- list()
  cancerdrugdb[['network']][['edges']] <-
    dplyr::bind_rows(target_drug_edges[['early_phase']],
                     target_drug_edges[['late_phase']]) |>
    dplyr::select(from, to, width) |>
    dplyr::distinct()

  cancerdrugdb[['network']][['nodes']] <-
    dplyr::bind_rows(target_drug_nodes[['early_phase']],
                     target_drug_nodes[['late_phase']]) |>
    dplyr::distinct()

  return(cancerdrugdb)


}

get_pathway_annotations <- function(
  raw_db_dir = NULL,
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
    raw_db_dir, "wikipathways",
    paste0("wikipathways-", wikipathways_version,
           "-gmt-Homo_sapiens.gmt")
  )

  netpath_mapping_fname <- file.path(
    raw_db_dir, "netpath",
    "id_mapping.tsv")

  keggdb_fname <- file.path(
    raw_db_dir,
    "kegg", "kegg.pathway.gene.tsv.gz"
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
        destpath = file.path(raw_db_dir,"wikipathways"),
        format = "gmt")
    }

    wikipathways <- clusterProfiler::read.gmt(
      wikipathways_gmt)
    wp2gene <- wikipathways |>
      tidyr::separate(term, c("name","version","wpid","org"), "%")
    wikipathwaydb <- list()
    wikipathwaydb[['VERSION']] <- wikipathways_version
    wikipathwaydb[['TERM2GENE']] <- wp2gene |>
      dplyr::select(wpid, gene) |>
      dplyr::rename(standard_name = wpid, entrez_gene = gene)
    wikipathwaydb[['TERM2NAME']] <- wp2gene |>
      dplyr::select(wpid, name) |>
      dplyr::rename(standard_name = wpid) |>
      dplyr::distinct()

    pathwaydb[['wikipathways']] <- wikipathwaydb
  }


  if(netpath == T & !is.null(omnipathdb) & !is.null(gene_info)){
      netpath_idmapping <-
      read.table(netpath_mapping_fname,
                 stringsAsFactors = F, header = F, sep = "\t",
                 col.names = c("name", "standard_name")) |>
      dplyr::mutate(standard_name = paste0("NetPath_",standard_name))

    netpath_pathway_data <- omnipathdb |>
      dplyr::filter(source == "NetPath") |>
      dplyr::select(genesymbol, value) |>
      dplyr::rename(name = value, symbol = genesymbol) |>
      dplyr::left_join(netpath_idmapping, by = "name",  multiple = "all") |>
      dplyr::mutate(name = paste0(name," signaling pathway")) |>
      dplyr::arrange(name, standard_name, symbol) |>
      dplyr::left_join(dplyr::select(
        gene_info, entrezgene, symbol),
        by = "symbol",  multiple = "all") |>
      dplyr::rename(entrez_gene = entrezgene) |>
      dplyr::mutate(entrez_gene = as.character(entrez_gene)) |>
      dplyr::select(standard_name, name, entrez_gene)
    netpathdb <- list()
    netpathdb[['VERSION']] <- netpath_version
    netpathdb[['TERM2GENE']] <- netpath_pathway_data |>
      dplyr::select(standard_name, entrez_gene)
    netpathdb[['TERM2NAME']] <- netpath_pathway_data |>
      dplyr::select(standard_name, name) |>
      dplyr::distinct()

    pathwaydb[['netpath']] <- netpathdb
  }


  if(kegg == T){
    kegg_pathway_genes <- as.data.frame(
      readr::read_tsv(file = keggdb_fname,
                 col_names = T, quote = "",
                 show_col_types = F))
    keggdb <- list()
    keggdb[['VERSION']] <- kegg_version
    keggdb[['TERM2NAME']] <- kegg_pathway_genes |>
      dplyr::select(name,pathway_id) |>
      dplyr::rename(standard_name = pathway_id) |>
      dplyr::select(standard_name, name) |>
      dplyr::distinct()

    keggdb[['TERM2GENE']] <- kegg_pathway_genes |>
      dplyr::select(pathway_id, gene_id) |>
      dplyr::rename(standard_name = pathway_id,
                    entrez_gene = gene_id) |>
      dplyr::select(standard_name,entrez_gene) |>
      dplyr::mutate(entrez_gene = as.character(
        entrez_gene)) |>
      dplyr::distinct()

    pathwaydb[['kegg']] <- keggdb
  }

  if(msigdb == T){

    msigdb <-
      get_msigdb_signatures(raw_db_dir = raw_db_dir,
                            db_version = msigdb_version)

    pathwaydb[['msigdb']] <- msigdb$db

  }

  return(pathwaydb)

}

get_protein_complexes <- function(
    raw_db_dir = NULL,
    update = F){

    rds_fname <- file.path(
      raw_db_dir,
      "protein_complexes", "omnipath_complexdb.rds")

    humap2_complexes_tsv <-file.path(
      raw_db_dir,
      "protein_complexes", "humap2_complexes.tsv")

    corum_complexes_tsv <- file.path(
      raw_db_dir,
      "protein_complexes", "coreComplexes.txt")

    humap_url <-
      "http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt"


    if(update == F & file.exists(rds_fname)){
      proteincomplexdb <- readRDS(rds_fname)
      return(proteincomplexdb)
    }

    ## Protein complex data from OmniPath - multiple underlying databases
    ## - CORUM, ComplexPortal, Compleat etc

    protein_complexes <- list()

    protein_complexes[['OmniPath']] <- OmnipathR::import_omnipath_complexes() |>
      dplyr::filter(!is.na(references) & !is.na(name)) |>
      dplyr::rename(pmid = references,
                    uniprot_acc = components,
                    complex_name = name) |>
      dplyr::select(complex_name, uniprot_acc,
                    pmid, sources) |>
      dplyr::mutate(complex_name = Hmisc::capitalize(complex_name)) |>
      dplyr::distinct() |>
      dplyr::filter(stringr::str_detect(uniprot_acc,"_")) |>
      tidyr::separate_rows(uniprot_acc, sep = "_") |>
      tidyr::separate_rows(sources, sep=";") |>
      tidyr::separate_rows(pmid, sep=";") |>
      dplyr::filter(nchar(pmid) > 3) |>
      dplyr::distinct() |>
      dplyr::mutate(complex_name = dplyr::case_when(
        complex_name == "26S proteasome" ~ "26S Proteasome",
        complex_name == "Abi1-Wasl complex" ~ "ABI1-WASL complex",
        complex_name == "Brm-associated complex" ~ "BRM-associated complex",
        complex_name == "Ski Complex" ~ "SKI complex",
        complex_name == "SF3b complex" ~ "SF3B complex",
        complex_name == "Sin3 complex" ~ "SIN3 complex",
        complex_name == "Sos1-Abi1-Eps8 complex" ~ "SOS1-ABI1-EPS8 complex",
        complex_name == "DCS complex (Ptbp1, Ptbp2, Hnrph1, Hnrpf)" ~
          "DCS complex (PTBP1, PTBP2, HNRPH1, HNRPF)",
        complex_name == "Tiam1-Efnb1-Epha2 complex" ~ "TIAM1-EFNB1-EPHA2 complex",
        complex_name == "Tip5-Dnmt-Hdac1 complex" ~ "TIP5-DNMT-HDAC1 complex",
        TRUE ~ as.character(complex_name)
      )) |>
      dplyr::group_by(complex_name) |>
      dplyr::summarise(
        uniprot_acc = paste(sort(unique(uniprot_acc)),
                            collapse = "_"),
        pmid = paste(sort(unique(pmid)),
                     collapse=";"),
        sources = paste(sort(unique(sources)),
                        collapse=";")) |>
      dplyr::distinct() |>
      dplyr::mutate(num_sources = stringr::str_count(sources,";") + 1) |>
      dplyr::mutate(complex_name = stringr::str_replace_all(
        complex_name,"\\\\u2013","-"
      )) |>
      dplyr::mutate(sources = stringr::str_replace_all(
        sources, ";", "; ")) |>
      dplyr::mutate(complex_id = dplyr::row_number())

     pmid2complex <- protein_complexes[['OmniPath']] |>
       dplyr::select(complex_id,
                     pmid) |>
       tidyr::separate_rows(pmid, sep=";") |>
       dplyr::filter(stringr::str_detect(pmid,"^[0-9]{4,}$"))

     protein_complex_citations <- get_citations_pubmed(
       pmid = unique(pmid2complex$pmid))

     pmid2complex <- pmid2complex |>
       dplyr::mutate(pmid = as.integer(pmid)) |>
       dplyr::left_join(protein_complex_citations, by = "pmid", multiple = "all") |>
       dplyr::group_by(complex_id) |>
       dplyr::summarise(
         complex_literature_support = paste(
           link, collapse= ", "),
         complex_literature = paste(
           citation, collapse="|"
         )
       )

     protein_complexes[['OmniPath']] <-
       protein_complexes[['OmniPath']] |>
       dplyr::left_join(pmid2complex,
                       by = "complex_id",  multiple = "all")

    ## CORUM annotations (complex comments, disease comments, purification methods)
    protein_complexes[['CORUM']] <- as.data.frame(
      readr::read_tsv(
        corum_complexes_tsv,
        show_col_types = F) |>
        janitor::clean_names() |>
        dplyr::rename(pmid = pub_med_id,
                      uniprot_acc = subunits_uni_prot_i_ds) |>
        dplyr::mutate(
          complex_name = stringr::str_replace_all(
            complex_name, "\\\\u2013", "-"
          )
        ) |>
        dplyr::mutate(
          complex_name = dplyr::if_else(
            complex_name == "VHL-ElonginB-ElonginC complex",
            "VHL-TCEB1-TCEB2 complex",
            as.character(complex_name)
          )
        ) |>
        dplyr::mutate(
          complex_name = dplyr::if_else(
            complex_name == "TGF-beta receptor II-TGF-beta3 complex",
            "TGF-betaR2-TGF-beta3 complex",
            as.character(complex_name)
          )
        ) |>
        dplyr::filter(organism == "Human") |>
        tidyr::separate_rows(uniprot_acc, sep = ";") |>
        tidyr::separate_rows(pmid, sep = ";") |>
        dplyr::filter(nchar(pmid) > 2) |>
        dplyr::group_by(complex_name) |>
        dplyr::summarise(
          purification_method =
            paste(unique(protein_complex_purification_method),
                  collapse = "; "),
          complex_comment =
            paste(unique(complex_comment),
                  collapse="; "),
          disease_comment =
            paste(unique(disease_comment),
                  collapse="; "),
          # pmid =
          #   paste(unique(pmid),
          #         collapse=";"),
          # uniprot_acc = paste(sort(unique(uniprot_acc)),
          #                     collapse = "_"),
          .groups = "drop"
        ) |>
        dplyr::distinct()
    )


    # missing_corum_entries <- protein_complexes[['CORUM']] |>
    #   dplyr::anti_join(protein_complexes[['OmniPath']], by = "complex_name") |>
    #   dplyr::anti_join(protein_complexes[['OmniPath']], by = "uniprot_acc") |>
    #   dplyr::filter(stringr::str_detect(
    #     uniprot_acc, "_"
    #   ))

    protein_complexes[['OmniPath']] <-
      protein_complexes[['OmniPath']] |>
      dplyr::left_join(
        protein_complexes[['CORUM']], by = "complex_name",  multiple = "all") |>
      #dplyr::select(-pmid) |>
      dplyr::mutate(confidence = NA) |>
      dplyr::mutate(complex_id = as.character(complex_id))

    if(!file.exists(humap2_complexes_tsv)){
      download.file(url = humap_url,
                    destfile = humap2_complexes_tsv)
    }

    complex_humap <- readr::read_csv(
      file = humap2_complexes_tsv,
      show_col_types = F) |>
      janitor::clean_names() |>
      dplyr::rename(complex_id = hu_map2_id,
                    uniprot_acc = uniprot_ac_cs) |>
      dplyr::select(-genenames) |>
      dplyr::mutate(uniprot_acc = stringr::str_replace_all(
        uniprot_acc, " ","_"
      )) |>
      dplyr::mutate(complex_name = complex_id) |>
      dplyr::mutate(sources = "hu.MAP 2.0") |>
      dplyr::anti_join(
        protein_complexes[['OmniPath']],
        by = "uniprot_acc")

    complex_all <- protein_complexes[['OmniPath']] |>
      dplyr::bind_rows(complex_humap)

    proteincomplexdb <- list()
    proteincomplexdb[["db"]] <- complex_all |>
      dplyr::select(-uniprot_acc) |>
      dplyr::distinct()

    proteincomplexdb[["up_xref"]] <- complex_all |>
      dplyr::select(complex_id, uniprot_acc) |>
      tidyr::separate_rows(uniprot_acc, sep="_") |>
      dplyr::distinct()


    saveRDS(proteincomplexdb,
            file = rds_fname)
    return(proteincomplexdb)
}

clean_project_score_cell_lines <- function(
  mat){

  ## Ignore cell lines that fail QC
  mat <- mat[
    , mat[3,] != "False"
  ]

  rownames(mat) <-
    mat[,"model_name"]
  mat <- mat[,-1]
  colnames(mat) <-
    mat["model_id",]
  mat <- mat[-1,]
  mat <- mat[-2,]
  mat <- mat[-2,]

  duplicate_cell_lines <-
    plyr::count(colnames(mat)) |>
    dplyr::filter(freq == 2)

  ## Mark cell lines included for analysis
  ## - if duplicates - keep Sanger cell line
  cell_line_identifiers <- colnames(mat)
  for(i in 1:length(cell_line_identifiers)){
    if(cell_line_identifiers[i] %in% duplicate_cell_lines$x){
      if(mat[1,i] == "Sanger"){
        mat[1,i] <- "Inc"
      }else{
        mat[1,i] <- "Ex"
      }
    }else{
      mat[1,i] <- "Inc"
    }
  }

  ## Exclude duplicate cell lines (Broad)
  mat <- mat[
    , mat[1,] != "Ex"
  ]
  mat <- mat[-1,]

  return(mat)

}

get_fitness_data_crispr <- function(
    raw_db_dir = NULL,
    gene_info = NULL){

    cell_lines <- read.csv2(
      file = file.path(
        raw_db_dir,
        "crispr_viability",
        "model_list.csv"),
      comment.char="",sep=",",
      header = T, stringsAsFactors = F,
      fill = T,na.strings = c("")) |>
      janitor::clean_names() |>
      dplyr::select(model_id, model_name, synonyms, model_type,
                    crispr_ko_data, tissue, cancer_type,
                    tissue_status, cancer_type_detail, cancer_type_ncit_id,
                    sample_site, gender, species, ethnicity, mutational_burden,
                    ploidy_wes, msi_status) |>
      dplyr::filter(model_type == "Cell Line" &
                      species == "Homo Sapiens") |>
      dplyr::filter(tissue != "Unknown" &
                      tissue != "Adrenal Gland" &
                      tissue != "Placenta" &
                      tissue != "Vulva") |>
      dplyr::mutate(tissue = dplyr::if_else(tissue == "Central Nervous System",
                                            "CNS/Brain",
                                            as.character(tissue))) |>
      dplyr::mutate(tissue = dplyr::if_else(tissue == "Large Intestine",
                                            "Colon/Rectum",
                                            as.character(tissue))) |>
      dplyr::mutate(tissue = dplyr::if_else(tissue == "Small Intestine",
                                            "Stomach",
                                            as.character(tissue)))

    gene_identifiers <-
      read.csv(file = file.path(
        raw_db_dir,
        "crispr_viability",
        "gene_identifiers.csv"),
        stringsAsFactors = F) |>
      dplyr::rename(gene_id_project_score = gene_id,
                    entrezgene = entrez_id) |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::left_join(
        dplyr::select(gene_info, symbol,
                      entrezgene), by = "entrezgene",  multiple = "all") |>
      dplyr::select(gene_id_project_score, entrezgene, symbol)

    bayes_factors <- as.matrix(
      readr::read_tsv(
        file = file.path(
          raw_db_dir,
          "crispr_viability",
          "scaledBayesianFactors.tsv.gz"),
        show_col_types = F)) |>
      clean_project_score_cell_lines()

    bayes_factors <-
      as.data.frame(setNames(reshape2::melt(bayes_factors, na.rm = T),
                             c('symbol', 'model_id', 'scaled_BF'))) |>
      dplyr::mutate(model_id = as.character(model_id),
                    symbol = as.character(symbol),
                    scaled_BF = as.numeric(scaled_BF) * -1) |>
      dplyr::filter(!is.na(scaled_BF)) |>
      dplyr::filter(scaled_BF < 0)

    ## Fitness scores - new
    cell_fitness_scores <- as.matrix(readr::read_tsv(
      file = file.path(
        raw_db_dir,
        "crispr_viability",
        "binaryFitnessScores.tsv.gz"),
        show_col_types = F)) |>
      clean_project_score_cell_lines()


    projectscoredb <- list()
    projectscoredb[['fitness_scores']] <-
      as.data.frame(setNames(reshape2::melt(cell_fitness_scores, na.rm = T),
                             c('symbol', 'model_id', 'loss_of_fitness'))) |>
      dplyr::mutate(model_id = as.character(model_id),
                    symbol = as.character(symbol),
                    loss_of_fitness = as.integer(loss_of_fitness)) |>
      dplyr::filter(loss_of_fitness == 1) |>
      dplyr::left_join(
        dplyr::select(cell_lines, model_id, model_name, tissue,
                      cancer_type, sample_site, tissue_status),
        by = c("model_id"),  multiple = "all") |>
      dplyr::filter(!is.na(model_name)) |>
      dplyr::left_join(gene_identifiers,by=c("symbol"),  multiple = "all") |>
      dplyr::filter(!is.na(gene_id_project_score)) |>
      dplyr::left_join(
        bayes_factors, by =
          c("symbol", "model_id"),  multiple = "all") |>
      dplyr::arrange(symbol, scaled_BF)

    ## Target priority scores
    projectscoredb[['target_priority_scores']] <-
      read.csv(file = file.path(
        raw_db_dir,
        "crispr_viability", "depmap-priority-scores.csv"),
        header = T, quote="", stringsAsFactors = F) |>
      purrr::set_names(
        c("tractability_bucket","gene_id","priority_score",
          "symbol","analysis_id","tumor_type")) |>
      dplyr::mutate(priority_score = as.numeric(
        stringr::str_replace_all(priority_score,"\\\"",""))) |>
      dplyr::mutate(symbol = as.character(
        stringr::str_replace_all(symbol,"\\\"",""))) |>
      dplyr::mutate(gene_id = as.character(
        stringr::str_replace_all(gene_id,"\\\"",""))) |>
      dplyr::mutate(tumor_type = as.character(
        stringr::str_replace_all(tumor_type,"\\\"",""))) |>
      dplyr::select(symbol, gene_id, tumor_type,
                    priority_score) |>
      dplyr::mutate(priority_score =
                      round(priority_score,
                            digits = 3)) |>
      dplyr::mutate(tumor_type = dplyr::case_when(
        tumor_type == "Pan-Cancer" ~ "Pancancer",
        tumor_type == "Haematopoietic and Lyphoid" ~
          "Haematopoietic and Lymphoid",
        TRUE ~ as.character(tumor_type)
      )) |>
      dplyr::distinct()

    ## order symbols by pancancer priority rank
    priority_order <-
      projectscoredb$target_priority_scores |>
      dplyr::filter(tumor_type == "Pancancer") |>
      dplyr::arrange(desc(priority_score))

    projectscoredb[['target_priority_scores']] <-
      projectscoredb[['target_priority_scores']] |>
      dplyr::mutate(symbol = factor(symbol,
                                    levels = priority_order$symbol))

    return(projectscoredb)

}


get_survival_associations <- function(
  gene_info = NULL,
  raw_db_dir = NULL){

  for(feature_type in
      c('CNAs','Mutations','Gene expression','Methylation',
        'miRNA expression','Protein expression')){
    destfile_fname <-
      file.path(raw_db_dir, "km_survival_cshl",
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

    fname <- file.path(raw_db_dir,
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
                             c('symbol', 'tcga_cohort', 'z_score'))) |>
      #dplyr::mutate(feature_type = feattype_brief) |>
      dplyr::mutate(z_score = as.numeric(stringr::str_trim(z_score))) |>
      dplyr::mutate(symbol = stringr::str_replace(
        symbol,"^'","")) |>
      dplyr::filter(!stringr::str_detect(tcga_cohort,"^Stouffer")) |>
      dplyr::left_join(dplyr::select(gene_info,
                                     symbol, gene_biotype),
                       by = "symbol",  multiple = "all") |>
      ## limit associations to protein-coding genes
      dplyr::filter(gene_biotype == "protein-coding") |>
      dplyr::select(-gene_biotype)

    pancancer_trend <- as.data.frame(
      project_survival_df |>
        dplyr::group_by(symbol) |>
        dplyr::summarise(tot_z_score = sum(z_score)) |>
        dplyr::arrange(desc(tot_z_score)))

    project_survival[[feattype_brief]] <- project_survival_df |>
      dplyr::mutate(
        symbol = factor(symbol, levels = pancancer_trend$symbol))
  }
  projectsurvivaldb <- project_survival

  return(projectsurvivaldb)

}

get_unique_transcript_xrefs <- function(
  raw_db_dir = NULL,
  gene_oncox = NULL,
  update = F){

  rds_fname <- file.path(
    raw_db_dir,
    "transcript_xref",
    "transcript_xref_db.rds"
  )

  if(update == F & file.exists(rds_fname)){
    transcript_xref_db <- readRDS(file = rds_fname)
    return(transcript_xref_db)
  }


  gene_oncox$alias$records <-
    gene_oncox$alias$records |>
    dplyr::filter(n_primary_map == 1 &
                    alias != symbol) |>
    dplyr::select(alias, entrezgene) |>
    dplyr::arrange(entrezgene) |>
    dplyr::mutate(property = "alias") |>
    dplyr::rename(value = alias) |>
    dplyr::mutate(entrezgene = as.integer(
      entrezgene
    ))

  gencode_merged <-
    gene_oncox$gencode$records[['grch37']] |>
    dplyr::bind_rows(gene_oncox$gencode$records[['grch38']]) |>
    dplyr::select(
      ensembl_gene_id,
      symbol,
      entrezgene,
      refseq_mrna,
      refseq_peptide,
      uniprot_acc,
      name,
      ensembl_transcript_id,
      ensembl_protein_id,
      gene_biotype
    ) |>
    dplyr::filter(
      !stringr::str_detect(entrezgene,"&")) |>
    dplyr::mutate(
      ensembl_protein_id = stringr::str_replace(
        ensembl_protein_id, "\\.[0-9]{1,}$",""
      ))

  transcript_xref_db <- data.frame()

  for(xref in c('ensembl_transcript_id',
                'uniprot_acc',
                'ensembl_gene_id',
                'ensembl_protein_id',
                'refseq_mrna',
                'symbol',
                'refseq_peptide')){

    tmp <- gencode_merged |>
      dplyr::select(entrezgene, !!rlang::sym(xref)) |>
      tidyr::separate_rows(!!rlang::sym(xref), sep ="&") |>
      dplyr::distinct() |>
      dplyr::filter(!is.na(!!rlang::sym(xref))) |>
      dplyr::group_by(!!rlang::sym(xref)) |>
      dplyr::summarise(n = dplyr::n(),
                       entrezgene = paste(unique(entrezgene),collapse=",")) |>
      dplyr::filter(n == 1) |>
      dplyr::mutate(property = xref,
                    value = !!rlang::sym(xref)) |>
      dplyr::select(entrezgene, property, value)

    transcript_xref_db <- dplyr::bind_rows(
      transcript_xref_db, tmp
    )

  }

  gene_oncox$alias$records$entrezgene <-
    as.character(gene_oncox$alias$records$entrezgene)

  transcript_xref_db <- dplyr::bind_rows(
    transcript_xref_db, gene_oncox$alias$records) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene))

  saveRDS(transcript_xref_db, file = rds_fname)
  return(transcript_xref_db)

}

get_omnipath_gene_annotations <- function(
  raw_db_dir = NULL,
  gene_info = NULL,
  update = F){

  rds_fname <- file.path(
    raw_db_dir,
    "omnipathdb",
    "omnipathdb.rds")

  if(update == F & file.exists(rds_fname)){
    omnipathdb <- readRDS(file = rds_fname)
    return(omnipathdb)
  }

  protein_coding_genes <- gene_info |>
    dplyr::filter(gene_biotype == "protein-coding") |>
    dplyr::select(symbol) |>
    dplyr::filter(!stringr::str_detect(symbol,"-")) |>
    dplyr::filter(symbol != "XYLB") |>
    dplyr::distinct()

  i <- 1
  omnipathdb <- data.frame()
  while(i <= nrow(protein_coding_genes) - 200){
    genes <-
      protein_coding_genes$symbol[i:min(nrow(protein_coding_genes),i + 199)]
    annotations <-
      OmnipathR::import_omnipath_annotations(proteins = genes) |>
      dplyr::filter(
        !stringr::str_detect(
          source, "^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$")) |>
      dplyr::rename(uniprot_acc = uniprot) |>
      dplyr::mutate(record_id = as.character(record_id))
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
      select_genes = genes) |>
    dplyr::filter(!stringr::str_detect(
      source, "^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$")) |>
    dplyr::rename(uniprot_acc = uniprot) |>
    dplyr::mutate(record_id = as.character(record_id))

  if(nrow(annotations) > 0){
    omnipathdb <-
      dplyr::bind_rows(omnipathdb, annotations)
  }

  saveRDS(omnipathdb, file = rds_fname)
  return(omnipathdb)
}

get_gene_go_terms <- function(
  raw_db_dir = NULL){


  goa_human_fname <-
    file.path(
      raw_db_dir,
      "gene_ontology",
      "goa_human.gaf.gz")

  go_terms <- data.frame()
  go_structure <- as.list(GO.db::GOTERM)
  for(n in names(go_structure)){
    term_df <-
      data.frame(
        'go_id' = as.character(n),
        'go_term' =
          as.character(AnnotationDbi::Term(go_structure[[n]])),
        stringsAsFactors = F)
    go_terms <- dplyr::bind_rows(go_terms, term_df)
  }

  go_annotations_qc <-
    read.table(gzfile(goa_human_fname),
               sep = "\t", comment.char = "!", header = F,
               stringsAsFactors = F, quote = "") |>
    dplyr::select(V3, V5, V7, V9, V14) |>
    dplyr::rename(symbol = V3, go_id = V5,
                  go_evidence_code = V7,
                  go_ontology = V9, go_annotation_date = V14) |>
    dplyr::distinct() |>
    dplyr::filter(go_evidence_code != "IEA" &
                    go_evidence_code != "ND") |>
    dplyr::filter(go_ontology != "C") |>
    dplyr::select(-go_annotation_date) |>
    dplyr::distinct() |>
    dplyr::left_join(
      go_terms, by = "go_id", multiple = "all"
    ) |>
    dplyr::mutate(
      go_term_link =
        paste0('<a href=\'http://amigo.geneontology.org/amigo/term/',
               .data$go_id,'\' target=\'_blank\'>',
               .data$go_term, '</a>'))

  go_function_terms_prgene <- go_annotations_qc |>
    dplyr::filter(go_ontology == "F") |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(
      num_ids_function = length(unique(go_id)),
      go_term_link_mf = paste(
        go_term_link, collapse = ", "
      ),
      .groups = "drop"
    )

  go_process_terms_prgene <- go_annotations_qc |>
    dplyr::filter(go_ontology == "P") |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(
      num_ids_process = length(unique(go_id)),
      go_term_link_proc = paste(
        go_term_link, collapse = ", "
      ),
      .groups = "drop"
    )

  go_terms <- as.data.frame(
    dplyr::full_join(go_function_terms_prgene,
                     go_process_terms_prgene,
                     by = "symbol",  multiple = "all") |>
      dplyr::filter(symbol != "") |>
      dplyr::mutate(
        num_ids_process = dplyr::if_else(
          is.na(num_ids_process),
          as.integer(0),
          as.integer(num_ids_process)
        )
      ) |>
      dplyr::mutate(
        num_ids_function = dplyr::if_else(
          is.na(num_ids_function),
          as.integer(0),
          as.integer(num_ids_function)
        )
      ) |>
      dplyr::mutate(num_go_terms = num_ids_function +
                      num_ids_process) |>
      dplyr::mutate(
        go_term_link = dplyr::if_else(
          num_go_terms <= 3 &
            num_go_terms > 0,
          paste0(go_term_link_proc,
                 ", ",
                 go_term_link_mf),
          as.character("")
        )
      ) |>
      dplyr::mutate(go_term_link = stringr::str_replace(
        go_term_link, "NA, |, NA",""
      )) |>
      dplyr::select(symbol, go_term_link, num_go_terms)
  )

  return(go_terms)
}
assign_unknown_function_rank <- function(
  gene_xref = NULL){

  gene_xref <- gene_xref |>
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
          TRUE ~ as.integer(7)
        )
    ) |>

    ## make exception for family of clustered histones (not properly mapped with GO)
    dplyr::mutate(unknown_function_rank = dplyr::if_else(
      stringr::str_detect(name, "clustered histone|H3.3 histone"),
      as.integer(8),
      as.integer(unknown_function_rank)
    ))


  return(gene_xref)

}


generate_gene_xref_df <- function(
  raw_db_dir = NULL,
  gene_info = NULL,
  transcript_xref_db = NULL,
  ts_oncogene_annotations = NULL,
  opentarget_associations = NULL,
  gene_oncox = NULL,
  cancerdrugdb = NULL,
  otdb = NULL,
  go_terms_pr_gene = NULL,
  update = T){

  rds_fname <- file.path(
    raw_db_dir,
    "gene_xref",
    "gene_xref.rds"
  )

  if(update == F & file.exists(rds_fname)){
    gene_xref <- readRDS(file = rds_fname)
  }


  invisible(assertthat::assert_that(!is.null(gene_info)))
  invisible(assertthat::assert_that(!is.null(transcript_xref_db)))
  invisible(assertthat::assert_that(!is.null(ts_oncogene_annotations)))
  invisible(assertthat::assert_that(!is.null(cancerdrugdb)))
  invisible(assertthat::assert_that(!is.null(otdb)))
  invisible(assertthat::assert_that(!is.null(go_terms_pr_gene)))
  invisible(assertthat::assert_that(!is.null(opentarget_associations)))

  invisible(assertthat::assert_that(is.data.frame(gene_info)))
  invisible(assertthat::assert_that(is.data.frame(transcript_xref_db)))
  invisible(assertthat::assert_that(is.data.frame(ts_oncogene_annotations)))
  invisible(assertthat::assert_that(typeof(cancerdrugdb) == "list"))
  invisible(assertthat::assert_that(typeof(otdb) == "list"))
  invisible(assertthat::assert_that(is.data.frame(go_terms_pr_gene)))
  invisible(assertthat::assert_that(is.data.frame(opentarget_associations)))

  invisible(assertthat::assert_that(
    "drug_per_target" %in% names(cancerdrugdb)))
  invisible(assertthat::assert_that(
    "early_phase" %in% names(cancerdrugdb[['drug_per_target']])))
  invisible(assertthat::assert_that(
    "late_phase" %in% names(cancerdrugdb[['drug_per_target']])))
  invisible(assertthat::assert_that(
    "max_site_rank" %in% names(otdb)))


  assertable::assert_colnames(
    otdb[['max_site_rank']],
    c('ensembl_gene_id','cancer_max_rank'), only_colnames = F,
    quiet = T)
  assertable::assert_colnames(
    gene_info, c('entrezgene', 'name', 'symbol',
                  'hgnc_id', 'gene_biotype'), only_colnames = F,
    quiet = T)
  assertable::assert_colnames(
    transcript_xref_db, c('entrezgene', 'property',
                 'value'), quiet = T)

  assertable::assert_colnames(
    go_terms_pr_gene,
    c('symbol','num_go_terms','go_term_link'),
    quiet = T)

  assertable::assert_colnames(
    ts_oncogene_annotations,
    c('entrezgene','tumor_suppressor','oncogene', 'cancer_driver'),
    quiet = T, only_colnames = F)

  assertable::assert_colnames(
    opentarget_associations,
    c('ensembl_gene_id','SM_tractability_category','SM_tractability_support',
      'AB_tractability_category','AB_tractability_support'),
    quiet = T, only_colnames = F)


  go2entrez <- transcript_xref_db |>
    dplyr::filter(property == "symbol") |>
    dplyr::rename(symbol = value) |>
    dplyr::left_join(go_terms_pr_gene, by = "symbol", multiple = "all") |>
    dplyr::filter(!is.na(num_go_terms)) |>
    dplyr::select(-c(property, symbol)) |>
    dplyr::distinct() |>
    dplyr::group_by(entrezgene) |>
    dplyr::summarise(
      num_go_terms = max(num_go_terms),
      go_term_link = paste(
        unique(go_term_link), collapse=", "),
      .groups = "drop"
  )

  upsummary2entrez <- gene_oncox$basic$records |>
    dplyr::select(entrezgene, dbnsfp_function_description) |>
    dplyr::mutate(gene_summary_uniprot = dplyr::if_else(
      !is.na(dbnsfp_function_description),
      paste0("<b>UniProt:</b> ",
             dbnsfp_function_description),
      as.character(NA))) |>
    dplyr::select(entrezgene, gene_summary_uniprot)

  ncbisummary2entrez <- gene_oncox$basic$records |>
    dplyr::select(entrezgene, ncbi_function_summary) |>
    dplyr::mutate(gene_summary_ncbi = dplyr::if_else(
      !is.na(ncbi_function_summary),
      paste0(
        "<b>NCBI/RefSeq/OMIM:</b> ",
        ncbi_function_summary
    ), as.character(NA))) |>
    dplyr::select(entrezgene, gene_summary_ncbi)

  ensembl2entrez <- transcript_xref_db |>
    dplyr::filter(property == "ensembl_gene_id") |>
    dplyr::rename(ensembl_gene_id = value) |>
    dplyr::select(entrezgene, ensembl_gene_id) |>
    dplyr::distinct()


  ## exclude ensembl gene id's which are mapped to
  ## an entrezgene, but where there is a different ensembl
  ## gene id that has a match in Open Targets
  opentargetEns2Entrez <- opentarget_associations |>
    dplyr::inner_join(ensembl2entrez, by = "ensembl_gene_id") |>
    dplyr::select(entrezgene, ensembl_gene_id) |>
    dplyr::mutate(opentargets = T)

  ensembl2ignore <- ensembl2entrez |>
    dplyr::left_join(opentargetEns2Entrez,
                     by = c("ensembl_gene_id",
                            "entrezgene")) |>
    dplyr::mutate(opentargets = dplyr::if_else(
      is.na(opentargets),
      as.logical(FALSE),
      as.logical(opentargets)
    )) |>
    dplyr::filter(opentargets == F) |>
    dplyr::inner_join(
      dplyr::select(opentargetEns2Entrez, entrezgene),
      by = "entrezgene",
      multiple = "all")

  gene_xref <- gene_info |>
    dplyr::select(symbol, hgnc_id, entrezgene,
                  name, gene_biotype) |>
    dplyr::left_join(
      ensembl2entrez, by = "entrezgene", multiple = "all") |>
    dplyr::left_join(
      go2entrez, by = "entrezgene", multiple = "all") |>
    dplyr::filter(!is.na(entrezgene) & !is.na(symbol)) |>
    dplyr::distinct() |>
    dplyr::mutate(
      genename =
        paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
               entrezgene,
               "' target='_blank'>",name,"</a>")) |>
    dplyr::left_join(opentarget_associations,
                     by = "ensembl_gene_id", multiple = "all") |>
    dplyr::anti_join(
      dplyr::select(
        ensembl2ignore, c("ensembl_gene_id", "entrezgene")),
      by = c("ensembl_gene_id", "entrezgene")
    ) |>
    dplyr::left_join(ts_oncogene_annotations,
                     by = "entrezgene", multiple = "all") |>
    dplyr::left_join(upsummary2entrez,
                     by = "entrezgene", multiple = "all") |>
    dplyr::left_join(ncbisummary2entrez,
                     by = "entrezgene", multiple = "all") |>
    dplyr::left_join(cancerdrugdb[['drug_per_target']][['early_phase']],
                     by = "entrezgene", multiple = "all") |>
    dplyr::left_join(cancerdrugdb[['drug_per_target']][['late_phase']],
                     by = "entrezgene", multiple = "all") |>
    dplyr::left_join(cancerdrugdb[['approved_per_target']],
                     by = "entrezgene", multiple = "all") |>
    dplyr::mutate(tumor_suppressor = dplyr::if_else(
      is.na(tumor_suppressor),
      FALSE,
      as.logical(tumor_suppressor))) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      is.na(oncogene),
      FALSE,
      as.logical(oncogene))) |>

    dplyr::mutate(tsg_confidence_level = dplyr::if_else(
      tumor_suppressor == FALSE,
      "NONE/LIMITED",
      as.character(tsg_confidence_level))) |>
    dplyr::mutate(oncogene_confidence_level = dplyr::if_else(
      oncogene == FALSE,
      "NONE/LIMITED",
      as.character(oncogene_confidence_level))) |>
    dplyr::mutate(cancer_driver = dplyr::if_else(
      is.na(cancer_driver),
      FALSE,
      as.logical(cancer_driver))) |>
    dplyr::mutate(num_go_terms = dplyr::if_else(
      is.na(num_go_terms),
      as.integer(0),
      as.integer(num_go_terms)
    )) |>
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
    )) |>
    dplyr::select(-c(gene_summary_ncbi, gene_summary_uniprot)) |>

    dplyr::mutate(has_gene_summary = dplyr::if_else(
      !is.na(gene_summary),TRUE,FALSE
    )) |>
    dplyr::left_join(otdb$max_site_rank,
                     by = "ensembl_gene_id", multiple = "all") |>
    dplyr::mutate(cancer_max_rank = dplyr::if_else(
      is.na(cancer_max_rank), 0, as.numeric(cancer_max_rank)
    )) |>
    assign_unknown_function_rank()
    #remove_duplicate_ensembl_genes()

  saveRDS(gene_xref, file = rds_fname)
  return(gene_xref)

}
#
get_tissue_celltype_specificity <- function(
  raw_db_dir = NULL){

    ## https://www.proteinatlas.org/about/assays+annotation#gtex_rna
    ##
    ## Consensus transcript expression levels summarized per gene in
    ## 55 tissues based on transcriptomics data from HPA and GTEx.
    ## The tab-separated file includes Ensembl
    ## gene identifier ("Gene"), analysed sample ("Tissue"),
    ## transcripts per million ("TPM"), protein-transcripts
    ## per million ("pTPM") and normalized expression ("nTPM").

    tissue_cell_expr <- list()
    tissue_cell_expr[['tissue']] <- list()
    tissue_cell_expr[['tissue']][['expr_df']] <- data.frame()
    tissue_cell_expr[['tissue']][['unit']] <- NULL
    tissue_cell_expr[['tissue']][['te_SE']] <- NULL
    tissue_cell_expr[['tissue']][['te_df']] <- data.frame()

    ## https://www.proteinatlas.org/about/assays+annotation#singlecell_rna
    ##
    ## Transcript expression levels summarized per gene in 79 cell types
    ## types from 30 studies. The tab-separated file includes Ensembl
    ## gene identifier ("Gene"), gene name ("Gene name"), analysed
    ## sample ("Cell type") and normalized expresion ("nTPM").

    tissue_cell_expr[['single_cell']] <- list()
    tissue_cell_expr[['single_cell']][['expr_df']] <- data.frame()
    tissue_cell_expr[['single_cell']][['unit']] <- NULL
    tissue_cell_expr[['single_cell']][['te_SE']] <- NULL
    tissue_cell_expr[['single_cell']][['te_df']] <- data.frame()

    for(t in c("rna_tissue_consensus",
               "rna_single_cell_type")){

      local_hpa_file <- file.path(
        raw_db_dir,
        "hpa",
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
      max_types_in_group <- 5

      data <- readr::read_tsv(
        local_hpa_file,
        show_col_types = F,
        col_names = T)

      num_tissues_per_gene <- data |>
        dplyr::group_by(Gene) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::filter(n == 1)

      data <- data |>
        dplyr::anti_join(num_tissues_per_gene,
                         by = "Gene")
      data <- as.data.frame(
        data[, c(1,3,4)] |>
          purrr::set_names(
            c("ensembl_gene_id","category","exp1")) |>
          dplyr::mutate(
            category =
              stringr::str_replace_all(category," |, ","_")) |>
          dplyr::group_by(ensembl_gene_id, category) |>
          dplyr::summarise(exp = mean(exp1), .groups = "drop") |>
          tidyr::pivot_wider(names_from = category,
                             values_from = exp)
      )
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

      if(t == "rna_tissue_consensus"){
        tissue_cell_expr[['tissue']][['expr_df']] <- data
        tissue_cell_expr[['tissue']][['te_SE']] <-
          te_gene_retrieval_se
        tissue_cell_expr[['tissue']][['unit']] <- 'nTPM'
        tissue_cell_expr[['tissue']][['te_df']] <-
          as.data.frame(
            as.data.frame(
              SummarizedExperiment::assay(
                tissue_cell_expr[['tissue']][['te_SE']])
            ) |>
              dplyr::rename(ensembl_gene_id = Gene,
                            category = Group) |>
              dplyr::group_by(ensembl_gene_id,
                              category) |>
              dplyr::summarise(
                tissue = paste(unique(Tissue), collapse=", "),
                .groups = "drop")
          ) |>
          dplyr::mutate(
            category = dplyr::case_when(
              category == "Expressed-In-All" ~ "Low tissue specificity",
              category == "Group-Enriched" ~ "Group enriched",
              category == "Tissue-Enhanced" ~ "Tissue enhanced",
              category == "Tissue-Enriched" ~ "Tissue enriched",
              category == "Not-Expressed" ~ "Not detected",
              TRUE ~ as.character("Mixed"))
          ) |>
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
            ) |>
              dplyr::rename(ensembl_gene_id = Gene,
                            category = Group) |>
              dplyr::group_by(ensembl_gene_id,
                              category) |>
              dplyr::summarise(
                cell_type = paste(unique(Tissue), collapse=", "),
                .groups = "drop")
          ) |>
          dplyr::mutate(
            category = dplyr::case_when(
              category == "Expressed-In-All" ~ "Low cell type specificity",
              category == "Group-Enriched" ~ "Group enriched",
              category == "Tissue-Enhanced" ~ "Cell type enhanced",
              category == "Tissue-Enriched" ~ "Cell type enriched",
              category == "Not-Expressed" ~ "Not detected",
              TRUE ~ as.character("Mixed"))
          ) |>
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
    cache_dir = NA,
    ot_associations = NULL){


  phenotype_auxiliary_maps <- phenOncoX::get_aux_maps(
    cache_dir = cache_dir
  )

  phenotype_cancer_map <- phenOncoX::get_terms(
    cache_dir = cache_dir
  )$records

  umls_concept_map <-
    phenotype_auxiliary_maps$records$umls$concept

  efo_name_map <-
    phenotype_auxiliary_maps$records$efo$efo2name

  phenotype_cancer_efo <-
    umls_concept_map |>
    dplyr::filter(main_term == T) |>
    dplyr::select(cui, cui_name) |>
    dplyr::distinct() |>
    dplyr::inner_join(
      dplyr::select(phenotype_cancer_map,
                    efo_id, cui,
                    cui_name, primary_site),
      by = c("cui", "cui_name"), multiple = "all") |>
    dplyr::rename(disease_efo_id = efo_id) |>
    dplyr::filter(!is.na(disease_efo_id)) |>
    dplyr::mutate(cancer_phenotype = TRUE) |>
    dplyr::distinct()

  phenotypes_non_primary <- phenotype_cancer_efo |>
    dplyr::filter(is.na(primary_site))

  phenotypes_primary <- phenotype_cancer_efo |>
    dplyr::filter(!is.na(primary_site))

  phenotypes_non_primary <- phenotypes_non_primary |>
    dplyr::anti_join(phenotypes_primary, by = c("disease_efo_id","cui")) |>
    dplyr::select(-primary_site)

  ### Hereditary breast/ovarian cancer syndrome
  df <- data.frame(primary_site = 'Breast', cui = 'C0677776', stringsAsFactors = F)
  df <- dplyr::bind_rows(df, data.frame(primary_site = 'Ovary/Fallopian Tube', cui = 'C0677776', stringsAsFactors = F))

  phenotypes_hered_brca_ov <- phenotypes_non_primary |>
    dplyr::left_join(df, by = "cui", multiple = "all") |>
    dplyr::filter(!is.na(primary_site))
  ###

  phenotypes_other <- phenotypes_non_primary |>
    dplyr::anti_join(phenotypes_hered_brca_ov, by = "cui") |>
    dplyr::mutate(primary_site = dplyr::case_when(
      stringr::str_detect(tolower(cui_name), "breast") ~ "Breast",
      stringr::str_detect(tolower(cui_name), "meningioma") ~ "CNS/Brain",
      stringr::str_detect(tolower(cui_name), "cancer-predisposing") ~ "Other/Unknown",
      stringr::str_detect(tolower(cui_name), "ovarian|fallopian tube") ~ "Ovary/Fallopian Tube",
      stringr::str_detect(tolower(cui_name), "wilms|papillary renal cell|renal cell carcinoma|renal cell cancer") ~ "Kidney",
      stringr::str_detect(tolower(cui_name), "pancreatic") ~ "Pancreas",
      stringr::str_detect(tolower(cui_name), "lynch|polyposis") ~ "Colon/Rectum",
      stringr::str_detect(tolower(cui_name), "hereditary gastric") ~ "Esophagus/Stomach",
      stringr::str_detect(tolower(cui_name), "prostate") ~ "Prostate",
      stringr::str_detect(tolower(cui_name), "xeroderma|melanoma") ~ "Skin",
      stringr::str_detect(tolower(cui_name), "retinoblastoma") ~ "Eye",
      stringr::str_detect(tolower(cui_name), "medullary thyroid|follicular thyroid|papillary thyroid") ~ "Thyroid",
      TRUE ~ as.character(NA)
    ))


  all_cancer_phenotypes_efo <-
    phenotypes_primary |>
    dplyr::bind_rows(phenotypes_hered_brca_ov) |>
    dplyr::bind_rows(phenotypes_other)

  otdb_tmp <- as.data.frame(
    dplyr::select(ot_associations,
                  ot_association,
                  #symbol,
                  ensembl_gene_id) |>
      dplyr::filter(!is.na(ot_association)) |>
      tidyr::separate_rows(ot_association, sep="&") |>
      tidyr::separate(
        ot_association,
        sep=":",
        c('disease_efo_id',
          'direct_ot_association',
          'ot_datatype_support',
          'ot_association_score')) |>
      dplyr::mutate(ot_association_score =
                      as.numeric(ot_association_score)) |>
      dplyr::filter(!is.na(disease_efo_id)) |>

      dplyr::mutate(
        disease_efo_id = stringr::str_replace(disease_efo_id,"_",":")) |>
      dplyr::left_join(
        dplyr::select(efo_name_map,
                      efo_id, efo_name),
        by = c("disease_efo_id" = "efo_id"), multiple = "all") |>
      ## exclude non-specific cancer associations
      dplyr::filter(
        !stringr::str_detect(
          tolower(efo_name),
          "^(((urogenital|mixed|papillary|hereditary|familial|susceptibility|syndrome|small cell|clear cell|squamous cell|digestive system|endocrine|metastatic|invasive|pharynx|vascular|intestinal|cystic|mucinous|epithelial|benign) )?(neoplasm|carcinoma|adenocarcinoma|cancer))$")
      ) |>
      dplyr::filter(
        !stringr::str_detect(
          tolower(efo_name)," neoplasm$"
        )
      ) |>
      dplyr::left_join(
        dplyr::select(all_cancer_phenotypes_efo,
                      disease_efo_id,
                      primary_site,
                      cancer_phenotype),
        by=c("disease_efo_id"), multiple = "all") |>
      dplyr::mutate(
        cancer_phenotype =
          dplyr::if_else(
            is.na(cancer_phenotype) &
              stringr::str_detect(
                tolower(efo_name),
                "carcinoma|tumor|cancer|neoplasm|melanoma|myeloma|hemangioma|astrocytoma|leiomyoma|leukemia|lymphoma|barrett|glioma|sarcoma|blastoma|teratoma|seminoma"),
            TRUE, as.logical(cancer_phenotype))) |>
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
            TRUE, as.logical(cancer_phenotype))) |>
      dplyr::mutate(
        ot_link = paste0(
          "<a href='https://platform.opentargets.org/evidence/",
          ensembl_gene_id,"/",
          stringr::str_replace(disease_efo_id,":","_"),
          "' target=\"_blank\">", stringr::str_to_title(efo_name),"</a>")) |>
      dplyr::distinct() |>
      dplyr::distinct()
  )

  set1 <- otdb_tmp |>
    dplyr::filter(!is.na(primary_site)) |>
    dplyr::distinct()

  set2 <- otdb_tmp |>
    dplyr::filter(is.na(primary_site) & cancer_phenotype == T) |>
    dplyr::select(-c(efo_name, primary_site)) |>
    dplyr::anti_join(
      dplyr::select(set1, disease_efo_id, ensembl_gene_id),
      by = c("disease_efo_id","ensembl_gene_id")) |>
    dplyr::inner_join(
      dplyr::select(efo_name_map,
                    efo_id, efo_name, primary_site),
      by = c("disease_efo_id" = "efo_id"), multiple = "all") |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(primary_site))

  set12 <- dplyr::bind_rows(set1, set2)

  set3 <- otdb_tmp |>
    dplyr::anti_join(
      dplyr::select(set12, disease_efo_id, ensembl_gene_id),
      by = c("disease_efo_id","ensembl_gene_id")) |>
    dplyr::distinct()

  otdb_all <- set12 |>
    dplyr::bind_rows(set3) |>
    dplyr::arrange(primary_site, ensembl_gene_id,
                   desc(ot_association_score))

  otdb_tissue_scores <- as.data.frame(
    otdb_all |>
      dplyr::filter(!is.na(primary_site)) |>
      dplyr::group_by(primary_site, ensembl_gene_id) |>
      dplyr::summarise(
        tissue_assoc_score =
          round(sum(ot_association_score), digits = 6),
        .groups = "drop"
      ) |>
      dplyr::arrange(
        primary_site, desc(tissue_assoc_score))
  )

  otdb_site_rank <- data.frame()
  for(s in unique(otdb_tissue_scores$primary_site)){
    tmp <- otdb_tissue_scores |>
      dplyr::filter(primary_site == s) |>
      dplyr::mutate(
        tissue_assoc_rank =
          round(dplyr::percent_rank(tissue_assoc_score),
                digits = 6)
      )
    otdb_site_rank <- otdb_site_rank |>
      dplyr::bind_rows(tmp)
  }

  otdb_global_rank <- as.data.frame(
    otdb_site_rank |>
      dplyr::group_by(.data$ensembl_gene_id) |>
      dplyr::summarise(
        global_assoc_score = sum(
          .data$tissue_assoc_rank),
        .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::arrange(
        dplyr::desc(.data$global_assoc_score)) |>
      dplyr::mutate(
        global_assoc_rank = round(
          dplyr::percent_rank(.data$global_assoc_score),
          digits = 5))
  )

  otdb_max_site_rank <- as.data.frame(
    otdb_site_rank |>
      dplyr::arrange(desc(tissue_assoc_rank)) |>
      dplyr::select(ensembl_gene_id, tissue_assoc_rank) |>
      dplyr::group_by(ensembl_gene_id) |>
      dplyr::summarise(
        cancer_max_rank = max(tissue_assoc_rank),
        .groups = "drop")
  )

  otdb <- list()
  otdb$all <- otdb_all
  otdb$gene_rank <- otdb_site_rank |>
    dplyr::left_join(
      otdb_global_rank, by = "ensembl_gene_id", multiple = "all")

  #otdb$site_rank <- otdb_site_rank
  otdb$max_site_rank <- otdb_max_site_rank

  return(otdb)

}

get_ligand_receptors <- function(
  raw_db_dir = NULL,
  keggdb = NULL,
  update = F){

  ligand_receptordb <- list()
  ligand_receptordb[['cellchatdb']] <- NULL
  ligand_receptordb[['celltalkdb']] <- NULL

  for(db in c('cellchatdb','celltalkdb')){

    if(db == 'cellchatdb'){
      rds_fname <- file.path(
        raw_db_dir,
        "cellchatdb",
        "cellchatdb_interactions.rds")

      if(update == F & file.exists(rds_fname)){
        ligand_receptordb[[db]] <- readRDS(file=rds_fname)
        next
      }

      ligand_receptor_db <- CellChat::CellChatDB.human$interaction |>
        magrittr::set_rownames(NULL) |>
        dplyr::rename(interaction_members = interaction_name) |>
        dplyr::rename(interaction_name = interaction_name_2) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC: 4393358", "PMID: 25772309"
        )) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC4571854", "PMID: 26124272"
        )) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC2194005", "PMID: 12208882"
        )) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC2194217", "PMID: 14530377"
        )) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC1237098", "PMID: 16093349"
        )) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,"PMC2431087", "PMID: 18250165"
        )) |>
        dplyr::mutate(evidence = stringr::str_squish(
          stringr::str_replace_all(
            evidence,
            ",","; "))
        ) |>
        dplyr::mutate(evidence = stringr::str_replace_all(
          evidence,
          " ","")) |>
        dplyr::mutate(interaction_id = dplyr::row_number()) |>
        dplyr::select(interaction_id, interaction_name,
                      annotation, pathway_name, dplyr::everything())


      ligand_receptor_kegg_support <- ligand_receptor_db |>
        dplyr::select(interaction_id, evidence) |>
        tidyr::separate_rows(evidence, sep=";") |>
        dplyr::filter(stringr::str_detect(evidence,"KEGG")) |>
        dplyr::mutate(evidence = stringr::str_replace(
          evidence, "KEGG:",""
        )) |>
        dplyr::left_join(keggdb$TERM2NAME,
                         by = c("evidence" = "standard_name"), multiple = "all") |>
        dplyr::mutate(ligand_receptor_kegg_support = paste0(
          "<a href='https://www.genome.jp/pathway/",
          stringr::str_replace(evidence,"hsa","map"),
          "' target='_blank'>",name,"</a>")) |>
        dplyr::select(-evidence)

      ligand_receptor_literature <- as.data.frame(
        ligand_receptor_db |>
          dplyr::select(interaction_id, evidence) |>
          tidyr::separate_rows(evidence, sep=";") |>
          dplyr::filter(stringr::str_detect(evidence,"PMID")) |>
          dplyr::mutate(evidence = stringr::str_replace(
            evidence, "PMID:",""
          )) |>
          dplyr::rename(pmid = evidence)
      )

      ligand_receptor_citations <- readRDS(
        file.path(
          raw_db_dir,
          "cellchatdb",
          "cellchatdb_citations.rds"))

      #ligand_receptor_citations <- get_citations_pubmed(
      #  pmid = unique(ligand_receptor_literature$pmid),
      #  chunk_size = 10)

      ligand_receptor_literature <- as.data.frame(
        ligand_receptor_literature |>
          dplyr::mutate(pmid = as.integer(pmid)) |>
          dplyr::left_join(ligand_receptor_citations,
                           by = c("pmid" = "pmid"), multiple = "all") |>
          dplyr::filter(!is.na(link)) |>
          dplyr::group_by(interaction_id) |>
          dplyr::summarise(
            ligand_receptor_literature_support = paste(
              link, collapse= ","),
            ligand_receptor_literature = paste(
              citation, collapse="|"
            )
          )
      )

      ligand_receptor_db_final <- ligand_receptor_db |>
        dplyr::left_join(dplyr::select(
          ligand_receptor_kegg_support, interaction_id,
          ligand_receptor_kegg_support), by = "interaction_id", multiple = "all") |>
        dplyr::left_join(dplyr::select(
          ligand_receptor_literature, interaction_id,
          ligand_receptor_literature_support), multiple = "all") |>
        dplyr::mutate(literature_support = dplyr::if_else(
          is.na(ligand_receptor_kegg_support),
          as.character(ligand_receptor_literature_support),
          paste(ligand_receptor_kegg_support,
                ligand_receptor_literature_support,
                sep = ", ")
        )) |>
        dplyr::mutate(literature_support = stringr::str_replace_all(
          literature_support,",?NA$",""
        )) |>
        dplyr::select(-c(ligand_receptor_literature_support,
                         ligand_receptor_kegg_support,
                         evidence))

      interaction2ligands <- as.data.frame(
        ligand_receptor_db_final |>
          dplyr::select(interaction_id, ligand) |>
          dplyr::rename(symbol = ligand) |>
          dplyr::mutate(class = "ligand")
      )

      interaction2receptor <- as.data.frame(
        ligand_receptor_db_final |>
          dplyr::select(interaction_id, receptor) |>
          dplyr::rename(symbol = receptor) |>
          tidyr::separate_rows(symbol, sep = "_") |>
          dplyr::mutate(symbol = dplyr::if_else(
            stringr::str_detect(symbol,"^(R2|TGFbR2)$"),
            "TGFBR2",
            as.character(symbol)
          )) |>
          dplyr::mutate(class = "receptor") |>
          dplyr::mutate(symbol = stringr::str_replace(
            symbol,"b","B"
          ))
      )

      cellchatdb <- list()
      cellchatdb[['db']] <- ligand_receptor_db_final
      cellchatdb[['xref']] <- interaction2ligands |>
        dplyr::bind_rows(interaction2receptor)

      saveRDS(cellchatdb, file = rds_fname)

      ligand_receptordb[[db]] <- cellchatdb
    }

    if(db == "celltalkdb"){

      rds_fname <- file.path(
        raw_db_dir,
        "celltalkdb",
        "celltalkdb_interactions.rds")

      if(update == F & file.exists(rds_fname)){
        celltalkdb <- readRDS(file=rds_fname)
        next
      }

      celltalkdb_raw <- readr::read_tsv(
        file = file.path(
          raw_db_dir,
          "celltalkdb",
          "human_lr_pair.txt"),
        col_types = "cccddccccc",
        show_col_types = F) |>
        dplyr::select(ligand_gene_id,
                      receptor_gene_id,
                      evidence) |>
        dplyr::mutate(interaction_id = dplyr::row_number()) |>
        dplyr::rename(entrezgene_ligand = ligand_gene_id,
                      entrezgene_receptor = receptor_gene_id,
                      pmid = evidence) |>
        dplyr::mutate(pmid = stringr::str_replace(
          pmid, ",321961154", ""))


      celltalkdb_literature <- as.data.frame(
        celltalkdb_raw |>
          dplyr::select(interaction_id, pmid) |>
          tidyr::separate_rows(pmid, sep=",") |>
        dplyr::filter(nchar(pmid) > 2)
      )

      ligand_receptor_citations <- readRDS(
        file.path(
          raw_db_dir,
          "celltalkdb",
          "celltalkdb_citations.rds"))
      # ligand_receptor_citations <- get_citations_pubmed(
      #   pmid = unique(celltalkdb_literature$pmid),
      #   chunk_size = 100)

      celltalkdb_literature <- as.data.frame(
        celltalkdb_literature |>
          dplyr::mutate(pmid = as.integer(pmid)) |>
          dplyr::left_join(ligand_receptor_citations,
                           by = c("pmid" = "pmid"), multiple = "all") |>
          dplyr::filter(!is.na(link)) |>
          dplyr::group_by(interaction_id) |>
          dplyr::summarise(
            ligand_receptor_literature_support = paste(
              link, collapse= ","),
            ligand_receptor_literature = paste(
              citation, collapse="|"
            )
          )
      )

      celltalkdb_final <- as.data.frame(
        celltalkdb_raw |>
          dplyr::select(entrezgene_ligand,
                        entrezgene_receptor,
                        interaction_id) |>
          dplyr::left_join(
            dplyr::select(
              celltalkdb_literature,
              interaction_id,
              ligand_receptor_literature_support,
              ligand_receptor_literature),
            by = "interaction_id", multiple = "all") |>
          dplyr::rename(
            literature_support = ligand_receptor_literature_support) |>
          dplyr::select(
            interaction_id,
            entrezgene_ligand,
            entrezgene_receptor,
            dplyr::everything()
          )
      )

      ligand_receptordb[[db]] <- celltalkdb_final

      saveRDS(celltalkdb_final, file = rds_fname)

    }
  }

  return(ligand_receptordb)
}

get_hpa_associations <- function(
  raw_db_dir = NULL,
  gene_xref = NULL,
  update = F){

  assertable::assert_colnames(
    gene_xref,
    c("symbol", "ensembl_gene_id"),
    only_colnames = F,
    quiet = T
  )

  rds_fname <- file.path(
    raw_db_dir,
    "hpa",
    "hpa.rds")

  if(update == F & file.exists(rds_fname)){
    hpa <- readRDS(file = rds_fname)
    return(hpa)
  }

  variables_json <- c('RNA tissue specificity',
                      'RNA tissue distribution',
                      'RNA tissue specific nTPM',
                      'RNA single cell type specificity',
                      'RNA single cell type distribution',
                      'RNA single cell type specific nTPM',
                      'RNA cancer specificity',
                      'RNA cancer distribution',
                      'RNA cancer specific FPKM',
                      'RNA cell line specificity',
                      'RNA cell line distribution',
                      'RNA cell line specific nTPM',
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
                    'rna_tissue_specific_nTPM',
                    'rna_single_cell_type_specificity',
                    'rna_single_cell_type_distribution',
                    'rna_single_cell_type_specificity_nTPM',
                    'rna_cancer_specificity',
                    'rna_cancer_distribution',
                    'rna_cancer_specific_FPKM',
                    'rna_cell_line_specificity',
                    'rna_cell_line_distribution',
                    'rna_cell_line_specific_nTPM',
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

  hpa <- data.frame()
  k <- 1

  for (i in 1:nrow(gene_xref)) {
    g <- gene_xref[i,]

    if(g$gene_biotype != "protein-coding"){
      next
    }
    if(is.na(g$ensembl_gene_id)){
      next
    }

    if(!startsWith(g$ensembl_gene_id, "ENSG") |
       stringr::str_detect(g$ensembl_gene_id, "PAR")){
      next
    }

    local_json <-
      file.path(
        raw_db_dir,
        "hpa",
        "json",
        paste0(g$ensembl_gene_id,".json"))

    # if(!file.exists(local_json)){
    #   protein_atlas_url <-
    #     paste0("https://www.proteinatlas.org/search/",
    #            g$ensembl_gene_id,"?format=json")
    #
    #   if(RCurl::url.exists(protein_atlas_url) == TRUE){
    #     cat(k, ' - ', protein_atlas_url, '\n')
    #     k <- k + 1
    #
    #     download.file(
    #       protein_atlas_url,
    #       destfile = local_json,
    #       quiet = T)
    #     Sys.sleep(4)
    #   }
    # }


    if(file.exists(local_json)){
      pa_data <- jsonlite::fromJSON(txt = local_json)

      for(j in 1:length(variables_json)){
        var_name_json <- variables_json[j]
        var_name_df <- variables_df[j]
        if(!is.null(pa_data[[var_name_json]])){
          if(!startsWith(var_name_df, "pathology") &
             !stringr::str_detect(var_name_df,"(TPM|FPKM|rrid)$")){
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

  hpa <- hpa |>
    dplyr::mutate(
      value =  dplyr::if_else(
        property == "antibody_rrid",
        stringr::str_replace_all(value, "^(NA\\|){1,}|(\\|NA){1,}$",""),
        as.character(value)
      )
    ) |>
    dplyr::mutate(
      value =  dplyr::if_else(
        property == "antibody_rrid",
        stringr::str_replace_all(value, "(\\|NA\\|)","|"),
        as.character(value)
      )
    ) |>
    dplyr::mutate(
      value =  dplyr::if_else(
        property == "antibody_rrid" & value == "NA",
        as.character(''),
        as.character(value)
      )
    )


  saveRDS(hpa, file = rds_fname)
  return(hpa)


}


get_mean_median_tpm <- function(
  tpm_matrix,
  detectable_tpm_level = 1,
  tumor_code = "BRCA_1",
  filter_on = "median"){

  tpm_pr_gene <- as.data.frame(
    data.frame(
      symbol = rownames(tpm_matrix),
      median_tpm_exp = round(matrixStats::rowMedians(
        tpm_matrix), digits = 4
      ),
      mean_tpm_exp = round(matrixStats::rowMeans2(
        tpm_matrix), digits = 4
      ),
      tumor = tumor_code,
      stringsAsFactors = F
    ) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(
        median_tpm_exp = round(mean(median_tpm_exp), digits = 4),
        mean_tpm_exp = round(mean(mean_tpm_exp), digits = 4),
        tumor = paste(unique(tumor), collapse=","),
        .groups = "drop") |>
      dplyr::filter(median_tpm_exp >= detectable_tpm_level)) |>
    dplyr::mutate(median_tpm_subtype = paste(
      tumor, median_tpm_exp, sep=":"
    )) |>
    dplyr::select(symbol, median_tpm_subtype)


  return(tpm_pr_gene)


}


get_tcga_db <- function(
  raw_db_dir = NULL,
  gene_xref = NULL,
  tcga_release = "release36_20221212",
  update = F){

  #tcga_release <- 'release34_20220727'

  rds_fname <- file.path(
    raw_db_dir,
    "tcga",
    "tcgadb.rds")

  rnaseq_path <- file.path(
    raw_db_dir,
    "tcga",
    "current_release",
    "rnaseq"
  )

  coexpression_tsv <- file.path(
    raw_db_dir,
    "tcga",
    paste0(
      "co_expression_strong_moderate.",
      tcga_release,
      ".tsv.gz"))

  tcga_clinical_rds <- file.path(
    raw_db_dir,
    "tcga",
    "current_release",
    "clinical",
    "tcga_clinical.rds"
  )

  tcga_aberration_rds <- file.path(
    raw_db_dir,
    "tcga",
    "current_release",
    "gene",
    "tcga_gene_aberration_rate.rds"
  )

  maf_codes_tsv <- file.path(
    raw_db_dir,
    "maf_codes.tsv"
  )
  maf_path <- file.path(
    raw_db_dir,
    "tcga",
    "current_release",
    "snv_indel")


  recurrent_variants_tsv <- file.path(
    raw_db_dir,
    "tcga",
    "current_release",
    "vcf",
    "tcga_recurrent_coding_gvanno_grch38.tsv.gz"
  )

  # recurrent_variants_tsv <-
  #   "~/project_data/analysis__tcga/tcga/output/release36_20221212/vcf/tcga_recurrent_coding_gvanno_grch38.tsv.gz"

  pfam_domains_tsv <- file.path(
    raw_db_dir,
    "pfam",
    "pfam.domains.tsv.gz"
  )

  if(update == F & file.exists(rds_fname)){
    tcgadb <- readRDS(file = rds_fname)
    return(tcgadb)
  }


  tcga_clinical <-
    readRDS(file = tcga_clinical_rds)
  tcga_aberration_stats <-
    readRDS(file = tcga_aberration_rds) |>
    dplyr::filter(
      clinical_strata == "site" |
        (clinical_strata == "site_diagnosis" &
           percent_mutated >= 1))

  tcga_aberration_stats$genomic_strata <- NULL
  tcga_aberration_stats$entrezgene <- NULL
  tcga_aberration_stats$consensus_calls <- NULL
  tcga_aberration_stats$fp_driver_gene <- NULL

  site_code <- tcga_aberration_stats |>
    dplyr::select(primary_site) |>
    dplyr::distinct() |>
    dplyr::arrange(primary_site) |>
    dplyr::mutate(site_code = dplyr::row_number())

  diagnosis_code <- tcga_aberration_stats |>
    dplyr::select(primary_diagnosis_very_simplified) |>
    dplyr::distinct() |>
    dplyr::rename(
      primary_diagnosis = primary_diagnosis_very_simplified) |>
    dplyr::arrange(primary_diagnosis) |>
    dplyr::mutate(
      diagnosis_code = dplyr::row_number())

  clinical_strata_code <- tcga_aberration_stats |>
    dplyr::select(clinical_strata) |>
    dplyr::distinct() |>
    dplyr::mutate(
      clinical_strata_code = dplyr::row_number())


  tcga_aberration_stats <- tcga_aberration_stats |>
    dplyr::rename(primary_diagnosis = primary_diagnosis_very_simplified) |>
    dplyr::left_join(clinical_strata_code, by = "clinical_strata",
                     multiple = "all") |>
    dplyr::left_join(diagnosis_code, by = "primary_diagnosis",  multiple = "all") |>
    dplyr::left_join(site_code, by = "primary_site",  multiple = "all") |>
    dplyr::select(-c(primary_site, primary_diagnosis,
                     clinical_strata))

  maf_codes <- read.table(file = maf_codes_tsv,
                          header = T, sep = "\t", quote ="")

  maf_datasets <- list()
  i <- 1
  while(i <= nrow(maf_codes)){
    primary_site <- maf_codes[i,]$primary_site
    maf_code <- maf_codes[i,]$code
    maf_file <- file.path(
      maf_path, paste0(
        "tcga_mutation_grch38_",
        maf_code,"_0.maf.gz"))
    if(file.exists(maf_file)){
      tmp <- read.table(gzfile(maf_file), quote="",
                        header = T, stringsAsFactors = F,
                        sep="\t", comment.char="#")
      tmp$primary_site <- NULL
      tmp$site_diagnosis_code <- NULL
      tmp$Tumor_Sample_Barcode <- stringr::str_replace(
        tmp$Tumor_Sample_Barcode,"-[0-9][0-9][A-Z]$","")

      clinical <- tcga_clinical$slim |>
        dplyr::filter(primary_site == primary_site) |>
        dplyr::select(bcr_patient_barcode, primary_diagnosis_very_simplified,
                      MSI_status, Gleason_score, ER_status,
                      PR_status, HER2_status,
                      pancan_subtype_selected) |>
        dplyr::rename(Diagnosis = primary_diagnosis_very_simplified,
                      Tumor_Sample_Barcode = bcr_patient_barcode,
                      PanCancer_subtype = pancan_subtype_selected) |>
        dplyr::semi_join(tmp, by = "Tumor_Sample_Barcode",  multiple = "all")

      write.table(tmp, file="tmp.maf", quote=F, row.names = F,
                  col.names = T, sep ="\t")
      system('gzip -f tmp.maf', intern = T)
      maf <- maftools::read.maf("tmp.maf.gz", verbose = F, clinicalData = clinical)
      maf_datasets[[maf_code]] <- maf
      #saveRDS(maf, file = file.path(raw_db_dir, "data-raw", "maf", paste0(maf_code,".maf.rds")))
    }
    i <- i + 1
    cat(primary_site,'\n')

  }
  system('rm -f tmp.maf.gz')


  pfam_domains <- as.data.frame(
    readr::read_tsv(
      pfam_domains_tsv,

      show_col_types = F
    )) |>
    dplyr::rename(PFAM_ID = pfam_id, PROTEIN_DOMAIN = url) |>
    dplyr::rename(PFAM_DOMAIN_NAME = name) |>
    dplyr::select(PFAM_ID, PFAM_DOMAIN_NAME) |>
    dplyr::mutate(PFAM_ID = stringr::str_replace(
      PFAM_ID,"\\.[0-9]{1,}$","")
    )

  ##TCGA recurrent SNVs/InDels
  recurrent_tcga_variants <- as.data.frame(readr::read_tsv(
    file = recurrent_variants_tsv,
    skip = 1, na = c("."), show_col_types = F) |>
      dplyr::select(TCGA_SITE_RECURRENCE,
                    CHROM, POS, REF, ALT,
                    TCGA_TOTAL_RECURRENCE,
                    PFAM_DOMAIN, HGVSc, LoF,
                    MUTATION_HOTSPOT,
                    #MUTATION_HOTSPOT_NEW,
                    #MUTATION_HOTSPOT_TTYPE_NEW,
                    MUTATION_HOTSPOT_MATCH,
                    HGVSp_short,
                    ENSEMBL_TRANSCRIPT_ID,
                    SYMBOL,
                    COSMIC_MUTATION_ID,
                    Consequence,
                    AMINO_ACID_START,
                    VEP_ALL_CSQ) |>
      dplyr::mutate(VAR_ID = paste(
        CHROM, POS, REF, ALT, sep = "_")
      ) |>
      dplyr::rename(PFAM_ID = PFAM_DOMAIN) |>
      dplyr::mutate(LOSS_OF_FUNCTION = FALSE) |>
      dplyr::mutate(LOSS_OF_FUNCTION = dplyr::if_else(
        !is.na(LoF) & LoF == "HC",
        as.logical(TRUE),
        as.logical(LOSS_OF_FUNCTION),
      )) |>
      dplyr::rename(CONSEQUENCE = Consequence,
                    TOTAL_RECURRENCE = TCGA_TOTAL_RECURRENCE,
                    PROTEIN_CHANGE = HGVSp_short) |>
      dplyr::select(SYMBOL,
                    VAR_ID,
                    CONSEQUENCE,
                    PROTEIN_CHANGE,
                    AMINO_ACID_START,
                    PFAM_ID,
                    HGVSc,
                    MUTATION_HOTSPOT,
                    #MUTATION_HOTSPOT_NEW,
                    #MUTATION_HOTSPOT_TTYPE_NEW,
                    MUTATION_HOTSPOT_MATCH,
                    LOSS_OF_FUNCTION,
                    ENSEMBL_TRANSCRIPT_ID,
                    COSMIC_MUTATION_ID,
                    TCGA_SITE_RECURRENCE,
                    TOTAL_RECURRENCE,
                    VEP_ALL_CSQ) |>
      tidyr::separate_rows(TCGA_SITE_RECURRENCE, sep=",") |>
      tidyr::separate(TCGA_SITE_RECURRENCE, into =
                        c("PRIMARY_SITE","SITE_RECURRENCE", "TCGA_SAMPLES"),
                      sep = ":",
                      remove = T) |>
      dplyr::select(-TCGA_SAMPLES) |>
      dplyr::distinct() |>
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
      ) |>
      dplyr::filter(!stringr::str_detect(
        CONSEQUENCE,"^(intron|intergenic|mature|non_coding|synonymous|upstream|downstream|3_prime|5_prime)"))
  )

  ## reduce number of properties in VEP_ALL_CSQ
  csq_all_slim <- as.data.frame(
    recurrent_tcga_variants |>
    dplyr::select(VAR_ID, VEP_ALL_CSQ) |>
    tidyr::separate_rows(VEP_ALL_CSQ, sep=",") |>
    tidyr::separate(
      VEP_ALL_CSQ, c("V1","V2","V3","V4",
                     "V5","V6","V7","V8"),
                     sep = ":") |>
    dplyr::mutate(VEP_ALL_CSQ = paste(
      V1,V7,V4,V5, sep=":"
    )) |>
    dplyr::group_by(VAR_ID) |>
    dplyr::summarise(VEP_ALL_CSQ = paste(
      unique(VEP_ALL_CSQ), collapse=", "
    ), .groups = "drop")
  )

  recurrent_tcga_variants$VEP_ALL_CSQ <- NULL
  recurrent_tcga_variants <- recurrent_tcga_variants |>
    dplyr::left_join(csq_all_slim, by = "VAR_ID",  multiple = "all") |>
    dplyr::select(
      SYMBOL, VAR_ID, CONSEQUENCE,
      PROTEIN_CHANGE, MUTATION_HOTSPOT,
      AMINO_ACID_START, PFAM_ID,
      LOSS_OF_FUNCTION, ENSEMBL_TRANSCRIPT_ID,
      MUTATION_HOTSPOT_MATCH,
      COSMIC_MUTATION_ID, PRIMARY_SITE,
      SITE_RECURRENCE, TOTAL_RECURRENCE,
      VEP_ALL_CSQ
    )


  ##TCGA co-expression data
  raw_coexpression <-
    readr::read_tsv(
      coexpression_tsv,
      col_names = c("symbol_A","symbol_B",
                    "r","p_value","tumor"),
      show_col_types = F)
  coexpression_genes1 <- raw_coexpression |>
    dplyr::rename(symbol = symbol_A,
                  symbol_partner = symbol_B) |>
    dplyr::left_join(
      dplyr::select(
        gene_xref, symbol,
        tumor_suppressor,
        oncogene, cancer_driver),
      by = "symbol",  multiple = "all") |>
    dplyr::filter(
      tumor_suppressor == T |
        oncogene == T |
        cancer_driver == T)

  coexpression_genes2 <- raw_coexpression |>
    dplyr::rename(symbol = symbol_B,
                  symbol_partner = symbol_A) |>
    dplyr::left_join(
      dplyr::select(
        gene_xref, symbol, tumor_suppressor,
        oncogene, cancer_driver),
      by = "symbol",  multiple = "all") |>
    dplyr::filter(
      tumor_suppressor == T |
        oncogene == T |
        cancer_driver == T)

  rm(raw_coexpression)

  tcga_coexp_db <- dplyr::bind_rows(
    coexpression_genes1,
    coexpression_genes2) |>
    dplyr::arrange(desc(r)) |>
    dplyr::mutate(p_value = signif(
      p_value, digits = 4)) |>
    dplyr::filter(p_value < 1e-6) |>
    dplyr::mutate(r = signif(r, digits = 3)) |>
    dplyr::filter((r <= -0.7 & r < 0) | (r >= 0.7))

  gdc_projects <- as.data.frame(TCGAbiolinks::getGDCprojects()) |>
    dplyr::filter(is.na(dbgap_accession_number) &
                    startsWith(id,"TCGA"))

  all_tcga_mean_median_tpm <- data.frame()

  for(i in 1:nrow(gdc_projects)){
    tumor_code <- gdc_projects[i,]$tumor
    rnaseq_rds_fname <-
      file.path(
        rnaseq_path,
        paste0("rnaseq_",
               tumor_code,
               ".rds")
      )
      # file.path("","Users","sigven", "project_data","analysis__tcga",
      #           "tcga", "output", "rnaseq",
      #           paste0("rnaseq_",tumor_code,"_",
      #                  tcga_release,".rds"))

    cat(tumor_code, '\n')
    if(file.exists(rnaseq_rds_fname)){
      exp_data <- readRDS(file = rnaseq_rds_fname)

      tpm_matrix <- exp_data$matrix$tpm$tumor
      rm(exp_data)

      colnames(tpm_matrix) <-
        stringr::str_replace(
          colnames(tpm_matrix),"-0[0-9][A-Z]$","")

      clinical_subtype_samples <- tcga_clinical$slim |>
        dplyr::filter(tumor == tumor_code) |>
        dplyr::filter(primary_diagnosis_very_simplified !=
                        "Other subtype(s)") |>
        dplyr::select(site_diagnosis_code, bcr_patient_barcode)

      all_tpm_data <-  get_mean_median_tpm(
        tpm_matrix,
        tumor_code = tumor_code
      )

      if(length(unique(clinical_subtype_samples$site_diagnosis_code)) > 1){

        subtypes <- unique(clinical_subtype_samples$site_diagnosis_code)

        for(subtype in subtypes){

          subtype_samples <- clinical_subtype_samples[
            clinical_subtype_samples$site_diagnosis_code == subtype,]$bcr_patient_barcode

          tpm_matrix_subtype <- tpm_matrix[
            , colnames(tpm_matrix) %in% subtype_samples
          ]

          tpm_data <- get_mean_median_tpm(
            tpm_matrix_subtype,
            tumor_code = subtype
          )

          all_tpm_data <- all_tpm_data |>
            dplyr::bind_rows(tpm_data)

        }

      }
      all_tcga_mean_median_tpm <- all_tcga_mean_median_tpm |>
        dplyr::bind_rows(all_tpm_data)
      rm(all_tpm_data)

    }
  }

  all_tcga_tpm_above_1 <- as.data.frame(
    all_tcga_mean_median_tpm |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(
        median_tpm_subtype = paste(
          sort(median_tpm_subtype), collapse="|"
        )
      )) |>
    dplyr::left_join(
      dplyr::select(gene_xref, symbol, gene_biotype),  multiple = "all"
    ) |>
    dplyr::filter(gene_biotype == "protein-coding")


  tcgadb <- list()
  tcgadb[['coexpression']] <- tcga_coexp_db
  tcgadb[['aberration']] <- tcga_aberration_stats
  tcgadb[['recurrent_variants']] <- recurrent_tcga_variants
  tcgadb[['median_ttype_expression']] <- all_tcga_tpm_above_1
  tcgadb[['pfam']] <- pfam_domains
  tcgadb[['maf_codes']] <- maf_codes
  tcgadb[['maf']] <- maf_datasets
  tcgadb[['site_code']] <- site_code
  tcgadb[['diagnosis_code']] <- diagnosis_code
  tcgadb[['clinical_strata_code']] <- clinical_strata_code

  saveRDS(tcgadb, file = rds_fname)

  return(tcgadb)

}

get_paralog_SL_predictions <- function(
  raw_db_dir = NULL,
  gene_info = NULL){


  predicted_SL_pairs_csv <- file.path(
    raw_db_dir,
    "synlethdb",
    "ryan_predicted_paralog_SL_pairs.csv"
  )

  predicted_SL_data <- readr::read_csv(
    file = predicted_SL_pairs_csv,
    show_col_types = F) |>
    dplyr::select(
      A1_entrez, A2_entrez, prediction_score,
      prediction_rank, prediction_percentile,
      min_sequence_identity, fet_ppi_overlap,
      conservation_score,
      family_size,
      has_pombe_ortholog,
      has_cerevisiae_ortholog
    ) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene),
      by = c("A1_entrez" = "entrezgene"),  multiple = "all"
    ) |>
    dplyr::rename(ppi_overlap = fet_ppi_overlap) |>
    dplyr::rename(symbol_A1 = symbol) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene),
      by = c("A2_entrez" = "entrezgene"),  multiple = "all"
    ) |>
    dplyr::rename(symbol_A2 = symbol) |>
    dplyr::mutate(
      sequence_identity_pct = round(
        min_sequence_identity * 100, digits = 3
      )
    ) |>
    dplyr::mutate(prediction_score = round(
      prediction_score, digits = 3
    )) |>
    dplyr::mutate(ppi_overlap = round(
      ppi_overlap, digits = 2
    )) |>
    dplyr::rename(entrezgene_A1 = A1_entrez,
                  entrezgene_A2 = A2_entrez) |>
    dplyr::select(entrezgene_A1,
                  symbol_A1,
                  entrezgene_A2,
                  symbol_A2,
                  prediction_score,
                  prediction_percentile,
                  sequence_identity_pct,
                  family_size,
                  conservation_score,
                  ppi_overlap,
                  has_pombe_ortholog,
                  has_cerevisiae_ortholog
                  )

  return(predicted_SL_data)


}

get_synthetic_lethality_pairs <- function(
  raw_db_dir = NULL){

  sl_pairs_csv <- file.path(
    raw_db_dir,
    "synlethdb",
    "Human_SL.csv"
  )

  sl_pairs_data_raw <- as.data.frame(
    readr::read_csv(
      sl_pairs_csv, show_col_types = F) |>
      janitor::clean_names() |>
      dplyr::mutate(entrezgene_a = as.integer(gene_a_identifier)) |>
      dplyr::mutate(entrezgene_b = as.integer(gene_b_identifier)) |>
      dplyr::rename(pmid = sl_pubmed_id) |>
      dplyr::select(entrezgene_a, entrezgene_b,
                    pmid, sl_source, sl_cell_line, sl_statistic_score) |>
      dplyr::mutate(sl_id = dplyr::row_number()) |>
      tidyr::separate_rows(pmid, sep=";") |>
      tidyr::separate_rows(pmid, sep="/") |>
      dplyr::filter(!stringr::str_detect(pmid, ",")) |>
      tidyr::separate_rows(sl_cell_line, sep=";") |>
      dplyr::distinct() |>
      dplyr::group_by(entrezgene_a, entrezgene_b, sl_id,
                      pmid, sl_source, sl_statistic_score) |>
      dplyr::summarise(sl_cell_line = paste(
        sl_cell_line, collapse=";"), .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(sl_statistic_score = round(
        as.numeric(sl_statistic_score), digits = 3)) |>
      dplyr::filter(stringr::str_detect(
        sl_source, "GenomeRNAi|CRISPR"))
  )

  pmid_data <-
    get_citations_pubmed(unique(sl_pairs_data_raw$pmid),
                         chunk_size = 100) |>
    dplyr::mutate(pmid = as.character(pmid)) |>
    dplyr::rename(sl_citation = citation,
                  sl_citation_link = link)

  sl_pairs_data <- as.data.frame(
    sl_pairs_data_raw |>
      dplyr::left_join(pmid_data, by = "pmid",  multiple = "all") |>
      dplyr::group_by(
        sl_id, entrezgene_a, entrezgene_b,
        sl_source, sl_cell_line, sl_statistic_score) |>
      dplyr::summarise(
        sl_pmid = paste(pmid, collapse=";"),
        sl_citation = paste(sl_citation, collapse="; "),
        sl_citation_link = paste(sl_citation_link, collapse= ", "),
        .groups = "drop")
  )

  return(sl_pairs_data)

}


get_subcellular_annotations <- function(
    raw_db_dir = NA,
    transcript_xref_db = NA,
    update = F){


  rds_fname <- file.path(
    raw_db_dir,
    "compartments",
    "compartmentdb.rds")

  if(update == F & file.exists(rds_fname)){
    compartmentdb <- readRDS(file = rds_fname)
    return(compartmentdb)
  }

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

  gganatogram_map_xlsx <- file.path(
    raw_db_dir,
    "compartments",
    "compartment_mapping_jw.xlsx")

  go_gganatogram_map <-
    openxlsx::read.xlsx(gganatogram_map_xlsx,
                        sheet = 2) |>
    janitor::clean_names() |>
    dplyr::select(organ, go_id) |>
    dplyr::mutate(organ = stringr::str_trim(organ)) |>
    dplyr::left_join(gganatogram::cell_key$cell, by = "organ",  multiple = "all") |>
    dplyr::select(-c(type,value,colour)) |>
    dplyr::rename(ggcompartment = organ)

  ensembl_protein_xref <- transcript_xref_db |>
    dplyr::filter(property == "ensembl_protein_id") |>
    dplyr::select(entrezgene, value) |>
    dplyr::rename(ensembl_protein_id = value) |>
    dplyr::mutate(
      ensembl_protein_id = stringr::str_replace(
        ensembl_protein_id, "\\.[0-9]{1,}$",""
      ))


  subcell_evidence <- list()
  for (e in c('experiments','knowledge','textmining')) {
    fname <- file.path(
      raw_db_dir, "compartments", paste0(
        "human_compartment_",
        e, "_full.tsv.gz"
      )
    )

    subcell_evidence[[e]] <- as.data.frame(
      readr::read_tsv(
        file = fname, col_names = F,
        show_col_types = F)
    )

    if (e == 'experiments') {
      subcell_evidence[[e]] <- subcell_evidence[[e]] |>
        dplyr::mutate(channel = "Experimental") |>
        dplyr::rename(
          ensembl_protein_id = X1,
          go_id = X3,
          source = X5,
          confidence = X7,
        ) |>
        dplyr::mutate(
          confidence = ceiling(as.numeric(confidence))
        ) |>
        dplyr::select(
          ensembl_protein_id,
          go_id,
          channel,
          confidence
        ) |>
        dplyr::filter(
          confidence > 2
        ) |>
        dplyr::distinct() |>
        dplyr::mutate(source = "HPA")
    }

    if (e == 'textmining') {
      subcell_evidence[[e]] <- subcell_evidence[[e]] |>
        dplyr::mutate(channel = "Text mining", source = "PubMed") |>
        dplyr::rename(
          ensembl_protein_id = X1,
          go_id = X3,
          confidence = X6,
          literature = X7,
        ) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          confidence = min(ceiling(as.numeric(confidence)), 4)
        ) |>
        dplyr::filter(
          stringr::str_detect(
            ensembl_protein_id,"ENSP"
          )
        ) |>
        dplyr::select(
          ensembl_protein_id,
          go_id,
          channel,
          confidence
        ) |>
        dplyr::filter(
          confidence > 2
        ) |>
        dplyr::mutate(source = "PubMed")
    }

    if (e == 'knowledge') {
      subcell_evidence[[e]] <- subcell_evidence[[e]] |>
        dplyr::mutate(channel = "Knowledge") |>
        dplyr::rename(
          ensembl_protein_id = X1,
          go_id = X3,
          confidence = X7,
          source = X5,
          source2 = X6) |>
        dplyr::mutate(source = paste(
          source, source2, sep =" - "
        )) |>
        dplyr::group_by(
          ensembl_protein_id,
          go_id,
          channel,
          confidence,
        ) |>
        dplyr::summarise(
          source = paste(unique(source), collapse=", "),
          .groups = "drop"
        ) |>
        dplyr::arrange(
          ensembl_protein_id, go_id, channel, confidence
        ) |>
        dplyr::group_by(
          ensembl_protein_id, go_id, channel
        ) |>
        dplyr::top_n(-1, confidence) |>
        dplyr::filter(
          confidence > 2
        )
    }

  }

  generic_locs_skip <- dplyr::bind_rows(
    data.frame('go_id' = 'GO:0005575', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0110165', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0005694', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0005622', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0043226', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0043227', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0043229', stringsAsFactors = F),
    data.frame('go_id' = 'GO:0043231', stringsAsFactors = F)
  )

  subcell_loc_compartments <-
    do.call(rbind, subcell_evidence) |>
    dplyr::anti_join(
      generic_locs_skip,  by = "go_id"
    ) |>
    dplyr::left_join(
      dplyr::select(
        go_terms, go_id, go_term
      ),
      by = "go_id",  multiple = "all") |>
    dplyr::left_join(
      ensembl_protein_xref,
      by = "ensembl_protein_id",  multiple = "all"
    ) |>
    dplyr::rename(
      annotation_source = source,
      annotation_channel = channel
    ) |>
    dplyr::mutate(
      entrezgene = as.integer(entrezgene)
    ) |>
    dplyr::select(
      entrezgene,
      go_id,
      go_term,
      annotation_source,
      annotation_channel,
      confidence
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

  compartmentdb <- list()
  compartmentdb[['compartments']] <- subcell_loc_compartments
  compartmentdb[['go_gganatogram_map']] <- go_gganatogram_map
  compartmentdb[['gganatogram_legend']] <- subcell_figure_legend

  saveRDS(compartmentdb, file = rds_fname)


  return(compartmentdb)

}


# get_subcellular_annotations <- function(
#   raw_db_dir = NULL,
#   transcript_xref_db = NULL){
#
#
#
#
#
#
#
#
#
#   uniprot_xref <- transcript_xref_db |>
#     dplyr::filter(property == "uniprot_acc") |>
#     dplyr::select(entrezgene, value) |>
#     dplyr::rename(uniprot_acc = value)
#
#   go_terms <- data.frame()
#   go_structure <- as.list(GO.db::GOTERM)
#   for(n in names(go_structure)){
#     term_df <-
#       data.frame(
#         'go_id' = as.character(n),
#         'go_term' =
#           as.character(AnnotationDbi::Term(go_structure[[n]])),
#         'go_ontology' =
#           AnnotationDbi::Ontology(go_structure[[n]]),
#         stringsAsFactors = F)
#     go_terms <- dplyr::bind_rows(go_terms, term_df)
#   }
#
#   comppi_fname <- file.path(
#     raw_db_dir,
#     "comppi",
#     "comppi_proteins.txt")
#
#   comppi_locations_large_minor_yaml <- file.path(
#     raw_db_dir,
#     "comppi",
#     "largelocs.yml"
#   )
#
#   comppi_locations_large_minor_manual <- file.path(
#     raw_db_dir,
#     "comppi",
#     "largelocs_manual.tsv"
#   )
#
#   gganatogram_map_xlsx <- file.path(
#     raw_db_dir,
#     "comppi",
#     "compartment_mapping_jw.xlsx")
#
#   go_gganatogram_map <-
#     openxlsx::read.xlsx(gganatogram_map_xlsx,
#                         sheet = 2) |>
#     janitor::clean_names() |>
#     dplyr::select(organ, go_id) |>
#     dplyr::mutate(organ = stringr::str_trim(organ)) |>
#     dplyr::left_join(gganatogram::cell_key$cell, by = "organ") |>
#     dplyr::select(-c(type,value,colour)) |>
#     dplyr::rename(ggcompartment = organ)
#
#
#   large_locs_yml <- yaml::read_yaml(file=comppi_locations_large_minor_yaml)
#   largeLoc2MinorLoc <- data.frame()
#   for(loc in names(large_locs_yml$largelocs)){
#     df <- data.frame('major_loc' = loc,
#                      'go_id' = large_locs_yml$largelocs[[loc]]$childrenIncluded,
#                      stringsAsFactors = F)
#     largeLoc2MinorLoc <- largeLoc2MinorLoc |>
#       dplyr::bind_rows(df)
#
#   }
#
#   large_locs_tsv <- as.data.frame(
#     readr::read_tsv(
#       file = comppi_locations_large_minor_manual, show_col_types = F
#     ) |>
#       dplyr::distinct()
#   )
#   largeLoc2MinorLoc <- largeLoc2MinorLoc |>
#     dplyr::bind_rows(large_locs_tsv) |>
#     dplyr::distinct()
#
#
#
#   comppidb <- as.data.frame(
#     read.table(
#       file = comppi_fname,
#       sep = "\t", header = T, quote = "",
#       stringsAsFactors = F,comment.char = "") |>
#       janitor::clean_names() |>
#       dplyr::filter(naming_convention == "UniProtKB/Swiss-Prot/P") |>
#       dplyr::select(
#         major_loc_with_loc_score,
#         protein_name,
#         minor_loc,
#         experimental_system_type,
#         localization_source_database,
#         pubmed_id) |>
#       tidyr::separate_rows(
#         minor_loc, pubmed_id,
#         localization_source_database,
#         experimental_system_type, sep="\\|") |>
#       dplyr::rename(
#         go_id = minor_loc,
#         uniprot_acc = protein_name,
#         annotation_source = localization_source_database) |>
#       dplyr::left_join(go_terms,by=c("go_id")) |>
#       dplyr::filter(!is.na(go_term)) |>
#       tidyr::separate_rows(
#         major_loc_with_loc_score,
#         sep = "\\|"
#       ) |>
#       tidyr::separate(
#         major_loc_with_loc_score,
#         into = c("major_loc","major_loc_score"),
#         sep = ":",
#         remove = T
#       ) |>
#       dplyr::mutate(
#         major_loc_score = as.numeric(major_loc_score)
#       ) |>
#       dplyr::inner_join(
#         largeLoc2MinorLoc, by = c("major_loc","go_id")
#       ) |>
#       dplyr::filter(major_loc != "N/A") |>
#       dplyr::filter(major_loc_score > 0.7) |>
#       dplyr::mutate(
#         annotation_type = dplyr::if_else(
#           stringr::str_detect(experimental_system_type,
#                               "Experimental"),
#           "Experimental","Predicted/Unknown")) |>
#       dplyr::select(-experimental_system_type) |>
#       dplyr::distinct() |>
#       dplyr::left_join(uniprot_xref,by=c("uniprot_acc")) |>
#       dplyr::filter(!is.na(uniprot_acc)) |>
#       dplyr::group_by(uniprot_acc,
#                       go_id,
#                       go_term) |>
#       dplyr::summarise(
#         annotation_source =
#           paste(sort(unique(annotation_source)), collapse="|"),
#         annotation_type =
#           paste(sort(unique(annotation_type)), collapse="|"),
#         .groups = "drop") |>
#       dplyr::mutate(
#         confidence =
#           stringr::str_count(annotation_source,
#                              pattern = "\\|") + 1)
#   )
#
#
#   subcell_figure_legend <- list()
#   cell_key[['cell']]$colour <- "#ffd700"
#   for (i in 1:nrow(cell_key[['cell']])) {
#     subcell_figure_legend[[i]] <-
#       gganatogram::gganatogram(
#         data=cell_key[['cell']][i,],
#         outline = T,
#         fillOutline='#a6bddb',
#         organism="cell",
#         fill="colour")  +
#       ggplot2::theme_void() +
#       ggplot2::ggtitle(cell_key[['cell']][i,]$organ) +
#       ggplot2::theme(
#         plot.title = ggplot2::element_text(
#           hjust=0.5, size=10)
#       ) +
#       ggplot2::coord_fixed()
#   }
#
#   subcelldb <- list()
#   subcelldb[['comppidb']] <- comppidb
#   subcelldb[['go_gganatogram_map']] <- go_gganatogram_map
#   subcelldb[['gganatogram_legend']] <- subcell_figure_legend
#
#   return(subcelldb)
#
# }
