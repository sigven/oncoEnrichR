library(magrittr)
library(TissueEnrich)
library(gganatogram)

msigdb_version <- 'v7.4 (April 2021)'
wikipathways_version <- "20211010"
netpath_version <- "2010"
opentargets_version <- "2021.09"
kegg_version <- "20211006"
gencode_version <- "38"
update_omnipathdb <- F
update_hpa <- F
update_ncbi_gene_summary <- F
update_project_score <- F
update_project_survival <- F
update_tcga <- F
update_cancer_hallmarks <- F
update_omnipath_regulatory <- F
update_omnipath_complexdb <- F
update_gencode <- F
update_ligand_receptor_db <- F


uniprot_release <- "2021_03"

software_db_version <-
  read.table(file="data-raw/RELEASE_NOTES",
             skip = 2, sep = "\t", stringsAsFactors = F,
             comment.char = "#",quote="")
colnames(software_db_version) <-
  c('name','url','description',
    'version','key','resource_type')
release_notes <- list()
i <- 1
while(i <= nrow(software_db_version)){
  release_notes[[software_db_version[i,]$key]] <-
    list('url' = software_db_version[i,]$url,
         'description' =
           software_db_version[i,]$description,
         'version' =
           software_db_version[i,]$version,
         'name' =
           software_db_version[i,]$name,
         'resource_type' =
           software_db_version[i,]$resource_type)
  i <- i + 1
}
rm(software_db_version)
usethis::use_data(release_notes, overwrite = T)


source('data-raw/data_utility_functions.R')

####---MSigDB signatures---####
msigdb <-
  get_msigdb_signatures(basedir = here::here(),
                        db_version = msigdb_version)

####---Gene info----####
gene_info <-
  get_gene_info_ncbi(
    basedir = here::here())

## gene synonyms
alias2primary <-
  get_gene_aliases(gene_info = gene_info)

####---Cancer hallmark annotations----####
cancer_hallmark_data <-
  get_cancer_hallmarks(
    opentargets_version = opentargets_version,
    basedir = here::here(),
    gene_info = gene_info,
    update = update_cancer_hallmarks)

####---dbNSFP: gene summary descriptions---####
uniprot_gene_summary <-
  get_dbnsfp_gene_annotations(
    basedir = here::here())

####----CancerMine/NCG---####
ts_oncogene_annotations <-
  get_ts_oncogene_annotations(
    basedir = here::here(),
    gene_info = gene_info,
    version = "42") %>%
  dplyr::select(
    symbol, tumor_suppressor,
    oncogene, citation_links_oncogene,
    citation_links_tsgene, cancergene_support,
    citation_links_cdriver, cancer_driver)

####--- OTP: Target-disease associations ----####
####--- OTP: Target tractability ----####
opentarget_associations <-
  get_opentarget_associations(
    basedir = here::here(),
    min_num_sources = 2,
    min_overall_score = 0.1,
    direct_associations_only = T)

####---UniProt/CORUM----####
uniprot_map <-
  get_uniprot_map(
    basedir = here::here(),
    uniprot_release = uniprot_release)

####---GENCODE----####

gencode <- get_gencode_data(
  basedir = here::here(),
  uniprot_map = uniprot_map,
  gene_info = gene_info,
  gencode_version = gencode_version,
  update = update_gencode
)

####---OmniPathDB - protein complexes ----####

if(update_omnipath_complexdb == T){

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
                     id2 = paste(id2, collapse=";")) %>%
    dplyr::mutate(num_sources = stringr::str_count(sources,";") + 1) %>%
    dplyr::arrange(complex_name, desc(nchar(uniprot_acc)), desc(num_sources)) %>%
    dplyr::mutate(id2 = dplyr::row_number()) %>%
    dplyr::filter(id2 == 1) %>%
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

  # protein_complex_citations <- get_citations_pubmed(
  #   pmid = unique(pmid2complex$pmid))

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
    readr::read_tsv(file = file.path(here::here(),"data-raw","corum","coreComplexes.txt"),
                    show_col_types = F) %>%
      janitor::clean_names() %>%
    dplyr::filter(organism == "Human") %>%
    dplyr::group_by(complex_name) %>%
    dplyr::summarise(purification_method =
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


  humap_url <-
    "http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt"

  if(!file.exists(file.path(
    here::here(),
    "data-raw","humap2","humap2_complexes.tsv"))){
    download.file(url = humap_url,
                  destfile = file.path(
                    here::here(),
                    "data-raw","humap2","humap2_complexes.tsv"))
  }

  complex_humap <- readr::read_csv(
    file = file.path(here::here(), "data-raw", "humap2", "humap2_complexes.tsv"),
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
          file = file.path(here::here(),"data-raw","omnipathdb","omnipath_complexdb.rds"))

}else{
  proteincomplexdb <- readRDS(
    file = file.path(here::here(),"data-raw","omnipathdb","omnipath_complexdb.rds"))
}




## Gene summaries from NCBI (provided by RefSeq/OMIM):
if(update_ncbi_gene_summary == T){
  i <- 0
  ncbi_gene_summary <- data.frame()
  if(is.null(gencode$entrezgene)){
    rlogging::warning("Column 'entrezgene' is missing from target set")
  }
  pcg <- gencode %>%
    dplyr::filter(gencode_gene_biotype == "protein_coding") %>%
    dplyr::select(entrezgene, symbol) %>%
    dplyr::distinct()

  for(s in unique(pcg$entrezgene)){
    if(is.na(s)){
      next
    }
    summary_json <-
      jsonlite::fromJSON(paste0('http://mygene.info/v3/query?q=',
                                s,'&fields=summary'))
    i <- i + 1
    df <- NULL
    if(!is.null(summary_json$hits$summary)){
      df <- data.frame('entrezgene' = s,
                       'gene_summary_ncbi' = summary_json$hits$summary[1],
                       stringsAsFactors = F)
    }else{
      df <- data.frame('entrezgene' = s,
                       'gene_summary_ncbi' = as.character(NA),
                       stringsAsFactors = F)
    }
    ncbi_gene_summary <-
      dplyr::bind_rows(ncbi_gene_summary, df)
    if(i %% 100 == 0){
      cat("Processed",i,"genes",sep=" " )
      cat("\n")
    }
  }

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

  saveRDS(ncbi_gene_summary,
          file=paste0(here::here(),
                      "/data-raw/ncbi_gene_summary.rds"))

}else{
  ncbi_gene_summary <- readRDS(file="data-raw/ncbi_gene_summary.rds")
}


## UniProt accessions to genes
uniprot_xref <- as.data.frame(
  gencode %>%
    dplyr::select(
      entrezgene, symbol, uniprot_acc,
      uniprot_reviewed) %>%
    dplyr::filter(
      !is.na(uniprot_acc) &
        !is.na(entrezgene) &
        uniprot_reviewed == T) %>%
    tidyr::separate_rows(uniprot_acc, sep="&") %>%
    dplyr::mutate(entrezgene = as.character(entrezgene)) %>%
    dplyr::select(-uniprot_reviewed) %>%
    dplyr::select(symbol, uniprot_acc) %>%
    dplyr::distinct()
)

refseq_mrna_xref <- as.data.frame(
  gencode %>%
    dplyr::select(symbol, refseq_mrna) %>%
    dplyr::filter(!is.na(refseq_mrna) &
                    !is.na(symbol)) %>%
    tidyr::separate_rows(refseq_mrna,sep="&") %>%
    dplyr::distinct()
)

refseq_protein_xref <- as.data.frame(
  gencode %>%
    dplyr::select(symbol, refseq_peptide) %>%
    dplyr::filter(!is.na(refseq_peptide) &
                    !is.na(symbol)) %>%
    tidyr::separate_rows(refseq_peptide,sep="&") %>%
    dplyr::distinct()
)

ensembl_protein_xref <- as.data.frame(
  gencode %>%
    dplyr::select(symbol, ensembl_protein_id) %>%
    dplyr::filter(!is.na(ensembl_protein_id) &
                    !is.na(symbol)) %>%
    tidyr::separate_rows(ensembl_protein_id, sep="&") %>%
    dplyr::mutate(ensembl_protein_id =
                    stringr::str_replace(ensembl_protein_id,
                                         "\\.[0-9]{1,}$","")) %>%
    dplyr::distinct()
)

ensembl_mrna_xref <- as.data.frame(
  gencode %>%
    dplyr::select(symbol, ensembl_transcript_id) %>%
    dplyr::distinct()
)



## Append all gene annotations to a single dataframe
## 1) Remove transcripts
## 2) Add gene name link
## 3) Add disease associations from Open Targets Platform
## 4) Add tumor suppressor/oncogene annotations from CancerMine
## 5) Add gene function summary descriptions from NCBI/UniProt
gene_xref <- gencode %>%
  dplyr::select(-c(ensembl_transcript_id,
                   ensembl_protein_id,
                   refseq_peptide,
                   refseq_mrna)) %>%
  dplyr::distinct() %>%
  resolve_duplicates() %>%
  dplyr::mutate(
    genename =
      paste0("<a href ='http://www.ncbi.nlm.nih.gov/gene/",
             entrezgene,
             "' target='_blank'>",name,"</a>")) %>%
  dplyr::left_join(opentarget_associations,
                   by = c("symbol" = "symbol")) %>%
  dplyr::left_join(ts_oncogene_annotations,
                   by = c("symbol" = "symbol")) %>%
  dplyr::left_join(uniprot_gene_summary,
                   by = c("symbol" = "symbol")) %>%
  dplyr::left_join(ncbi_gene_summary,
                   by = c("entrezgene" = "entrezgene")) %>%
  dplyr::mutate(tumor_suppressor =
                  dplyr::if_else(is.na(tumor_suppressor),
                                 FALSE,
                                 as.logical(tumor_suppressor))) %>%
  dplyr::mutate(oncogene =
                  dplyr::if_else(is.na(oncogene),
                                 FALSE,
                                 as.logical(oncogene))) %>%
  dplyr::mutate(cancer_driver =
                  dplyr::if_else(is.na(cancer_driver),
                                 FALSE,
                                 as.logical(cancer_driver))) %>%
  dplyr::mutate(gene_summary =
                  dplyr::if_else(
                    !is.na(gene_summary_ncbi) &
                      !is.na(gene_summary_uniprot),
                    paste0(gene_summary_ncbi," ",
                           gene_summary_uniprot),
                    as.character(NA)
                  )
  ) %>%
  dplyr::mutate(gene_summary =
                  dplyr::if_else(
                    !is.na(gene_summary_ncbi) &
                      is.na(gene_summary_uniprot),
                    gene_summary_ncbi,
                    as.character(gene_summary)
                  )
  ) %>%
  dplyr::mutate(gene_summary =
                  dplyr::if_else(
                    is.na(gene_summary_ncbi) &
                      !is.na(gene_summary_uniprot),
                    gene_summary_uniprot,
                    as.character(gene_summary)
                  )
  )


corum_complexes <- uniprot_map$corum_complexes

protein_coding_genes <- gene_xref %>%
  dplyr::filter(gencode_gene_biotype == "protein_coding") %>%
  dplyr::select(symbol) %>%
  dplyr::filter(!stringr::str_detect(symbol,"-")) %>%
  dplyr::filter(symbol != "XYLB") %>%
  dplyr::distinct()

####----OmniPathDB - gene annotations----####
omnipathdb <- data.frame()

if(update_omnipathdb == T){
  i <- 1
  omnipathdb <- data.frame()
  while(i <= nrow(protein_coding_genes) - 200){
    genes <-
      protein_coding_genes$symbol[i:min(nrow(protein_coding_genes),i + 199)]
    annotations <-
      OmnipathR::import_omnipath_annotations(proteins = genes) %>%
      dplyr::filter(
        !stringr::str_detect(
          source, "^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$"))
    omnipathdb <-
      dplyr::bind_rows(
        omnipathdb, annotations)
    cat(i,min(nrow(protein_coding_genes),i+199),'\n')
    i <- i + 200
  }
  genes <-
    protein_coding_genes$symbol[i:nrow(protein_coding_genes)]
  annotations <-
    OmnipathR::import_omnipath_annotations(
      select_genes = genes) %>%
    dplyr::filter(!stringr::str_detect(
      source,"^(HPA_tissue|DisGeNet|KEGG|MSigDB|ComPPI|DGIdb|LOCATE|Vesiclepedia|Ramilowski2015)$"))
  omnipathdb <-
    dplyr::bind_rows(omnipathdb, annotations)

  saveRDS(omnipathdb,
          file=paste0(here::here(),
                      "/data-raw/omnipathdb/omnipathdb.rds"))
}else{
  omnipathdb <- readRDS(
    file = paste0(here::here(),
                  "/data-raw/omnipathdb/omnipathdb.rds"))
}



####---OmniPathDB - TF interactions ----####

if(update_omnipath_regulatory == T){
  tf_target_interactions <- get_tf_target_interactions()

  saveRDS(tf_target_interactions,
          file = file.path(here::here(), "data-raw",
                           "omnipathdb", "omnipath_tf_interactions.rds"))

}else{
  tf_target_interactions <-
    readRDS(file.path(here::here(), "data-raw",
                      "omnipathdb", "omnipath_tf_interactions.rds"))
}






####--- UMLS/EFO/DiseaseOntology----####

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


####---- OTP - Target-disease associations---####
otdb_tmp <- as.data.frame(
  #dplyr::select(opentarget_associations2, ot_association, symbol) %>%
  dplyr::select(gene_xref,
                ot_association,
                symbol,
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

set1 <- otdb_tmp %>% dplyr::filter(!is.na(primary_site)) %>% dplyr::distinct()
set2 <- otdb_tmp %>% dplyr::filter(is.na(primary_site) & cancer_phenotype == T) %>%
  dplyr::select(-c(efo_name, primary_site)) %>%
  dplyr::anti_join(
    dplyr::select(set1, disease_efo_id, symbol),
    by = c("disease_efo_id","symbol")) %>%
  dplyr::inner_join(
    dplyr::select(oncoPhenoMap::auxiliary_maps$efo$efo2name,
                  efo_id, efo_name, primary_site),
    by = c("disease_efo_id" = "efo_id")) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(primary_site))

set12 <- dplyr::bind_rows(set1, set2)

set3 <- otdb_tmp %>%
  dplyr::anti_join(dplyr::select(set12, disease_efo_id, symbol),
                   by = c("disease_efo_id","symbol")) %>%
  dplyr::distinct()

otdb_all <- set12 %>%
  dplyr::bind_rows(set3) %>%
  dplyr::arrange(primary_site, symbol,
                 desc(ot_association_score))

rm(otdb_tmp)
rm(set12)
rm(set3)
rm(set1)
rm(set2)

otdb_tissue_scores <- as.data.frame(
  otdb_all %>%
    dplyr::filter(!is.na(primary_site)) %>%
    dplyr::group_by(primary_site, symbol) %>%
    dplyr::summarise(
      tissue_assoc_score =
        round(sum(ot_association_score), digits = 3)
    ) %>%
    dplyr::arrange(primary_site, desc(tissue_assoc_score))
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

otdb <- list()
otdb$all <- otdb_all
otdb$site_rank <- otdb_site_rank

gene_xref <- gene_xref %>%
  dplyr::select(-c(ot_association))


####--- Ligand-Receptor interactions ----####

if(update_ligand_receptor_db == T){

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
          file = file.path(here::here(),"data-raw",
                           "cellchatdb","cellchatdb_interactions.rds"))
  rm(ligand_receptor_citations)
  rm(ligand_receptor_db_final)
  rm(ligand_receptor_db)
  rm(ligand_receptor_kegg_support)
  rm(ligand_receptor_literature)
  rm(interaction2ligands)
  rm(interaction2receptor)

}else{
  ligandreceptordb <-
    readRDS(file.path(here::here(),"data-raw","cellchatdb","cellchatdb_interactions.rds"))
}



####----Gene Ontology----####
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

####----ComPPI - subcellular compartments---####
go_gganatogram_map <-
  openxlsx::read.xlsx("data-raw/compartment_mapping_jw.xlsx",
                      sheet = 2) %>%
  janitor::clean_names() %>%
  dplyr::select(organ, go_id) %>%
  dplyr::mutate(organ = stringr::str_trim(organ)) %>%
  dplyr::left_join(gganatogram::cell_key$cell, by = "organ") %>%
  dplyr::select(-c(type,value,colour)) %>%
  dplyr::rename(ggcompartment = organ)

comppidb <- as.data.frame(
  read.table(
    file=paste0(here::here(),
                '/data-raw/compppi/comppi_proteins.txt'),
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
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::group_by(symbol,
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


if(update_project_score == T){
  ####---Project Score - CRISPR/Cas9---####

  cell_lines <- read.table(file="data-raw/project_score/model_list.csv",quote="",
                           comment.char="",sep=",", header = F, stringsAsFactors = F,
                           fill = T,na.strings = c("")) %>%
    dplyr::select(V1, V2, V3, V4, V18, V22, V24, V25, V32, V40, V39, V41) %>%
    dplyr::filter(V4 == "Cell Line" & V39 == "Homo Sapiens") %>%
    dplyr::filter(!stringr::str_detect(V1,"model_id")) %>%
    dplyr::mutate(V2 = stringr::str_replace_all(V2,"-",".")) %>%
    dplyr::filter(V18 == "True") %>%
    magrittr::set_colnames(c('model_id','model_name','synonyms','model_type',
                             'crispr_ko_data','tissue','cancer_type',
                            'cancer_type_detail','sample_site',
                            'gender','species','ethnicity')) %>%
    dplyr::select(-c(gender,species,ethnicity,cancer_type_detail,cancer_type,
                     crispr_ko_data,synonyms)) %>%
    dplyr::filter(tissue != "Eye" & tissue != "Biliary Tract" &
                    tissue != "Prostate" & tissue != "Soft Tissue") %>%
    dplyr::mutate(tissue = dplyr::if_else(tissue == "Central Nervous System",
                                          "CNS/Brain",
                                          as.character(tissue))) %>%
    dplyr::mutate(tissue = dplyr::if_else(tissue == "Large Intestine",
                                          "Colon/Rectum",
                                          as.character(tissue)))

  gene_identifiers <- read.csv(file="data-raw/project_score/gene_identifiers.csv",
                               stringsAsFactors = F) %>%
    dplyr::rename(gene_id_project_score = gene_id) %>%
    dplyr::rename(entrezgene = entrez_id) %>%
    dplyr::mutate(entrezgene = as.character(entrezgene)) %>%
    dplyr::filter(!is.na(entrezgene)) %>%
    dplyr::left_join(
      dplyr::select(gene_xref, symbol,
                    entrezgene), by = "entrezgene") %>%
    dplyr::select(-c(cosmic_gene_symbol, hgnc_id, hgnc_symbol,
                     refseq_id, uniprot_id, ensembl_gene_id))

  bayes_factors <-
    read.table(file="data-raw/project_score/03_scaledBayesianFactors.tsv",header=T,
               quote="",sep="\t", stringsAsFactors = F)

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
    dplyr::left_join(gene_identifiers,by=c("symbol")) %>%
    dplyr::filter(!is.na(gene_id_project_score))


  ## Fitness scores
  cell_fitness_scores <-
    read.table(file="data-raw/project_score/binaryDepScores.tsv",header=T,
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
    read.csv(file="data-raw/project_score/depmap-priority-scores.csv",
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
    dplyr::mutate(symbol = factor(symbol, levels = priority_order$symbol))

  usethis::use_data(projectscoredb, overwrite = T)

  rm(cell_fitness_scores)
  rm(cell_lines)
  rm(priority_order)
  rm(gene_identifiers)
}


####---OTP - Targeted cancer drugs ---####
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
  dplyr::filter(!is.na(cancer_drug) & cancer_drug == T) %>%
  dplyr::select(target_symbol,
                target_entrezgene,
                drug_name,
                primary_site,
                molecule_chembl_id,
                drug_max_ct_phase,
                cancer_drug) %>%
  dplyr::distinct()


cancer_drugs[['early_phase']] <-
  oncoPharmaDB::get_onco_drugs(
    drug_is_targeted = T,
    source_opentargets_only = T,
    drug_minimum_phase_any_indication = 0) %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::anti_join(
    cancer_drugs[['late_phase']], by = "molecule_chembl_id") %>%
  dplyr::filter(!is.na(cancer_drug) & cancer_drug == T) %>%
  #dplyr::filter(!is.na(primary_site))
  dplyr::select(target_symbol,
                target_entrezgene,
                drug_name,
                primary_site,
                molecule_chembl_id,
                drug_max_ct_phase,
                cancer_drug) %>%
  dplyr::distinct()

for(p in c('early_phase','late_phase')){

  cancer_drugs[[p]] <- cancer_drugs[[p]] %>%
    dplyr::mutate(to = paste0("s",molecule_chembl_id)) %>%
    dplyr::mutate(from = paste0("s",target_entrezgene)) %>%
    dplyr::distinct()

  drugs_per_gene[[p]] <- as.data.frame(
    cancer_drugs[[p]] %>%
      dplyr::select(target_symbol,
                    drug_name,
                    drug_max_ct_phase,
                    molecule_chembl_id) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        drug_link =
          paste0("<a href = 'https://platform.opentargets.org/drug/",
                 molecule_chembl_id,"' target='_blank'>",drug_name,"</a>")) %>%
      dplyr::arrange(target_symbol, desc(drug_max_ct_phase)) %>%
      dplyr::group_by(target_symbol) %>%
      dplyr::summarise(
        targeted_cancer_drugs_lp =
          paste(unique(drug_link),
                collapse = ", "), .groups = "drop")
  )
  if(p == 'early_phase'){
    drugs_per_gene[[p]] <- drugs_per_gene[[p]] %>%
      dplyr::rename(
        targeted_cancer_drugs_ep = targeted_cancer_drugs_lp)
  }

  gene_xref <- dplyr::left_join(
    gene_xref,
    drugs_per_gene[[p]],
    by = c("symbol" = "target_symbol"))

  target_drug_edges[[p]] <- cancer_drugs[[p]] %>%
    dplyr::select(drug_name,
                  from,
                  to,
                  target_symbol,
                  target_entrezgene,
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

####--- Cancer-KM-Survival (CSHL) ---####
if(update_project_survival == T){

  for(feature_type in
      c('CNAs','Mutations','Gene expression','Methylation',
        'miRNA expression','Protein expression')){
    destfile_fname <-
      paste0("data-raw/km_survival_cshl/",
             tolower(
               stringr::str_replace(
                 stringr::str_replace(feature_type," ","_"),
                 "s$","")),
             ".xlsx")
    if(!file.exists(destfile_fname)){
      download.file(url = paste0(
        "https://storage.googleapis.com/cancer-survival-km-2021/cancer-survival-km/full-downloads/",
        stringr::str_replace(feature_type," ","%20"),".xlsx"),
        destfile = destfile_fname,
        quiet = T)
    }

  }

  project_survival <- list()
  for(feature_type in c('cna','mutation','gene_expression','methylation')){

    fname <- file.path("data-raw","km_survival_cshl",
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

  usethis::use_data(projectsurvivaldb, overwrite = T)

  rm(pancancer_trend)
  rm(project_survival_df)
  rm(projectsurvivaldb)
  rm(data)

}




####--- Human Protein Atlas ---####

if(update_hpa == T){
  hpa <- data.frame()

  protein_coding_genes <- gene_xref %>%
    dplyr::filter(gencode_gene_biotype == "protein_coding" &
                    !is.na(ensembl_gene_id))

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

  for (i in 1:nrow(protein_coding_genes)) {
    g <- protein_coding_genes[i,]

    local_json <-
      paste0("data-raw/hpa/json/",
             g$ensembl_gene_id,".json")
    if(!file.exists(local_json)){
      protein_atlas_url <-
        paste0("https://www.proteinatlas.org/",
               g$ensembl_gene_id,".json")
      if(RCurl::url.exists(protein_atlas_url)){
        download.file(
          protein_atlas_url,
          destfile=paste0("data-raw/hpa/json/",
                          g$ensembl_gene_id, ".json"),
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
            df <- data.frame('ensembl_gene_id' = g$ensembl_gene_id,
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

    if(i %% 100 == 0){
      cat('Completed ',i,' genes..','\n')
      saveRDS(hpa, file=paste0(here::here(),"/data-raw/hpa/hpa.rds"))
    }

  }
  saveRDS(hpa,
          file=paste0(here::here(),"/data-raw/hpa/hpa.rds"))
}else{
  hpa <- readRDS(
    file = paste0(here::here(),"/data-raw/hpa/hpa.rds"))
}

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
  if(!file.exists(
    paste0(
      here::here(),
      "/data-raw/hpa/",
      t,
      ".tsv"))){
    download.file(
      url = paste0("https://www.proteinatlas.org/download/",
                   t, ".tsv.zip"),
      destfile = paste0(here::here(),
                        "/data-raw/hpa/",
                        t,
                        ".tsv.zip")
    )
    system(paste0('unzip ',
                  here::here(),
                  "/data-raw/hpa/",
                  t,
                  ".tsv.zip",
                  " -d ",
                  here::here(),
                  "/data-raw/hpa"
                  )
           )
  }

  ##set max number of tissues/cell_types for
  ##determination of groups in group-enriched genes
  max_types_in_group <- 10
  if(t == "rna_tissue_gtex"){
    max_types_in_group <- 5
  }

  data <- read.table(file = paste0(here::here(),
                                   "/data-raw/hpa/",
                                   t,
                                   ".tsv"),
                     header = T, stringsAsFactors = F,
                     sep = "\t")
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
    tissue_cell_expr[['tissue']][['unit']] <- 'NX'
    tissue_cell_expr[['tissue']][['te_df']] <-
      as.data.frame(
        as.data.frame(
          assay(tissue_cell_expr[['tissue']][['te_SE']])
        ) %>%
          dplyr::rename(ensembl_gene_id = Gene,
                        category = Group) %>%
          dplyr::left_join(
            dplyr::select(gene_xref, symbol,
                          name,
                          entrezgene,
                          ensembl_gene_id),
            by = "ensembl_gene_id") %>%
          dplyr::rename(genename = name) %>%
          dplyr::filter(!is.na(entrezgene)) %>%
          dplyr::mutate(
            genename = paste0(
              "<a href='https://gtexportal.org/home/gene/",
              ensembl_gene_id,"' target='_blank'>",
              genename,"</a>")
          ) %>%
          dplyr::group_by(ensembl_gene_id,
                          symbol,
                          genename,
                          entrezgene,
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
    tissue_cell_expr[['single_cell']][['unit']] <- 'TPM'
    tissue_cell_expr[['single_cell']][['te_SE']] <-
      te_gene_retrieval_se
    tissue_cell_expr[['single_cell']][['te_df']] <-
      as.data.frame(
        as.data.frame(
          assay(tissue_cell_expr[['single_cell']][['te_SE']])
        ) %>%
          dplyr::rename(ensembl_gene_id = Gene,
                        category = Group) %>%
          dplyr::left_join(
            dplyr::select(gene_xref, symbol,
                          name,
                          ensembl_gene_id,
                          entrezgene),
            by = "ensembl_gene_id") %>%
          dplyr::rename(genename = name) %>%
          dplyr::filter(!is.na(entrezgene)) %>%
          dplyr::mutate(
            genename = paste0(
              "<a href='https://www.proteinatlas.org/",
              ensembl_gene_id,
              "-",
              symbol,
              "/celltype' target='_blank'>",
              genename,"</a>")
          ) %>%
          dplyr::group_by(ensembl_gene_id,
                          symbol,
                          genename,
                          entrezgene,
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

cancerdrugdb <- list()
cancerdrugdb[['tractability']] <- list()
cancerdrugdb[['tractability']][['ab']] <-
  opentarget_associations %>%
  dplyr::select(target_tractability_symbol,
                target_ensembl_gene_id,
                AB_tractability_category,
                AB_tractability_support) %>%
  dplyr::rename(symbol = target_tractability_symbol,
                ensembl_gene_id = target_ensembl_gene_id)
cancerdrugdb[['tractability']][['sm']] <-
  opentarget_associations %>%
  dplyr::select(target_tractability_symbol,
                target_ensembl_gene_id,
                SM_tractability_category,
                SM_tractability_support) %>%
  dplyr::rename(symbol = target_tractability_symbol,
                ensembl_gene_id = target_ensembl_gene_id)
cancerdrugdb[['ppi']] <- list()
cancerdrugdb[['ppi']][['edges']] <-
  dplyr::bind_rows(target_drug_edges[['early_phase']],
                   target_drug_edges[['late_phase']]) %>%
  dplyr::select(from,to,width) %>%
  dplyr::distinct()

cancerdrugdb[['ppi']][['nodes']] <-
  dplyr::bind_rows(target_drug_nodes[['early_phase']],
                   target_drug_nodes[['late_phase']]) %>%
  dplyr::distinct()


####----CORUM protein complexes----####
corumdb <- readRDS(file="data-raw/corum/corumdb.rds")

####----TCGA aberration data ----####

if(update_tcga == T){
  tcga_clinical <-
    readRDS(file="data-raw/tcga/tcga_clinical.rds")
  tcga_aberration_stats <-
    readRDS(file="data-raw/tcga/tcga_gene_aberration_stats.rds") %>%
    dplyr::filter(
      clinical_strata == "site" |
        (clinical_strata == "site_diagnosis" &
           percent_mutated >= 1))

  tcga_aberration_stats$genomic_strata <- NULL
  tcga_aberration_stats$entrezgene <- NULL
  tcga_aberration_stats$consensus_calls <- NULL
  tcga_aberration_stats$fp_driver_gene <- NULL

  maf_codes <- read.table(file="data-raw/maf_codes.tsv",
                          header = T, sep = "\t", quote ="")
  maf_path <- 'data-raw/tcga/maf_path'
  i <- 1
  while(i <= nrow(maf_codes)){
    primary_site <- maf_codes[i,]$primary_site
    maf_code <- maf_codes[i,]$code
    maf_file <- paste0(maf_path,"/tcga_mutation_grch38_release30_20210923.",maf_code,"_0.maf.gz")
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

      write.table(tmp, file="tmp.maf", quote=F, row.names = F, col.names = T, sep ="\t")
      system('gzip -f tmp.maf', intern = T)
      maf <- maftools::read.maf("tmp.maf.gz", verbose = F, clinicalData = clinical)
      saveRDS(maf, paste0("maf/",maf_code,".maf.rds"))
    }
    i <- i + 1
    cat(primary_site,'\n')

  }
  system('rm -f tmp.maf.gz')


  pfam_domains <- as.data.frame(
    readr::read_tsv(
      file="data-raw/pfam/pfam.domains.tsv.gz",
      show_col_types = F
    )) %>%
    dplyr::rename(PFAM_ID = pfam_id, PROTEIN_DOMAIN = url) %>%
    dplyr::rename(PFAM_DOMAIN_NAME = name) %>%
    dplyr::select(PFAM_ID, PFAM_DOMAIN_NAME) %>%
    dplyr::mutate(PFAM_ID = stringr::str_replace(
      PFAM_ID,"\\.[0-9]{1,}$","")
    )

  ####----TCGA recurrent SNVs/InDels-----####
  recurrent_tcga_variants <- as.data.frame(readr::read_tsv(
    file="data-raw/tcga/tcga_recurrent_coding_gvanno_grch38.tsv.gz",
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

  ####----TCGA co-expression data----####
  raw_coexpression <-
    readr::read_tsv("data-raw/tcga/co_expression_strong_moderate.release30_20210923.tsv.gz",
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

  tcga_coexp_db <- dplyr::bind_rows(co_expression_genes1, co_expression_genes2) %>%
    dplyr::arrange(desc(r)) %>%
    dplyr::mutate(p_value = signif(p_value, digits = 4)) %>%
    dplyr::filter(p_value < 1e-6) %>%
    dplyr::mutate(r = signif(r, digits = 3)) %>%
    dplyr::filter((r <= -0.6 & r < 0) | (r >= 0.6))


  tcgadb <- list()
  tcgadb[['co_expression']] <- tcga_coexp_db
  tcgadb[['aberration']] <- tcga_aberration_stats
  tcgadb[['recurrent_variants']] <- recurrent_tcga_variants

  usethis::use_data(tcgadb, overwrite = T)
  usethis::use_data(maf_codes, overwrite = T)
  usethis::use_data(pfam_domains, overwrite = T)

  rm(co_expression_genes1)
  rm(co_expression_genes2)
  rm(raw_coexpression)
  rm(tcga_clinical)
  rm(pfam_domains)
  rm(maf_codes)
  rm(recurrent_tcga_variants)
  rm(tcga_coexp_db)
  rm(tcga_aberration_stats)
  rm(tcgadb)
  rm(maf_path)
  rm(maf_file)

}



####----Poorly defined genes ----####
go_annotations_qc <-
  read.table(gzfile("data-raw/gene_ontology/goa_human.gaf.gz"),
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
    go_ids_function =
      paste(unique(sort(go_id)),
            collapse = "|"),
    num_ids_function = length(unique(go_id)),
    .groups = "drop"
    )

go_process_terms_prgene <- go_annotations_qc %>%
  dplyr::filter(go_ontology == "P") %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    go_ids_process =
      paste(unique(sort(go_id)),
            collapse = "|"),
    num_ids_process = length(unique(go_id)),
    .groups = "drop"
  )

genes_function_knowledge_level <- gene_xref %>%
  dplyr::filter(gencode_gene_biotype == "protein_coding") %>%
  dplyr::select(symbol, ensembl_gene_id,
                name, entrezgene,
                gene_summary) %>%
  dplyr::left_join(go_process_terms_prgene, by = "symbol") %>%
  dplyr::left_join(go_function_terms_prgene, by = "symbol") %>%
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

poorly_defined_genes <- genes_function_knowledge_level %>%
  dplyr::filter(!is.na(symbol) & !stringr::str_detect(symbol,"-")) %>%
  dplyr::filter(is.na(name) | !stringr::str_detect(name,"readthrough")) %>%
  dplyr::filter(!is.na(unknown_function_rank)) %>%
  dplyr::rename(genename = name) %>%
  dplyr::select(symbol, genename, num_go_terms, entrezgene,
                unknown_function_rank, gene_summary) %>%
  dplyr::mutate(has_gene_summary = dplyr::if_else(
    is.na(gene_summary),FALSE,TRUE,TRUE)
  ) %>%
  dplyr::arrange(unknown_function_rank, symbol) %>%
  dplyr::mutate(genename = paste0(
    '<a href="https://www.ncbi.nlm.nih.gov/gene/',
    entrezgene,'" target="_blank">',genename,'</a>')
  ) %>%
  dplyr::select(-entrezgene)



gene_xref$gene_summary_uniprot <- NULL
gene_xref$gene_summary_ncbi <- NULL

msigdb$df <- NULL
tmp <- msigdb$db
rm(msigdb)
msigdb <- tmp
rm(tmp)

pathwaydb <- list()
pathwaydb[['wikipathway']] <- wikipathwaydb
pathwaydb[['netpath']] <- netpathdb
pathwaydb[['kegg']] <- keggdb
pathwaydb[['msigdb']] <- msigdb
rm(wikipathwaydb)
rm(netpathdb)
rm(keggdb)

genedb <- list()
genedb[['all']] <- gene_xref
genedb[['alias2primary']] <- alias2primary
genedb[['poorly_defined']] <- poorly_defined_genes
genedb[['uniprot_xref']] <- uniprot_xref
genedb[['refseq_mrna_xref']] <- refseq_mrna_xref
genedb[['ensembl_mrna_xref']] <- ensembl_mrna_xref
genedb[['refseq_protein_xref']] <- refseq_protein_xref
genedb[['ensembl_protein_xref']] <- ensembl_protein_xref
genedb[['protein_complex']] <- proteincomplexdb
genedb[['cancer_hallmark']] <- cancer_hallmark_data


usethis::use_data(pathwaydb, overwrite = T)
usethis::use_data(ligandreceptordb, overwrite = T)
usethis::use_data(genedb, overwrite = T)
usethis::use_data(otdb, overwrite = T)
usethis::use_data(tissue_cell_expr, overwrite = T)
usethis::use_data(hpa, overwrite = T)
usethis::use_data(release_notes, overwrite = T)
usethis::use_data(subcelldb, overwrite = T)
usethis::use_data(cancerdrugdb, overwrite = T)
usethis::use_data(tf_target_interactions, overwrite = T)

rm(gene_xref)
rm(gencode)
rm(te_gene_retrieval_se)
rm(se)
rm(msigdb)
rm(hpa)
rm(omnipathdb)
rm(wikipathways)
rm(protein_coding_genes)
rm(genedb)
rm(kegg_pathway_genes)
rm(otdb)
rm(corumdb)
rm(pathwaydb)
rm(release_notes)
rm(cancerdrugdb)
rm(go_structure)
rm(uniprot_xref)
rm(target_drug_edges)
rm(target_drug_nodes)
rm(comppidb)
rm(cancer_drugs)
rm(wp2gene)
rm(go_terms)
rm(term_df)
rm(drugs_per_gene)
rm(gene_info)

rm(corum_complexes)
rm(alias2primary)
rm(genes_function_knowledge_level)
rm(ncbi_gene_summary)
rm(uniprot_gene_summary)
rm(tissue_cell_expr)
rm(opentarget_associations)
rm(ts_oncogene_annotations)
rm(uniprot_map)
rm(go_annotations_qc)
rm(go_function_terms_prgene)
rm(go_process_terms_prgene)
rm(poorly_defined_genes)
rm(subcell_figure_legend)
rm(cell_key)
rm(subcelldb)
rm(t)
rm(n)
rm(p)
rm(phenotype_cancer_efo)
rm(go_gganatogram_map)
rm(cancer_hallmark_data)
rm(tf_target_interactions)
rm(ligandreceptordb)
rm(proteincomplexdb)
rm(otdb_all)
rm(otdb_site_rank)
rm(otdb_tissue_scores)
rm(ensembl_mrna_xref)
rm(refseq_mrna_xref)
rm(refseq_protein_xref)
rm(ensembl_protein_xref)
rm(netpath_idmapping)
rm(netpath_pathway_data)

