library(magrittr)
library(TissueEnrich)
library(gganatogram)

msigdb_version <- 'v7.5.1 (Jan 2022)'
wikipathways_version <- "20220210"
netpath_version <- "2010"
opentargets_version <- "2021.11"
kegg_version <- "20211223"
gencode_version <- "39"
update_omnipathdb <- F
update_hpa <- F
update_ncbi_gene_summary <- F
update_project_score <- F
update_project_survival <- F
update_tcga <- F
update_cancer_hallmarks <- F
update_omnipath_regulatory <- F
update_omnipath_complexdb <- F
update_gencode <- T
update_ligand_receptor_db <- F


uniprot_release <- "2021_04"

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

####---Gene info----####
gene_info <-
  get_gene_info_ncbi(
    basedir = here::here())

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
    entrezgene, tumor_suppressor,
    oncogene, citation_links_oncogene,
    citation_links_tsgene, cancergene_support,
    citation_links_cdriver, cancer_driver)

####--- OTP: Target-disease associations ----####
opentarget_associations <-
  get_opentarget_associations(
    basedir = here::here(),
    min_num_sources = 2,
    min_overall_score = 0.1,
    direct_associations_only = T)

otdb <- quantify_gene_cancer_relevance(
  opentarget_associations)

####---GENCODE----####

gencode <- list()
for(build in c('grch37','grch38')){

  gencode_v <- gencode_version
  if(build == 'grch37'){
    gencode_v <- "19"
  }

  gencode[[build]] <- get_gencode_data(
    basedir = here::here(),
    gene_info = gene_info,
    build = build,
    gencode_version = gencode_v,
    update = update_gencode
  )

}

transcript_xref_db <- get_unique_transcript_xrefs(
  basedir = here::here(),
  gene_info = gene_info,
  gencode = gencode,
  update = T
)


####---Pathway annotations ----####
pathwaydb <- get_pathway_annotations(
  basedir = here::here(),
  gene_info = gene_info,
  wikipathways_version = wikipathways_version,
  kegg_version = kegg_version,
  netpath_version = netpath_version,
  msigdb_version = msigdb_version,
  msigdb = T,
  netpath = T,
  wikipathways = T,
  kegg = T
)

####---OmniPathDB - protein complexes ----####
proteincomplexdb <- get_protein_complexes(
  basedir = here::here(),
  update = update_omnipath_complexdb
)

## Gene summaries from NCBI (provided by RefSeq/OMIM):
ncbi_gene_summary <- get_function_summary_ncbi(
  basedir = here::here(),
  gene_df = gene_info,
  update = update_ncbi_gene_summary)

####---OTP - Targeted cancer drugs ---####
cancerdrugdb <- get_cancer_drugs()

####--- Gene Ontology ---####
go_terms_pr_gene <- get_gene_go_terms(
  basedir = here::here()
)

## Append all gene annotations to a single dataframe
## 1) Remove transcripts
## 2) Add gene name link
## 3) Add targeted drugs (early/late phase) and tractability category
## 4) Add tumor suppressor/oncogene annotations from CancerMine
## 5) Add gene function summary descriptions from NCBI/UniProt
## 6) Assign unknown function rank

gene_xref <- generate_gene_xref_df(
  gene_info = gene_info,
  transcript_xref_db = transcript_xref_db,
  ts_oncogene_annotations = ts_oncogene_annotations,
  go_terms_pr_gene = go_terms_pr_gene,
  ncbi_gene_summary = ncbi_gene_summary,
  uniprot_gene_summary = uniprot_gene_summary,
  otdb = otdb,
  cancerdrugdb = cancerdrugdb
)

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
  #resolve_duplicates() %>%
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

####---OmniPathDB - gene annotations ----####
omnipathdb <- get_omnipath_gene_annotations(
  basedir = here::here(),
  gene_info = gene_info,
  update = update_omnipathdb,
)

####---OmniPathDB - TF interactions ----####
tf_target_interactions <- get_tf_target_interactions(
  basedir = here::here(),
  update = update_omnipath_regulatory
)


####--- Ligand-Receptor interactions ----####
ligandreceptordb <- get_ligand_receptors(
  basedir = NULL,
  keggdb = pathwaydb$kegg,
  update = update_ligand_receptor_db
)

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
subcelldb <- get_subcellular_annotations(
  basedir = here::here()
)

####---- Project Score/CRISPR ----####
projectscoredb <- get_crispr_scores(
  basedir = here::here(),
  gene_info = gene_info
)

####--- Cancer-KM-Survival (CSHL) ---####
projectsurvivaldb <- get_survival_associations(
  gene_info = gene_info,
  basedir = here::here()
)

####--- Human Protein Atlas ---####
hpa <- get_hpa_associations(
  basedir = here::here(),
  gene_xref = gene_xref,
  update = update_hpa
)

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

  tcga_coexp_db <- dplyr::bind_rows(co_expression_genes1,
                                    co_expression_genes2) %>%
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
# go_annotations_qc <-
#   read.table(gzfile("data-raw/gene_ontology/goa_human.gaf.gz"),
#              sep = "\t", comment.char = "!", header = F,
#              stringsAsFactors = F, quote = "") %>%
#   dplyr::select(V3, V5, V7, V9, V14) %>%
#   dplyr::rename(symbol = V3, go_id = V5,
#                 go_evidence_code = V7,
#                 go_ontology = V9, go_annotation_date = V14) %>%
#   dplyr::distinct() %>%
#   dplyr::filter(go_evidence_code != "IEA" &
#                   go_evidence_code != "ND") %>%
#   dplyr::filter(go_ontology != "C")
#
# go_function_terms_prgene <- go_annotations_qc %>%
#   dplyr::filter(go_ontology == "F") %>%
#   dplyr::group_by(symbol) %>%
#   dplyr::summarise(
#     go_ids_function =
#       paste(unique(sort(go_id)),
#             collapse = "|"),
#     num_ids_function = length(unique(go_id)),
#     .groups = "drop"
#     )
#
# go_process_terms_prgene <- go_annotations_qc %>%
#   dplyr::filter(go_ontology == "P") %>%
#   dplyr::group_by(symbol) %>%
#   dplyr::summarise(
#     go_ids_process =
#       paste(unique(sort(go_id)),
#             collapse = "|"),
#     num_ids_process = length(unique(go_id)),
#     .groups = "drop"
#   )
#
# genes_function_knowledge_level <- gene_xref %>%
#   dplyr::filter(gencode_gene_biotype == "protein_coding") %>%
#   dplyr::select(symbol, ensembl_gene_id,
#                 name, entrezgene,
#                 gene_summary) %>%
#   dplyr::left_join(go_process_terms_prgene, by = "symbol") %>%
#   dplyr::left_join(go_function_terms_prgene, by = "symbol") %>%
#   dplyr::mutate(num_ids_process =
#                   dplyr::if_else(
#                     is.na(num_ids_process),as.integer(0),
#                     as.integer(num_ids_process)
#                   )
#   ) %>%
#   dplyr::mutate(num_ids_function =
#                   dplyr::if_else(
#                     is.na(num_ids_function),as.integer(0),
#                     as.integer(num_ids_function)
#                   )
#   ) %>%
#   dplyr::mutate(num_go_terms = num_ids_function +
#                   num_ids_process) %>%
#   dplyr::mutate(
#     unknown_function_rank =
#       dplyr::case_when(
#         stringr::str_detect(tolower(name),
#                             "uncharacterized|open reading frame") &
#           is.na(gene_summary) &
#           num_go_terms == 0 ~ as.integer(1),
#         stringr::str_detect(tolower(name),
#                             "uncharacterized|open reading frame") &
#           is.na(gene_summary) &
#           num_go_terms > 0 &
#           num_go_terms <= 3 ~ as.integer(2),
#         is.na(gene_summary) &
#           num_go_terms == 0 &
#           !stringr::str_detect(tolower(name),
#                                "uncharacterized|open reading frame") ~ as.integer(3),
#         is.na(gene_summary) &
#           num_go_terms == 1 &
#           !stringr::str_detect(tolower(name),
#                                "uncharacterized|open reading frame") ~ as.integer(4),
#         !is.na(gene_summary) &
#           num_go_terms == 0 &
#           !stringr::str_detect(tolower(name),
#                                "uncharacterized|open reading frame") ~ as.integer(5),
#         is.na(gene_summary) &
#           num_go_terms == 2 &
#           !stringr::str_detect(tolower(name),
#                                "uncharacterized|open reading frame") ~ as.integer(6),
#         TRUE ~ as.integer(NA)
#       )
#   )
#
# poorly_defined_genes <- genes_function_knowledge_level %>%
#   dplyr::filter(!is.na(symbol) & !stringr::str_detect(symbol,"-")) %>%
#   dplyr::filter(is.na(name) | !stringr::str_detect(name,"readthrough")) %>%
#   dplyr::filter(!is.na(unknown_function_rank)) %>%
#   dplyr::rename(genename = name) %>%
#   dplyr::select(symbol, genename, num_go_terms, entrezgene,
#                 unknown_function_rank, gene_summary) %>%
#   dplyr::mutate(has_gene_summary = dplyr::if_else(
#     is.na(gene_summary),FALSE,TRUE,TRUE)
#   ) %>%
#   dplyr::arrange(unknown_function_rank, symbol) %>%
#   dplyr::mutate(genename = paste0(
#     '<a href="https://www.ncbi.nlm.nih.gov/gene/',
#     entrezgene,'" target="_blank">',genename,'</a>')
#   ) %>%
#   dplyr::select(-entrezgene)
#
#
#
# gene_xref$gene_summary_uniprot <- NULL
# gene_xref$gene_summary_ncbi <- NULL
#
# msigdb$df <- NULL
# tmp <- msigdb$db
# rm(msigdb)
# msigdb <- tmp
# rm(tmp)
#
# pathwaydb <- list()
# pathwaydb[['wikipathway']] <- wikipathwaydb
# pathwaydb[['netpath']] <- netpathdb
# pathwaydb[['kegg']] <- keggdb
# pathwaydb[['msigdb']] <- msigdb
# rm(wikipathwaydb)
# rm(netpathdb)
# rm(keggdb)

for(n in c('cancerdrugdb',
           'genedb',
           'hpa',
           'ligandreceptordb',
           'otdb',
           'pathwaydb',
           'pfamdomaindb',
           'projectscoredb',
           'projectsurvivaldb',
           'subcelldb',
           'tcgadb',
           'tftargetdb',
           'tissuecelldb')){

}

genedb <- list()
genedb[['all']] <- gene_xref
genedb[['alias2primary']] <- alias2primary
#genedb[['poorly_defined']] <- poorly_defined_genes
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


for(n in c('cancerdrugdb',
           'genedb',
           'hpa',
           'ligandreceptordb',
           'otdb',
           'pathwaydb',
           'pfamdomaindb',
           'projectscoredb',
           'projectsurvivaldb',
           'subcelldb',
           'tcgadb',
           'tftargetdb',
           'tissuecelldb')){



  }


rm(gene_xref)
rm(gencode)
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
rm(comppidb)
rm(cancer_drugs)
rm(gene_info)

rm(corum_complexes)
rm(alias2primary)
rm(genes_function_knowledge_level)
rm(ncbi_gene_summary)
rm(uniprot_gene_summary)
rm(tissue_cell_expr)
rm(opentarget_associations)
rm(ts_oncogene_annotations)
rm(go_annotations_qc)
rm(go_function_terms_prgene)
rm(go_process_terms_prgene)
rm(poorly_defined_genes)
rm(subcell_figure_legend)
rm(cell_key)
rm(subcelldb)
rm(cancer_hallmark_data)
rm(tf_target_interactions)
rm(ligandreceptordb)
rm(proteincomplexdb)

