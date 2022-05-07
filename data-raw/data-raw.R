library(magrittr)
library(TissueEnrich)
library(gganatogram)

msigdb_version <- 'v7.5.1 (Jan 2022)'
wikipathways_version <- "20220410"
netpath_version <- "2010"
opentargets_version <- "2022.04"
kegg_version <- "20220324"
gencode_version <- "40"
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
update_ligand_receptor_db <- T


uniprot_release <- "2022_01"

software_db_version <-
  read.table(file="data-raw/RELEASE_NOTES.txt",
             skip = 1, sep = "\t", stringsAsFactors = F,
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
#usethis::use_data(release_notes, overwrite = T)
genedb <- list()

source('data-raw/data_utility_functions.R')

####---Gene info----####
gene_info <-
  get_gene_info_ncbi(
    basedir = here::here())

####---Cancer hallmark annotations----####
genedb[['cancer_hallmark']] <-
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
    version = "45") %>%
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

genedb[['transcript_xref']] <- get_unique_transcript_xrefs(
  basedir = here::here(),
  gene_info = gene_info,
  gencode = gencode,
  update = T
)


####---OmniPathDB - gene annotations ----####
omnipathdb <- get_omnipath_gene_annotations(
  basedir = here::here(),
  gene_info = gene_info,
  update = update_omnipathdb
)

####---OmniPathDB - TF interactions ----####
tf_target_interactions <- get_tf_target_interactions(
  basedir = here::here(),
  update = update_omnipath_regulatory
)

####---Pathway annotations ----####
pathwaydb <- get_pathway_annotations(
  basedir = here::here(),
  gene_info = gene_info,
  wikipathways_version = wikipathways_version,
  kegg_version = kegg_version,
  netpath_version = netpath_version,
  msigdb_version = msigdb_version,
  omnipathdb = omnipathdb,
  msigdb = T,
  netpath = T,
  wikipathways = T,
  kegg = T
)

####---OmniPathDB - protein complexes ----####
genedb[['proteincomplexdb']] <- get_protein_complexes(
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

genedb[['all']] <- generate_gene_xref_df(
  basedir = here::here(),
  gene_info = gene_info,
  transcript_xref_db = genedb[['transcript_xref']],
  ts_oncogene_annotations = ts_oncogene_annotations,
  opentarget_associations = opentarget_associations,
  go_terms_pr_gene = go_terms_pr_gene,
  ncbi_gene_summary = ncbi_gene_summary,
  uniprot_gene_summary = uniprot_gene_summary,
  otdb = otdb,
  cancerdrugdb = cancerdrugdb,
  update = T
)

####--- Ligand-Receptor interactions ----####
ligandreceptordb <- get_ligand_receptors(
  basedir = here::here(),
  keggdb = pathwaydb$kegg,
  update = update_ligand_receptor_db
)

####----ComPPI - subcellular compartments---####
subcelldb <- get_subcellular_annotations(
  basedir = here::here(),
  transcript_xref_db = genedb[['transcript_xref']]
)

####---- Project Score/CRISPR ----####
projectscoredb <- get_fitness_data_crispr(
  basedir = here::here(),
  gene_info = gene_info
)

####--- Cancer-KM-Survival (CSHL) ---####
projectsurvivaldb <- get_survival_associations(
  gene_info = gene_info,
  basedir = here::here()
)

###--- Predicted SL paralogs ---####
slparalogdb <- get_paralog_SL_predictions(
  basedir = here::here(),
  gene_info = gene_info
)

# synlethdb <- get_synthetic_lethality_pairs(
#   basedir = here::here()
# )


####--- Human Protein Atlas ---####
hpa <- get_hpa_associations(
  basedir = here::here(),
  gene_xref = genedb[['all']],
  update = update_hpa
)

####----TCGA aberration data ----####
tcgadb <- get_tcga_db(
  basedir = here::here(),
  update = update_tcga,
  gene_xref = genedb[['all']]
)

tissuecelldb <- get_tissue_celltype_specificity(
  basedir = here::here()
)

db_props <- data.frame()

oedb <- list()
oedb[['cancerdrugdb']] <- cancerdrugdb
oedb[['release_notes']] <- release_notes
oedb[['subcelldb']] <- subcelldb
oedb[['ligandreceptordb']] <- ligandreceptordb
oedb[['genedb']] <- genedb
oedb[['otdb']] <- otdb
oedb[['tftargetdb']] <- tf_target_interactions
oedb[['tissuecelldb']] <- tissuecelldb
oedb[['hpa']] <- hpa
oedb[['projectsurvivaldb']] <- projectsurvivaldb
oedb[['projectscoredb']] <- projectscoredb
oedb[['tcgadb']] <- tcgadb
oedb[['pathwaydb']] <- pathwaydb
oedb[['slparalogdb']] <- slparalogdb
#oedb[['synlethdb']] <- synlethdb


for(n in c('cancerdrugdb',
           'genedb',
           'hpa',
           'ligandreceptordb',
           'otdb',
           'pathwaydb',
           'projectscoredb',
           'projectsurvivaldb',
           'subcelldb',
           'slparalogdb',
           #'synlethdb',
           'tcgadb',
           'tftargetdb',
           'tissuecelldb',
           'release_notes')){


  checksum_db <- NA
  size <- NA
  hsize <- NA
  if(n != "subcelldb"){
    checksum_db <- R.cache::getChecksum(oedb[[n]])
    size <- utils::object.size(oedb[[n]])
    hsize <- R.utils::hsize.object_size(size)
  }else{
    checksum_db <- R.cache::getChecksum(oedb[[n]][['comppidb']])
    size <- utils::object.size(oedb[[n]][['comppidb']])
    hsize <- R.utils::hsize.object_size(size)
  }

  db_entry <- data.frame('name' = n,
                          'size' = as.character(size),
                          'hsize' = as.character(hsize),
                          'checksum' = checksum_db,
                          'version' = '1.0.9',
                          'date' = Sys.Date(),
                         stringsAsFactors = F
                        )
  db_props <- db_props %>%
    dplyr::bind_rows(db_entry)

  saveRDS(oedb[[n]], file = paste0('db/',n,'.rds'))

}

save(oedb, file="inst/internal_db/oedb.rda")

usethis::use_data(db_props, overwrite = T)

# tumor_types <- data.frame(
#   code = 1, name = 'Adrenal Gland', tcga_code = 'ACC', project_score_code = NA, hpa_code = "adrenal_gland", stringsAsFactors = F
# )
# tumor_types <- tumor_types %>%
#   dplyr::bind_rows(
#     data.frame(code = 2, name = 'Biliary Tract', tcga_code = "CHOL", project_score_code = NA, hpa_code = NA)) %>%
#   dplyr::bind_rows(
#     data.frame(code = 3, name = 'Bladder/Urinary Tract', tcga_code = 'BLCA', hpa_code = 'urinary_bladder')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 4, name = 'Bone', tcga_code = NA,  project_score_code = 'Bone')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 5, name = 'Breast', tcga_code = 'BRCA', project_score_code = 'Breast', hpa_code = 'breast')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 6, name = 'Cervix', tcga_code = 'CESC', project_score_code = NA, hpa_code = 'cervix')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 7, name = 'CNS/Brain', tcga_code = 'GBM|LGG', project_score_code = 'CNS/Brain',
#                hpa_code = 'cerebral_cortex|cerebellum|hypothalamus|substantia_nigra|hippocampus|amygdala|nucleus_accumbens|caudate|putamen|')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 8, name = 'Colon/Rectum', tcga_code = 'COAD|READ', project_score_code = 'Colon/Rectum',
#                hpa_code = 'colon')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 9, name = 'Esophagus/Stomach', tcga_code = 'ESCA|STAD', project_score_code = 'Esophagus|Stomach',
#                hpa_code = 'small_intestine|esophagus|stomach')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 10, name = 'Eye', tcga_code = 'UVM', hpa_code = 'retina')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 11, name = 'Head and Neck', tcga_code = 'HNSC', project_score_code = "Head and Neck")) %>%
#   dplyr::bind_rows(
#     data.frame(code = 12, name = 'Kidney', tcga_code = 'KICH|KIRC|KIRP', project_score_code = "Kidney", hpa_code = 'kidney')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 13, name = 'Liver', tcga_code = 'LIHC', hpa_code = 'liver')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 14, name = 'Lung', tcga_code = 'LUAD|LUSC', project_score_code = "Lung", hpa_code = 'lung')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 15, name = 'Lymphoid', tcga_code = 'DLBC', project_score_code = "Haematopoietic and Lymphoid")) %>%
#   dplyr::bind_rows(
#     data.frame(code = 16, name = 'Myeloid', tcga_code = 'LAML', project_score_code = "Haematopoietic and Lymphoid")) %>%
#   dplyr::bind_rows(
#     data.frame(code = 17, name = 'Ovary/Fallopian Tube', tcga_code = 'OV', project_score_code = "Ovary",
#                hpa_code = 'ovary|fallopian_tube')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 18, name = 'Pancreas', tcga_code = 'PRAD', project_score_code = 'Pancreas',
#                hpa_code = 'pancreas')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 19, name = 'Peripheral Nervous System', tcga_code = NA, project_score_code = 'Peripheral Nervous System')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 20, name = 'Pleura', tcga_code = 'MESO')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 21, name = 'Prostate', tcga_code = 'PRAD', hpa_code = 'prostate')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 22, name = 'Skin', tcga_code = 'SKCM', project_score_code = "Skin", hpa_code = 'skin')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 23, name = 'Soft Tissue', tcga_code = 'SARC')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 24, name = 'Testis', tcga_code = 'TGCT')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 25, name = 'Thymus', tcga_code = 'THYM')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 26, name = 'Thyroid', tcga_code = 'THCA', hpa_code = 'thyroid_gland')) %>%
#   dplyr::bind_rows(
#     data.frame(code = 27, name = 'Uterus', tcga_code = 'UCEC', project_score_code = "Endometrium", hpa_code = "endometrium")) %>%
#   dplyr::bind_rows(
#     data.frame(code = 28, name = 'Vulva/Vagina', hpa_code = 'vagina'))
#

rm(tcgadb)
rm(gencode)
rm(hpa)
rm(omnipathdb)
rm(otdb)
rm(pathwaydb)
rm(release_notes)
rm(cancerdrugdb)
rm(genedb)
rm(gene_info)
rm(ncbi_gene_summary)
rm(uniprot_gene_summary)
rm(projectscoredb)
rm(projectsurvivaldb)
rm(opentarget_associations)
rm(ts_oncogene_annotations)
rm(subcelldb)
rm(tf_target_interactions)
rm(ligandreceptordb)
rm(tissuecelldb)
rm(go_terms_pr_gene)
rm(oedb)
rm(slparalogdb)


