library(TissueEnrich)
library(gganatogram)

source('data_processing_code/data_utility_functions.R')

msigdb_version <- '2022.1'
wikipathways_version <- "20220810"
netpath_version <- "2010"
opentargets_version <- "2022.06"
kegg_version <- "20220809"
gencode_version <- "41"
uniprot_release <- "2022_03"

## Which databases to update or retrieve from last updated state
update_omnipathdb <- F
update_hpa <- F
update_ncbi_gene_summary <- T
update_project_score <- F
update_project_survival <- F
update_tcga <- F
update_cancer_hallmarks <- F
update_omnipath_regulatory <- F
update_omnipath_complexdb <- F
update_gencode <- F
update_ligand_receptor_db <- T

oe_version <- "1.3.0"

data_raw_dir <- "/Users/sigven/project_data/package__oncoEnrichR/db/raw"
data_output_dir <- "/Users/sigven/project_data/package__oncoEnrichR/db/output"

software_db_version <-
  read.table(file="data_processing_code/RELEASE_NOTES.txt",
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
genedb <- list()


####---Gene info----####
gene_info <-
  get_gene_info_ncbi(
    raw_db_dir = data_raw_dir)

####---Protein domains----####
pfamdb <- as.data.frame(
  readr::read_tsv(
    file.path(data_raw_dir,
              "pfam",
              "pfam.uniprot.tsv.gz"),
    show_col_types = F) |>
    dplyr::select(-uniprot_acc) |>
    dplyr::rename(uniprot_acc = uniprot_acc_noversion) |>
    ## reviewed accessions only
    dplyr::filter(nchar(uniprot_acc) == 6) |>
    dplyr::group_by(uniprot_acc,
                    pfam_id,
                    pfam_short_name,
                    pfam_long_name) |>
    dplyr::summarise(
      domain_freq = dplyr::n(),
      .groups = "drop")
)

####---Cancer hallmark annotations----####
genedb[['cancer_hallmark']] <-
  get_cancer_hallmarks(
    opentargets_version = opentargets_version,
    raw_db_dir = data_raw_dir,
    gene_info = gene_info,
    update = update_cancer_hallmarks)

####---dbNSFP: gene summary descriptions---####
uniprot_gene_summary <-
  get_dbnsfp_gene_annotations(
    raw_db_dir = data_raw_dir)

####----CancerMine/NCG---####
ts_oncogene_annotations <-
  get_ts_oncogene_annotations(
    raw_db_dir = data_raw_dir,
    gene_info = gene_info,
    version = "47") |>
  dplyr::select(
    entrezgene, tumor_suppressor,
    oncogene, citation_links_oncogene,
    citation_links_tsgene, cancergene_support,
    citation_links_cdriver, cancer_driver)

####--- OTP: Target-disease associations ----####
opentarget_associations <-
  get_opentarget_associations(
    raw_db_dir = data_raw_dir,
    min_num_sources = 2,
    min_overall_score = 0.05,
    direct_associations_only = T)

otdb <- quantify_gene_cancer_relevance(
  ot_associations = opentarget_associations)

####---GENCODE----####
gencode <- list()
for(build in c('grch37','grch38')){

  gencode_v <- gencode_version
  if(build == 'grch37'){
    gencode_v <- "19"
  }

  gencode[[build]] <- get_gencode_data(
    raw_db_dir = data_raw_dir,
    gene_info = gene_info,
    build = build,
    gencode_version = gencode_v,
    update = update_gencode
  )

}

genedb[['transcript_xref']] <- get_unique_transcript_xrefs(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info,
  gencode = gencode,
  update = T
)


####---OmniPathDB - gene annotations ----####
omnipathdb <- get_omnipath_gene_annotations(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info,
  update = update_omnipathdb
)

####---OmniPathDB - TF interactions ----####
tf_target_interactions <- get_tf_target_interactions(
  raw_db_dir = data_raw_dir,
  update = update_omnipath_regulatory
)

####---Pathway annotations ----####
pathwaydb <- get_pathway_annotations(
  raw_db_dir = data_raw_dir,
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
  raw_db_dir = data_raw_dir,
  update = update_omnipath_complexdb
)

## Gene summaries from NCBI (provided by RefSeq/OMIM):
ncbi_gene_summary <- get_function_summary_ncbi(
  raw_db_dir = data_raw_dir,
  gene_df = gene_info,
  update = update_ncbi_gene_summary)

####---OTP - Targeted cancer drugs ---####
cancerdrugdb <- get_cancer_drugs(
  raw_db_dir = data_raw_dir
)

####--- Gene Ontology ---####
go_terms_pr_gene <- get_gene_go_terms(
  raw_db_dir = data_raw_dir
)

## Append all gene annotations to a single dataframe
## 1) Remove transcripts
## 2) Add gene name link
## 3) Add targeted drugs (early/late phase) and tractability category
## 4) Add tumor suppressor/oncogene annotations from CancerMine
## 5) Add gene function summary descriptions from NCBI/UniProt
## 6) Assign unknown function rank

genedb[['all']] <- generate_gene_xref_df(
  raw_db_dir = data_raw_dir,
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
  raw_db_dir = data_raw_dir,
  keggdb = pathwaydb$kegg,
  update = update_ligand_receptor_db
)

####----ComPPI - subcellular compartments---####
subcelldb <- get_subcellular_annotations(
  raw_db_dir = data_raw_dir,
  transcript_xref_db = genedb[['transcript_xref']]
)

####---- Project Score/CRISPR ----####
projectscoredb <- get_fitness_data_crispr(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info
)

####--- Cancer-KM-Survival (CSHL) ---####
projectsurvivaldb <- get_survival_associations(
  gene_info = gene_info,
  raw_db_dir = data_raw_dir
)

###--- Predicted SL paralogs ---####
slparalogdb <- get_paralog_SL_predictions(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info
)

# synlethdb <- get_synthetic_lethality_pairs(
#   raw_db_dir = data_raw_dir
# )


####--- Human Protein Atlas ---####
hpa <- get_hpa_associations(
  raw_db_dir = data_raw_dir,
  gene_xref = genedb[['all']],
  update = update_hpa
)

####----TCGA aberration data ----####
tcgadb <- get_tcga_db(
  raw_db_dir = data_raw_dir,
  update = update_tcga,
  gene_xref = genedb[['all']]
)

tissuecelldb <- get_tissue_celltype_specificity(
  raw_db_dir = data_raw_dir
)

#db_props <- data.frame()

oedb <- list()
oedb[['cancerdrugdb']] <- cancerdrugdb
oedb[['release_notes']] <- release_notes
oedb[['subcelldb']] <- subcelldb
oedb[['ligandreceptordb']] <- ligandreceptordb
oedb[['genedb']] <- genedb
oedb[['otdb']] <- otdb
oedb[['pfamdb']] <- pfamdb
oedb[['tftargetdb']] <- tf_target_interactions
oedb[['tissuecelldb']] <- tissuecelldb
oedb[['hpa']] <- hpa
oedb[['projectsurvivaldb']] <- projectsurvivaldb
oedb[['projectscoredb']] <- projectscoredb
oedb[['tcgadb']] <- tcgadb
oedb[['pathwaydb']] <- pathwaydb
oedb[['slparalogdb']] <- slparalogdb
#oedb[['synlethdb']] <- synlethdb


googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))
gd_records <- list()
db_id_ref <- data.frame()


for(elem in c('cancerdrugdb',
           'genedb',
           'hpa',
           'ligandreceptordb',
           'otdb',
           'pfamdb',
           'pathwaydb',
           'projectscoredb',
           'projectsurvivaldb',
           'subcelldb',
           'slparalogdb',
           'tcgadb',
           'tftargetdb',
           'tissuecelldb',
           'release_notes')){

  if(!dir.exists(
    file.path(data_output_dir, paste0("v",oe_version)))){
    dir.create(
      file.path(
        data_output_dir,
        paste0("v",oe_version))
      )
  }

  local_rds_fpath <- file.path(data_output_dir, paste0("v",oe_version),
                           paste0(elem,"_v", oe_version, ".rds"))
  saveRDS(oedb[[elem]],
          file = local_rds_fpath)

  (gd_records[[elem]] <- googledrive::drive_upload(
    local_rds_fpath,
    paste0("oncoEnrichR_DB/", elem, "_v", oe_version,".rds")
  ))

  google_rec_df <-
    dplyr::select(
      as.data.frame(gd_records[[elem]]), name, id) |>
    dplyr::rename(
      gid = id,
      filename = name) |>
    dplyr::mutate(
      name =
        stringr::str_replace(filename,"_v\\S+$",""),
      date = as.character(Sys.Date()),
      pVersion = oe_version) |>
    dplyr::mutate(
      md5Checksum =
        gd_records[[elem]]$drive_resource[[1]]$md5Checksum)

  size <- NA
  hsize <- NA

  if(elem != "subcelldb"){
    size <- utils::object.size(oedb[[elem]])
    hsize <- R.utils::hsize.object_size(size)
  }else{
    size <- utils::object.size(oedb[[elem]][['comppidb']])
    hsize <- R.utils::hsize.object_size(size)
  }

  google_rec_df$size <- as.character(size)
  google_rec_df$hsize <- hsize


  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)

}

usethis::use_data(db_id_ref, internal = T, overwrite = T)



# for(n in c('cancerdrugdb',
#            'genedb',
#            'hpa',
#            'ligandreceptordb',
#            'otdb',
#            'pfamdb',
#            'pathwaydb',
#            'projectscoredb',
#            'projectsurvivaldb',
#            'subcelldb',
#            'slparalogdb',
#            'tcgadb',
#            'tftargetdb',
#            'tissuecelldb',
#            'release_notes')){
#
#
#   checksum_db <- NA
#   size <- NA
#   hsize <- NA
#   if(n != "subcelldb"){
#     checksum_db <- R.cache::getChecksum(oedb[[n]])
#     size <- utils::object.size(oedb[[n]])
#     hsize <- R.utils::hsize.object_size(size)
#   }else{
#     checksum_db <- R.cache::getChecksum(oedb[[n]][['comppidb']])
#     size <- utils::object.size(oedb[[n]][['comppidb']])
#     hsize <- R.utils::hsize.object_size(size)
#   }
#
#   db_entry <- data.frame('name' = n,
#                           'size' = as.character(size),
#                           'hsize' = as.character(hsize),
#                           'checksum' = checksum_db,
#                           'version' = oe_version,
#                           'date' = Sys.Date(),
#                          stringsAsFactors = F
#                         )
#   db_props <- db_props |>
#     dplyr::bind_rows(db_entry)
#
#
#   if(!dir.exists(
#     file.path(data_output_dir, paste0("v",oe_version)))){
#     system(paste0('mkdir ', file.path(
#       data_output_dir,
#       paste0("v",oe_version)))
#     )
#   }
#
#   saveRDS(
#     oedb[[n]],
#     file = file.path(
#       data_output_dir,
#       paste0("v",oe_version),
#       paste0(n,'.rds'))
#     )
#
# }

save(oedb, file="inst/internal_db/oedb.rda")

####--- Clean-up ----####
rm(pfamdb)
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

####---Zenodo upload ----####

# zenodo_files_for_upload <-
#   list.files(file.path(data_output_dir, paste0("v",oe_version)),
#              full.names = T)
# zenodo <- zen4R::ZenodoManager$new(
#   token = Sys.getenv("ZENODO_TOKEN"),
#   logger = "INFO"
# )
#
# oedb_rec <- zenodo$getDepositionByDOI("10.5281/zenodo.6700715")
# oedb_rec$setPublicationDate(publicationDate = Sys.Date())
# oedb_rec$setVersion(version = paste0("v",oe_version))
# oedb_rec <- zenodo$depositRecordVersion(
#   oedb_rec,
#   delete_latest_files = TRUE,
#   files = zenodo_files_for_upload,
#   publish = FALSE)
#
# #zenodo$discardChanges(oedb_rec$id)
# db_props$zenodo_doi <- oedb_rec$metadata$doi

#usethis::use_data(db_props, overwrite = T)

#oedb_rec <- zenodo$publishRecord(oedb_rec$id)

db_packages <-
  c('CellChat',
    'rWikiPathways',
    'ComplexHeatmap',
    'Rtsne',
    'expm',
    'irlba',
    'pbapply',
    'cowplot',
    'reticulate',
    'RSpectra',
    'sna',
    'FNN',
    'tweenr',
    'reticulate',
    'TCGAbiolinksGUI.data',
    'TCGAbiolinks',
    'rlogging',
    'zen4R',
    'biomaRt',
    'redland',
    'rdflib',
    'magrittr',
    'log4r',
    'keyring',
    'pharmaOncoX',
    'oncoPhenoMap')

renv_packages <- jsonlite::fromJSON("renv.lock")
for(c in db_packages){
  renv_packages$Packages[[c]] <- NULL
}
docker_renv <- jsonlite::toJSON(
  renv_packages, flatten = T, auto_unbox = T) |>
  jsonlite::prettify()
sink(file = "docker/renv.lock")
cat(docker_renv)
sink()

devtools::build(path = "docker", vignettes = F, manual = F)
