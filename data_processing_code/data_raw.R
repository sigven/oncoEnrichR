library(TissueEnrich)
library(gganatogram)

source('data_processing_code/data_utility_functions.R')

msigdb_version <- 'v2022.1.Hs'
wikipathways_version <- "20230210"
netpath_version <- "2010"
opentargets_version <- "2022.11"
kegg_version <- "20230101"
gencode_version <- "42"
uniprot_release <- "2022_05"

## Which databases to update or retrieve from last updated state
update_omnipathdb <- F
update_hpa <- F
update_tcga <- F
update_cancer_hallmarks <- F
update_omnipath_regulatory <- F
update_omnipath_complexdb <- F
update_subcelldb <- F
update_ligand_receptor_db <- F

oe_version <- "1.4.0"

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
gOncoX <- list()
gOncoX[['basic']] <- geneOncoX::get_basic(
  cache_dir = data_raw_dir)
gOncoX[['alias']] <- geneOncoX::get_alias(
  cache_dir = data_raw_dir)

gOncoX[['basic']]$records <-
  gOncoX[['basic']]$records |>
  dplyr::mutate(entrezgene = as.integer(
    entrezgene
  ))


gOncoX[['gencode']] <- geneOncoX::get_gencode(
  cache_dir = data_raw_dir
)

for(b in c('grch37','grch38')){
  gOncoX[['gencode']]$records[[b]] <-
    gOncoX[['gencode']]$records[[b]] |>
    dplyr::mutate(entrezgene = as.integer(
      entrezgene
    ))
}

gene_info <- gOncoX[['basic']]$records |>
  dplyr::select(
    entrezgene,
    symbol,
    name,
    gene_biotype,
    hgnc_id
  )



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

####---BIOGRID----####
biogrid <- as.data.frame(
  readr::read_tsv(
    file.path(
      data_raw_dir,
      "biogrid",
      "BIOGRID-Physical.tsv.gz"),
    col_names =
    c("entrezgene_A","entrezgene_B",
      "method","pmid",
      "throughput"), show_col_types = F)) |>
  dplyr::arrange(entrezgene_A, entrezgene_B) |>
  dplyr::mutate(pmid =
    stringr::str_replace(
      stringr::str_trim(pmid),"PUBMED:","")) |>
  dplyr::mutate(pmid = dplyr::if_else(
    !stringr::str_detect(pmid,"^([0-9]{5,})$"),
    as.numeric(NA),
    as.numeric(pmid)
  )) |>
  dplyr::mutate(tmp_A = dplyr::if_else(
    .data$entrezgene_B < .data$entrezgene_A,
    .data$entrezgene_B,
    .data$entrezgene_A
  )) |>
  dplyr::mutate(tmp_B = dplyr::if_else(
    .data$entrezgene_B < .data$entrezgene_A,
    .data$entrezgene_A,
    .data$entrezgene_B
  )) |>
  dplyr::select(
    -c("entrezgene_A","entrezgene_B")
  ) |>
  dplyr::rename(
    entrezgene_A = "tmp_A",
    entrezgene_B = "tmp_B"
  ) |>
  dplyr::filter(
    entrezgene_A != entrezgene_B
  ) |>
  dplyr::distinct()


####---Cancer hallmark annotations----####
genedb[['cancer_hallmark']] <-
  get_cancer_hallmarks(
    opentargets_version = opentargets_version,
    raw_db_dir = data_raw_dir,
    gene_info = gene_info,
    update = update_cancer_hallmarks)

####----TSG/Oncogene/driver annotations---####
ts_oncogene_annotations <-
  geneOncoX:::assign_cancer_gene_roles(
    gox_basic = gOncoX[['basic']],
    min_sources_driver = 2) |>
  dplyr::select(
    entrezgene,
    tsg,
    tsg_confidence_level,
    tsg_support,
    oncogene,
    oncogene_confidence_level,
    oncogene_support,
    driver,
    driver_support,
    cancergene_evidence
  ) |>
  dplyr::rename(tumor_suppressor = tsg,
                cancer_driver = driver)

####--- OTP: Target-disease associations ----####
opentarget_associations <-
  get_opentarget_associations(
    raw_db_dir = data_raw_dir,
    min_num_sources = 2,
    min_overall_score = 0.05,
    direct_associations_only = T)

otdb <- quantify_gene_cancer_relevance(
  ot_associations = opentarget_associations,
  cache_dir = data_raw_dir)

genedb[['transcript_xref']] <- get_unique_transcript_xrefs(
  raw_db_dir = data_raw_dir,
  gene_oncox = gOncoX,
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
  gene_oncox = gOncoX,
  go_terms_pr_gene = go_terms_pr_gene,
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

####----COMPARTMENTS - subcellular compartments---####
subcelldb <- get_subcellular_annotations(
  raw_db_dir = data_raw_dir,
  transcript_xref_db = genedb[['transcript_xref']],
  update = update_subcelldb
)

####---- Project Score/DepMap ----####
depmapdb <- get_fitness_data_crispr(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info
)

####--- Cancer-KM-Survival (CSHL) ---####
survivaldb <- get_survival_associations(
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
oedb[['survivaldb']] <- survivaldb
oedb[['depmapdb']] <- depmapdb
oedb[['tcgadb']] <- tcgadb
oedb[['pathwaydb']] <- pathwaydb
oedb[['slparalogdb']] <- slparalogdb
oedb[['biogrid']] <- biogrid
#oedb[['synlethdb']] <- synlethdb

save(oedb, file="inst/internal_db/oedb.rda")

#googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))
gd_records <- list()
db_id_ref <- data.frame()


for(elem in c('cancerdrugdb',
           'genedb',
           'hpa',
           'ligandreceptordb',
           'otdb',
           'pfamdb',
           'pathwaydb',
           'depmapdb',
           'survivaldb',
           'subcelldb',
           'slparalogdb',
           'tcgadb',
           'tftargetdb',
           'tissuecelldb',
           'biogrid',
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
    size <- utils::object.size(oedb[[elem]][['compartments']])
    hsize <- R.utils::hsize.object_size(size)
  }

  google_rec_df$size <- as.character(size)
  google_rec_df$hsize <- hsize


  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)

}

usethis::use_data(db_id_ref, internal = T, overwrite = T)


####--- Clean-up ----####
rm(pfamdb)
rm(tcgadb)
rm(hpa)
rm(omnipathdb)
rm(otdb)
rm(pathwaydb)
rm(release_notes)
rm(cancerdrugdb)
rm(genedb)
rm(gene_info)
rm(depmapdb)
rm(survivaldb)
rm(tcga_maf_datasets)
rm(opentarget_associations)
rm(ts_oncogene_annotations)
rm(subcelldb)
rm(tf_target_interactions)
rm(ligandreceptordb)
rm(tissuecelldb)
rm(go_terms_pr_gene)
rm(slparalogdb)
rm(gOncoX)
rm(biogrid)

rm(oedb)

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
