
source('../oncoEnrichR/data_processing_code/data_utility_functions.R')

####--- Database versions and update flags ----####
msigdb_version <- 'v2025.1.Hs'
wikipathways_version <- "20250910"
netpath_version <- "2010"
opentargets_version <- "2025.09"
kegg_version <- "20250603"
gencode_version <- "48"
uniprot_release <- "2025_03"
biogrid_release <- "4.4.249"

## Which databases to update or retrieve from last updated state
db_updates <- list()
db_updates[['omnipathdb']] <- F
db_updates[['hpa']] <- F
db_updates[['ot']] <- F
db_updates[['tcga']] <- F
db_updates[['cancer_hallmarks']] <- F
db_updates[['omnipath_complexdb']] <- F
db_updates[['tftargetdb']] <- F
db_updates[['subcelldb']] <- F
db_updates[['ligand_receptor_db']] <- F
db_updates[['cellmodeldb']] <- F
db_updates[['biogrid']] <- F
db_updates[['pfamdb']] <- F

oe_version <- "1.6.0"

data_raw_dir <-
  "/Users/sigven/project_data/packages/package__oncoEnrichR/db/raw"
data_output_dir <-
  "/Users/sigven/project_data/packages/package__oncoEnrichR/db/output"

## oncoenrichr db object
oedb <- list()

## Initialize sub-lists for gene elements
oedb[['genedb']] <- list()
oedb[['genedb']][['cancer_hallmark']] <- list()
oedb[['genedb']][['proteincomplexdb']] <- list()
oedb[['genedb']][['transcript_xref']] <- list()
oedb[['genedb']][['all']] <- list()


## Release notes
oedb[['release_notes']] <- list()
software_db_version <-
  readr::read_tsv(
    file = "../oncoEnrichR/data_processing_code/RELEASE_NOTES.txt",
    skip = 1, col_names = T, show_col_types = F,
    comment = "#") |>
  dplyr::mutate(license_url = dplyr::case_when(
    license == "CC-BY 4.0" ~ "https://creativecommons.org/licenses/by/4.0/",
    license == "CC0 1.0" ~ "https://creativecommons.org/publicdomain/zero/1.0/",
    license == "MIT" ~ "https://en.wikipedia.org/wiki/MIT_License",
    license == "Apache 2.0" ~ "https://www.apache.org/licenses/LICENSE-2.0",
    license == "GPL v3.0" ~ "https://www.gnu.org/licenses/gpl-3.0.en.html",
    license == "CC-BY-SA 3.0" ~ "https://creativecommons.org/licenses/by-sa/3.0/",
    license == "Artistic-2.0" ~ "https://opensource.org/license/artistic-2-0/",
    license == "Open Access" ~ "Open Access",
    TRUE ~ as.character(".")
  )) |>
  dplyr::filter(
    name != "EFO" &
      name != "DiseaseOntology" &
      name != "CellTalkDB" &
      name != "GENCODE"
  )


i <- 1
while(i <= nrow(software_db_version)){
  oedb[['release_notes']][[software_db_version[i,]$key]] <-
    list('url' = software_db_version[i,]$url,
         'description' =
           software_db_version[i,]$description,
         'version' =
           software_db_version[i,]$version,
         'name' =
           software_db_version[i,]$name,
         'resource_type' =
           software_db_version[i,]$resource_type,
         'license' =
           software_db_version[i,]$license,
         'license_url' =
           software_db_version[i,]$license_url)
  i <- i + 1
}
rm(software_db_version)




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

####---GENCODE transcripts-----####
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
oedb[['pfamdb']] <- get_pfam_annotations(
  raw_db_dir = data_raw_dir,
  update = db_updates[['pfamdb']]
)

####---BIOGRID----####
oedb[['biogrid']] <- get_biogrid_physical_interactions(
  raw_db_dir = data_raw_dir,
  biogrid_release = biogrid_release,
  update = db_updates[['biogrid']]
)

####---Cancer hallmark annotations----####
oedb[['genedb']][['cancer_hallmark']] <-
  get_cancer_hallmarks(
    opentargets_version = opentargets_version,
    raw_db_dir = data_raw_dir,
    gene_info = gene_info,
    update = db_updates[['cancer_hallmarks']])

####----TSG/Oncogene/driver annotations---####
ts_oncogene_annotations <-
  geneOncoX:::assign_cancer_gene_roles(
    gox_basic = gOncoX[['basic']],
    min_sources_driver = 1) |>
  dplyr::filter(
    !(driver == T &
        driver_support == "CancerMine")
  ) |>
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
    min_overall_score = 0.02,
    release = opentargets_version,
    direct_associations_only = T,
    update = db_updates[['ot']])

oedb[['otdb']] <- quantify_gene_cancer_relevance(
  ot_associations = opentarget_associations,
  cache_dir = data_raw_dir)

####---Transcript xrefs ----####
oedb[['genedb']][['transcript_xref']] <-
  get_unique_transcript_xrefs(
    raw_db_dir = data_raw_dir,
    gene_oncox = gOncoX,
    update = T
  )


####---OmniPathDB - gene annotations ----####
omnipathdb <- get_omnipath_gene_annotations(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info,
  update = db_updates[['omnipathdb']]
)

####---OmniPathDB/Collectri - TF-target interactions ----####
oedb[['tftargetdb']] <- get_regulatory_collectri(
  raw_db_dir = data_raw_dir,
  update = db_updates[['tftargetdb']]
)

####---Pathway annotations ----####
oedb[['pathwaydb']] <- get_pathway_annotations(
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
oedb[['genedb']][['proteincomplexdb']] <- get_protein_complexes(
  raw_db_dir = data_raw_dir,
  update = db_updates[['omnipath_complexdb']]
)

####---OTP - Targeted cancer drugs ---####
oedb[['cancerdrugdb']] <- get_cancer_drugs(
  raw_db_dir = data_raw_dir
)

## Append all gene annotations to a single dataframe
## 1) Remove transcripts
## 2) Add gene name link
## 3) Add targeted drugs (early/late phase) and tractability category
## 4) Add tumor suppressor/oncogene annotations from CancerMine
## 5) Add gene function summary descriptions from NCBI/UniProt
## 6) Assign unknown function rank

oedb[['genedb']][['all']] <- generate_gene_xref_df(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info,
  transcript_xref_db = oedb[['genedb']][['transcript_xref']],
  ts_oncogene_annotations = ts_oncogene_annotations,
  opentarget_associations = opentarget_associations,
  gene_oncox = gOncoX,
  otdb = oedb[['otdb']],
  cancerdrugdb = oedb[['cancerdrugdb']],
  update = T
)

#oedb[['genedb']] <- genedb

####--- Ligand-Receptor interactions ----####
oedb[['ligandreceptordb']] <- get_ligand_receptors(
   raw_db_dir = data_raw_dir,
   keggdb = oedb[['pathwaydb']][['kegg']],
   update = db_updates[['ligand_receptor_db']]
)

####----COMPARTMENTS - subcellular compartments---####
oedb[['subcelldb']] <- get_subcellular_annotations(
  raw_db_dir = data_raw_dir,
  transcript_xref_db = oedb[['genedb']][['transcript_xref']],
  update = db_updates[['subcelldb']]
)

####---- Cell Model Passports/CRISPR Fitness ----####
oedb[['cellmodeldb']] <- get_fitness_data_CMP(
  raw_db_dir = data_raw_dir,
  update = db_updates[['cellmodeldb']],
  gene_info = gene_info
)

####--- Cancer-KM-Survival (CSHL) ---####
oedb[['survivaldb']] <- get_survival_associations(
  gene_info = gene_info,
  raw_db_dir = data_raw_dir
)

###--- Predicted SL paralogs ---####
oedb[['slparalogdb']] <- get_paralog_SL_predictions(
  raw_db_dir = data_raw_dir,
  gene_info = gene_info
)

####--- Human Protein Atlas ---####
oedb[['hpa']] <- get_hpa_associations(
  raw_db_dir = data_raw_dir,
  gene_xref = oedb[['genedb']][['all']],
  update = db_updates[['hpa']]
)

####----TCGA aberration data ----####
oedb[['tcgadb']] <- get_tcga_db(
  raw_db_dir = data_raw_dir,
  update = db_updates[['tcga']],
  gene_xref = oedb[['genedb']][['all']]
)

save(oedb, file = file.path(
  here::here(), "inst", "internal_db", "oedb.rda"))

####--- Clean-up ----####
rm(omnipathdb)
rm(gene_info)
rm(opentarget_associations)
rm(ts_oncogene_annotations)
rm(gOncoX)


####--- Color palettes ----####
color_palette <- list()

color_palette[['tissue']] <-
  colorspace::qualitative_hcl(28, palette = "Dark 2")
  #pals::stepped(30)

subcell_compartment_colors <-
  c("#F7FCB9",
    "#FFEDA0",
    "#FEC966",
    "#FBAD4C",
    "#F98F3B",
    "#F2672E",
    "#D8432B",
    "#B92C2E",
    "#981B2E",
    "#800026")

color_palette[['subcell_compartments']] <-
  data.frame(
    'bin' = seq(1:10),
    'fill' = subcell_compartment_colors)
usethis::use_data(color_palette, overwrite = T)

####--- Animal cell compartments ----####
# This is a TSV file with subcellular compartment names, IDs, and
# that are visible in the bscui map figure (https://www.swissbiopics.org/name/Animal_cell)
visible_compartments_swissbio <-
  readr::read_tsv(
    "data_processing_code/animal_cell_compartments.tsv",
    show_col_types = F)

####--- Subcellular compartment map ----####
bscui_map <- list()
bscui_map[['figure']] <-
  bscui::bscui(
    xml2::read_xml(system.file(
      "examples",
      "Animal_cells.svg.gz",
      package = "bscui"
    )))

bscui_map[['map']] <- readr::read_tsv(system.file(
  "examples",
  "uniprot_cellular_locations.txt.gz",
  package = "bscui"),
  col_types=strrep("c", 6)) |>
  dplyr::mutate(
    id = stringr::str_remove(
      `Subcellular location ID`, "-")) |>
  janitor::clean_names() |>
  dplyr::mutate(
    gene_ontologies = stringr::str_replace(
      gene_ontologies,"GO:","GO_")) |>
  tidyr::separate(
    gene_ontologies, c("go_id","go_term"),
    sep=":", remove = T) |>
  dplyr::mutate(go_id = stringr::str_replace(
    go_id,"_",":")) |>
  dplyr::filter(!is.na(go_id)) |>

  dplyr::semi_join(
    visible_compartments_swissbio,
    by = "name"
  ) |>

  ## ignore extracellular space and secreted, as
  ## well as very generic compartmental concepts
  dplyr::filter(
    id != "SL0112" & # Extracellular space
      id != "SL0243" & # Secreted
      id != "SL0086" & # Cytoplasm
      id != "SL0162" & # Membrane
      id != "SL0191" & # Nucleus
      id != "SL0209" & # Plastid
      id != "SL0038" & # Cell junction
      id != "SL0280" & # Cell projection
      #id != "" & # Cytoplasmic vesicle
      id != "SL0244" # Secretory vesicle
      #id != "" # Host vacuole

  ) |>
  dplyr::select(-c("go_term","category"))

bscui_map[['ui_elements']] <-
  bscui_map[['map']] |>
  dplyr::mutate(
    ui_type = "selectable",
    title = glue::glue(
      '<div style="width:350px; height:200px; ',
      'overflow:auto; padding:5px;',
      'font-size:75%;',
      'border:black 1px solid; background:#FFFFF0AA;">',
      "<strong>{name}</strong>: {description}",
      "</div>",
      .sep=" "
    )
  ) |>
  dplyr::select(id, ui_type, title)

subcell_map <- bscui_map
usethis::use_data(subcell_map, overwrite = T)

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
           'cellmodeldb',
           'survivaldb',
           'subcelldb',
           'slparalogdb',
           'tcgadb',
           'tftargetdb',
           #'tissuecelldb',
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

  local_rds_fpath <- file.path(
    data_output_dir, paste0("v",oe_version),
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

  size <- utils::object.size(oedb[[elem]])
  hsize <- R.utils::hsize.object_size(size)

  google_rec_df$size <- as.character(size)
  google_rec_df$hsize <- hsize

  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)

}

usethis::use_data(db_id_ref, internal = T, overwrite = T)

cp_output_cols <-
  c('db',
    'description',
    'standard_name',
    'gene_ratio',
    'background_ratio',
    'enrichment_factor',
    'rich_factor',
    'z_score',
    'pvalue',
    'p.adjust',
    'qvalue',
    'gene_id',
    'gene_symbol',
    'gene_symbol_link',
    'setting_p_value_cutoff',
    'setting_q_value_cutoff',
    'setting_p_value_adj_method',
    'setting_min_geneset_size',
    'setting_max_geneset_size',
    'description_link',
    'exact_source',
    'external_url',
    'url')

usethis::use_data(cp_output_cols, overwrite = T)

base_urls <- list()
base_urls[['otp']] <- "https://platform.opentargets.org"
base_urls[['ncbi']] <- "https://www.ncbi.nlm.nih.gov"
base_urls[['wpway']] <- "https://www.wikipathways.org"
base_urls[['kegg']] <- "https://www.genome.jp"
base_urls[['cmp']] <- "https://cellmodelpassports.sanger.ac.uk"
base_urls[['ensembl']] <- "https://www.ensembl.org/Homo_sapiens"
base_urls[['string']] <- "https://string-db.org"
base_urls[['hpa']] <- "https://www.proteinatlas.org"
base_urls[['pfam']] <- "https://www.ebi.ac.uk/interpro"

usethis::use_data(base_urls, overwrite = T)

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
#

renv_packages <- jsonlite::fromJSON("renv.lock")
#for(c in db_packages){
#  renv_packages$Packages[[c]] <- NULL
#}
docker_renv <- jsonlite::toJSON(
  renv_packages, flatten = T, auto_unbox = T) |>
  jsonlite::prettify()
sink(file = "docker/renv.lock")
cat(docker_renv)
sink()

devtools::build(path = "docker", vignettes = F, manual = F)
