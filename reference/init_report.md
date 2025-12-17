# Function that initiates oncoEnrichR report structure

Function that initiates oncoEnrichR report structure

## Usage

``` r
init_report(
  oeDB,
  project_title = "_Project title_",
  project_owner = "_Project owner_",
  query_id_type = "symbol",
  ignore_id_err = TRUE,
  project_description = "_Project description_",
  ppi_string_min_score = 0.9,
  ppi_string_network_type = "functional",
  ppi_biogrid_min_evidence = 3,
  ppi_add_nodes = 30,
  ppi_show_isolated_nodes = FALSE,
  ppi_node_shadow = TRUE,
  ppi_show_drugs = T,
  bgset_description = "All annotated protein-coding genes",
  bgset_id_type = "symbol",
  enrichment_p_value_cutoff = 0.05,
  enrichment_p_value_adj = "BH",
  enrichment_q_value_cutoff = 0.2,
  enrichment_min_geneset_size = 10,
  enrichment_max_geneset_size = 500,
  enrichment_plot_num_terms = 20,
  enrichment_simplify_go = T,
  subcellcomp_min_confidence = 3,
  subcellcomp_min_channels = 1,
  regulatory_min_resources = 2,
  fitness_max_score = -2,
  show_ppi = T,
  show_disease = T,
  show_top_diseases_only = T,
  show_cancer_hallmarks = T,
  show_drug = T,
  show_enrichment = T,
  show_aberration = T,
  show_coexpression = T,
  show_subcell_comp = T,
  show_fitness = T,
  show_cell_tissue = F,
  show_ligand_receptor = T,
  show_regulatory = T,
  show_unknown_function = T,
  show_prognostic = T,
  show_synleth = T,
  show_complex = T,
  show_domain = T
)
```

## Arguments

- oeDB:

  oncoEnrichR annotation database object

- project_title:

  project title (title of report)

- project_owner:

  name of project owner \#param html_floating_toc logical - float the
  table of contents to the left of the main document content. The
  floating table of contents will always be visible even when the
  document is scrolled \#param html_report_theme Bootswatch theme for
  HTML report (any of
  "bootstrap","cerulean","cosmo","default","flatly","journal","lumen","paper","sandstone","simplex","spacelab","united","yeti")

- query_id_type:

  character indicating source of query (one of "uniprot_acc", "symbol",
  "entrezgene", or
  "ensembl_gene","ensembl_mrna","refseq_transcript_id","ensembl_protein","refseq_protein")

- ignore_id_err:

  logical indicating if analysis should continue when uknown query
  identifiers are encountered

- project_description:

  project background information

- ppi_string_min_score:

  minimum score (between 0 and 1) for confidence of retrieved
  protein-protein interactions (STRING)

- ppi_string_network_type:

  STRING network type ('physical' or 'functional')

- ppi_biogrid_min_evidence:

  minimum number of evidence items for protein-protein interactions
  shown (BIOGRID)

- ppi_add_nodes:

  number of nodes to add to target set when computing the
  protein-protein interaction network (STRING/BIOGRID)

- ppi_show_isolated_nodes:

  logical indicating if targets/nodes without any interactions should be
  displayed in the protein-protein interaction network

- ppi_node_shadow:

  show shadow for nodes in the displayed PPI network (STRING/BIOGRID)

- ppi_show_drugs:

  logical indicating if targeted drugs (\> phase 3) should be displayed
  in protein-protein interaction network

- bgset_description:

  character indicating type of background (e.g. "All lipid-binding
  proteins (n = 200)")

- bgset_id_type:

  character indicating source of query (one of "uniprot_acc", "symbol",
  "entrezgene", or
  "ensembl_gene","ensembl_mrna","refseq_transcript_id","ensembl_protein","refseq_protein")

- enrichment_p_value_cutoff:

  cutoff p-value for enrichment/over-representation analysis

- enrichment_p_value_adj:

  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
  "none"

- enrichment_q_value_cutoff:

  cutoff q-value for enrichment analysis

- enrichment_min_geneset_size:

  minimal size of geneset annotated by term for testing in
  enrichment/over-representation analysis

- enrichment_max_geneset_size:

  maximal size of geneset annotated by term for testing in
  enrichment/over-representation analysis

- enrichment_plot_num_terms:

  number of top enriched Gene Ontology terms (max) to show in enrichment
  barplot

- enrichment_simplify_go:

  remove highly similar GO terms in results from GO
  enrichment/over-representation analysis

- subcellcomp_min_confidence:

  minimun confidence level for subcellular compartment annotation in
  COMPARTMENTS (min = 3, max = 5)

- subcellcomp_min_channels:

  minimum number of channels that support a subcellular annotation in
  COMPARTMENTS

- regulatory_min_resources:

  minimum resources supporting regulatory interactions (TF-target)
  retrieved from CollecTRI (min = 1, max = 6)

- fitness_max_score:

  maximum loss-of-fitness score (scaled Bayes factor from BAGEL) for
  genes retrieved from Project Score

- show_ppi:

  logical indicating if report should contain protein-protein
  interaction data (STRING)

- show_disease:

  logical indicating if report should contain disease associations (Open
  Targets Platform, association_score \>= 0.05, support from at least
  two data types)

- show_top_diseases_only:

  logical indicating if report should contain top (n = 20) disease
  associations only pr. query gene (Open Targets Platform)

- show_cancer_hallmarks:

  logical indicating if report should contain annotations/evidence of
  cancer hallmarks per query gene (COSMIC/Open Targets Platform)

- show_drug:

  logical indicating if report should contain targeted cancer drug
  information

- show_enrichment:

  logical indicating if report should contain functional
  enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME,
  NetPath, WikiPathways)

- show_aberration:

  logical indicating if report should contain TCGA aberration plots
  (amplifications/deletions)

- show_coexpression:

  logical indicating if report should contain TCGA co-expression data
  (RNAseq) of query set with oncogenes/tumor suppressor genes

- show_subcell_comp:

  logical indicating if report should provide subcellular compartment
  annotations (COMPARTMENTS)

- show_fitness:

  logical indicating if report should provide fitness scores and target
  priority scores from CRISPR/Cas9 loss-of-fitness screens (Project
  Score)

- show_cell_tissue:

  logical indicating if report should contain tissue-specificity and
  single cell-type specificity assessments (Human Protein Atlas) of
  target genes, using data from the Human Protein Atlas

- show_ligand_receptor:

  logical indicating if report should contain ligand-receptor
  interactions (CellChatDB)

- show_regulatory:

  logical indicating if report should contain data on transcription
  factor (TF) - target interactions relevant for the query set
  (DoRothEA)

- show_unknown_function:

  logical indicating if report should highlight target genes with
  unknown or poorly defined functions (GO/Uniprot KB/NCBI)

- show_prognostic:

  logical indicating if mRNA-based (single-gene) prognostic associations
  to cancer types should be listed (Human Protein Atlas/TCGA)

- show_synleth:

  logical indicating if report should list overlap with predicted
  synthetic lethality interactions (gene paralogs only, De Kegel et al.,
  Cell Systems, 2021)

- show_complex:

  logical indicating if report should provide target memberships in
  known protein complexes (ComplexPortal/Compleat/PDB/CORUM)

- show_domain:

  logical indicating if report should provide target memberships in
  known protein domains (Pfam)
