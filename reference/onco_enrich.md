# Interrogate a gene list for cancer relevance

Function that interrogates a list of human gene identifiers for cancer
relevance. Multiple perspectives are offered, including tumor aberration
and co-expression patterns, druggability, protein-protein interactions,
gene fitness effects, regulatory interactions, subcellular compartment
enrichment, pathway enrichment, synthetic lethality interactions,
prognostic associations, and more.

## Usage

``` r
onco_enrich(
  query = NULL,
  oeDB = NULL,
  query_id_type = "symbol",
  ignore_id_err = TRUE,
  project_title = "_Project Title_",
  project_owner = "_Project owner_",
  project_description = "_Project description_",
  bgset = NULL,
  bgset_id_type = "symbol",
  bgset_description = "All protein-coding genes",
  enrichment_p_value_cutoff = 0.05,
  enrichment_p_value_adj = "BH",
  enrichment_q_value_cutoff = 0.2,
  enrichment_min_geneset_size = 10,
  enrichment_max_geneset_size = 500,
  enrichment_plot_num_terms = 20,
  enrichment_simplify_go = TRUE,
  subcellcomp_min_confidence = 3,
  subcellcomp_min_channels = 1,
  regulatory_min_resources = 1,
  fitness_max_score = -2,
  ppi_add_nodes = 30,
  ppi_string_min_score = 0.9,
  ppi_string_network_type = "functional",
  ppi_biogrid_min_evidence = 3,
  ppi_node_shadow = TRUE,
  ppi_show_drugs = TRUE,
  ppi_show_isolated_nodes = FALSE,
  show_ppi = TRUE,
  show_disease = TRUE,
  show_top_diseases_only = TRUE,
  show_cancer_hallmarks = TRUE,
  show_drug = TRUE,
  show_enrichment = TRUE,
  show_aberration = FALSE,
  show_coexpression = FALSE,
  show_ligand_receptor = FALSE,
  show_regulatory = FALSE,
  show_unknown_function = TRUE,
  show_prognostic = TRUE,
  show_subcell_comp = FALSE,
  show_synleth = FALSE,
  show_fitness = TRUE,
  show_complex = TRUE,
  show_domain = FALSE,
  ...
)
```

## Arguments

- query:

  character vector with gene/query identifiers (minimum 2, maximum 1000)

- oeDB:

  oncoEnrichR data repository object - as returned from
  [`load_db()`](https://sigven.github.io/oncoEnrichR/reference/load_db.md)

- query_id_type:

  character indicating source of query (one of "uniprot_acc",
  "symbol","entrezgene", or "ensembl_gene", "ensembl_mrna",
  "refseq_transcript_id", "ensembl_protein", "refseq_protein")

- ignore_id_err:

  logical indicating if analysis should continue when uknown query
  identifiers are encountered

- project_title:

  project title (title of report)

- project_owner:

  name of project owner

- project_description:

  project background information

- bgset:

  character vector with gene identifiers, used as reference/background
  for enrichment/over-representation analysis

- bgset_id_type:

  character indicating source of background ("uniprot_acc", "symbol",
  "entrezgene", "ensembl_gene", "ensembl_mrna", "refseq_transcript_id",
  "ensembl_protein", "refseq_protein"), default: "symbol"

- bgset_description:

  character indicating type of background (e.g. "All lipid-binding
  proteins (n = 200)")

- enrichment_p_value_cutoff:

  cutoff p-value for enrichment/ over-representation analysis (default:
  0.05)

- enrichment_p_value_adj:

  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
  "none" (clusterProfiler, default: "BH")

- enrichment_q_value_cutoff:

  cutoff q-value for enrichment analysis (clusterProfiler, default: 0.2)

- enrichment_min_geneset_size:

  minimal size of geneset annotated by term for testing in
  enrichment/over-representation analysis (clusterProfiler, default: 10)

- enrichment_max_geneset_size:

  maximal size of geneset annotated by term for testing in
  enrichment/over-representation analysis (clusterProfiler, default:
  500)

- enrichment_plot_num_terms:

  number of top enriched Gene Ontology terms (max) to show in enrichment
  barplot (default: 15)

- enrichment_simplify_go:

  remove highly similar GO terms in results from GO
  enrichment/over-representation analysis (default: TRUE)

- subcellcomp_min_confidence:

  minimum confidence level for subcellular compartment annotation in
  COMPARTMENTS (min = 3, max = 5, default: 3)

- subcellcomp_min_channels:

  minimum number of channels that support a subcellular compartment
  annotation in COMPARTMENTS (min = 1, max = 3, default: 1)

- regulatory_min_resources:

  minimum resource level for regulatory interactions (TF-target)
  retrieved from Collectri (min = 1, max = 6, default: 1)

- fitness_max_score:

  maximum loss-of-fitness score (scaled Bayes factor from BAGEL) for
  genes retrieved from DepMap/Project Score, default:-2

- ppi_add_nodes:

  number of nodes to add to target set when computing the
  protein-protein interaction network (STRING/BioGRID, default: 30)

- ppi_string_min_score:

  minimum score (between 0 and 1) for confidence of retrieved
  protein-protein interactions (STRING, default: 0.9)

- ppi_string_network_type:

  type of network to show for interactions in STRING ('functional' or
  'physical', default: 'functional')

- ppi_biogrid_min_evidence:

  minimum number of evidence items required for protein-protein
  interactions retrieved (BioGRID, default: 3)

- ppi_node_shadow:

  show shadow for nodes in the displayed PPI network (default: TRUE)

- ppi_show_drugs:

  logical indicating if targeted drugs (\>= phase 3) should be displayed
  in protein-protein interaction networks (default: TRUE)

- ppi_show_isolated_nodes:

  logical indicating if targets/nodes without any interactions should be
  displayed in the protein-protein interaction networks (default: FALSE)

- show_ppi:

  logical indicating if report should contain protein-protein
  interaction views of the query set (STRING and BioGRID, default: TRUE)

- show_disease:

  logical indicating if report should contain disease associations (Open
  Targets Platform, association_score \>= 0.05, support from at least
  two data types), and tumor suppressor/oncogene annotations ( default:
  TRUE)

- show_top_diseases_only:

  logical indicating if report should contain top (n = 20) disease
  associations only pr. query gene default: TRUE)

- show_cancer_hallmarks:

  logical indicating if report should contain annotations/evidence of
  cancer hallmarks per query gene (COSMIC/Open Targets Platform,
  default: TRUE

- show_drug:

  logical indicating if report should contain targeted cancer drug
  information (default: TRUE)

- show_enrichment:

  logical indicating if report should contain functional
  enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME,
  NetPath, WikiPathways, default: TRUE)

- show_aberration:

  logical indicating if report should contain TCGA aberration plots
  (amplifications/deletions, default: TRUE)

- show_coexpression:

  logical indicating if report should contain TCGA co-expression data
  (RNAseq) of query set with oncogenes/tumor suppressor genes (default:
  TRUE) \#param show_cell_tissue logical indicating if report should
  contain tissue-specificity and single cell-type specificity
  assessments (Human Protein Atlas) of target genes (default: FALSE)

- show_ligand_receptor:

  logical indicating if report should contain ligand-receptor
  interactions (CellChatDB, default: TRUE)

- show_regulatory:

  logical indicating if report should contain data on transcription
  factor (TF) - target interactions relevant for the query set
  (DoRothEA, default: TRUE)

- show_unknown_function:

  logical indicating if report should highlight target genes with
  unknown or poorly defined functions (GO/Uniprot KB/NCBI, default:
  TRUE)

- show_prognostic:

  logical indicating if mRNA-based (single-gene) prognostic associations
  to cancer types should be listed (Human Protein Atlas/TCGA, default:
  TRUE

- show_subcell_comp:

  logical indicating if report should provide subcellular compartment
  annotations (COMPARTMENTS, default: TRUE)

- show_synleth:

  logical indicating if report should list overlap with predicted
  synthetic lethality interactions (gene paralogs only, De Kegel et al.,
  Cell Systems, 2021). Default: TRUE

- show_fitness:

  logical indicating if report should provide fitness scores and target
  priority scores from CRISPR/Cas9 loss-of-fitness screens
  (DepMap/Project Score, default: TRUE)

- show_complex:

  logical indicating if report should provide target memberships in
  known protein complexes (ComplexPortal/Compleat/hu.MAP2/PDB/CORUM,
  default: TRUE)

- show_domain:

  logical indicating if report should provide target memberships in
  known protein domains (Pfam, default: TRUE)

- ...:

  arguments for Galaxy/web-based processing

## Value

An oncoEnrichR report list object, with two main elements, `data` and
`config`. `data` contains data that goes into each section of the output
reports, `config` contains all metadata for annotation resources used.
