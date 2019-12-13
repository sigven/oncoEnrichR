# `onco_enrich`: Function that interrogates and analyzes a list of human protein-coding genes for cancer relevance

## Description


 Function that interrogates and analyzes a list of human protein-coding genes for cancer relevance


## Usage

```r
onco_enrich(
  query,
  query_source = "symbol",
  ignore_unknown = FALSE,
  p_title = "Project Title",
  p_owner = "Project Owner",
  background_fname = NULL,
  background_enrichment = NULL,
  background_enrichment_source = "symbol",
  background_enrichment_description = "All protein-coding genes",
  p_value_cutoff_enrichment = 0.05,
  p_value_adjustment_method = "BH",
  q_value_cutoff_enrichment = 0.2,
  min_geneset_size = 10,
  max_geneset_size = 500,
  simplify_go = F,
  ppi_add_nodes = 50,
  ppi_score_threshold = 900,
  show_ppi = T,
  show_drugs_in_ppi = F,
  show_disease = T,
  show_enrichment = T,
  show_tcga_aberration = T,
  show_tcga_coexpression = T,
  show_subcell_comp = T,
  show_crispr_lof = T,
  show_complex = T
)
```


## Arguments

Argument      |Description
------------- |----------------
```query```     |     character vector with gene/query identifiers
```query_source```     |     character indicating source of query (one of 'uniprot_acc', 'symbol', 'entrezgene', or 'ensembl_gene_id')
```ignore_unknown```     |     logical indicating if analysis should continue when uknown query identifiers are encountered
```p_title```     |     title of report
```p_owner```     |     name of project owner
```background_fname```     |     filename for simple text file with project background information, one line per background item
```background_enrichment```     |     character vector with gene identifiers, used as reference/background for enrichment/over-representation analysis
```background_enrichment_source```     |     character indicating source of background ('uniprot_acc','symbol','entrezgene','ensembl_gene_id')
```background_enrichment_description```     |     character indicating type of background (e.g. 'All lipid-binding proteins (n = 200)')
```p_value_cutoff_enrichment```     |     cutoff p-value for enrichment/over-representation analysis
```p_value_adjustment_method```     |     one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
```q_value_cutoff_enrichment```     |     cutoff q-value for enrichment analysis
```min_geneset_size```     |     minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
```max_geneset_size```     |     maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
```simplify_go```     |     remove highly similar GO terms in results from GO enrichment/over-representation analysis
```ppi_add_nodes```     |     number of nodes to add to query set when computing the protein-protein interaction network (STRING)
```ppi_score_threshold```     |     minimum score (0-1000) for retrieval of protein-protein interactions (STRING)
```show_ppi```     |     logical indicating if report should contain protein-protein interaction data (STRING)
```show_drugs_in_ppi```     |     logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
```show_disease```     |     logical indicating if report should contain disease associations (Open Targets Platform)
```show_enrichment```     |     logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME etc.)
```show_tcga_aberration```     |     logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
```show_tcga_coexpression```     |     logical indicating if report should contain TCGA co-expression data (RNAseq) of queryset with oncogenes/tumor suppressor genes
```show_subcell_comp```     |     logical indicating if report should list subcellular compartment annotations (ComPPI)
```show_crispr_lof```     |     logical indicating if report should list results from CRISPR/Cas9 loss-of-fitness screens (Project Score)
```show_complex```     |     logical indicating if report should list proteins in known protein complexes (CORUM)

