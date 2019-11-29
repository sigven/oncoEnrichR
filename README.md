## oncoEnrichR

In short: An R package for functional interrogation of human genesets in the context of cancer.

The **oncoEnrichR** package is intended for exploratory analysis and prioritization of a gene list (referred to as **target set** below) from high-throughput cancer biology experiments (e.g. siRNA knowndown, proteomics etc). The tool queries a number of high-quality resources (i.e. [Open Targets Platform](https://targetvalidation.org), [The Cancer Genome Atlas](https://portal.gdc.cancer.gov/), [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp), [OmnipathDB](http://omnipathdb.org), and [STRING](https://string-db.org)) in order to assemble useful gene annotations and analyses in an interactive report. The contents of the final report attempts to shed light on the following questions:

  * Which diseases/tumor types are associated with genes in the target set? Are there druggable proteins in the target set?
  * Which protein complexes are relevant for proteins in the target set?
  * Which subcellular compartments (nucleus, cytosol, plasma membrane etc) are dominant for proteins in the target set?
  * Which protein-protein interactions are known within the target set? Are there interactions between members of the target set and other cancer-relevant proteins (e.g. proto-oncogenes, tumor-suppressors or predicted cancer drivers)? Which proteins constitute hubs in the protein-protein interaction network?
  * Are there specific pathways or functional properties that are enriched within the target set, as compared to a reference/background set?
  * Which members of the target set are frequently mutated in tumor sample cohorts (TCGA, homozygous deletions/amplifications)?
  * Which members of the target set are co-expressed (strong negative or positive correlations) with proto-oncogenes or tumor suppressors in tumor sample cohorts (TCGA, RNAseq)?

**To appear soon**: A more detailed documentation of the underlying routines/analyses performed in _oncoEnrichR_. Example use case and analysis.

### Installation

`devtools::install_github('sigven/oncoEnrichR')`

### Usage

_oncoEnrichR_ works in two basic steps through the following two methods:

   __1.__ `onco_enrich`
  * Takes an input list of human gene/protein identifiers and performs various types of annotations, frequency calculations, and enrichment analyses. Currently acceptable input formats are HGNC symbols, UniProt accessions, Ensembl and Entrez gene identifiers. Returns a *list object* with all contents of the analyses performed

  __2.__ `write`

* takes the contents of the report object retrieved in _1)_ to assemble a structured and interactive _oncoEnrichR_ HTML report
