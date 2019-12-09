## oncoEnrichR

<img align="center" width = "560" height = "700" src="docs/oncoenrichr2.gif">


__oncoEnrichR__ is an R package for functional interrogation of human genesets in the context of cancer.

The package is intended for exploratory analysis and prioritization of a gene list (referred to as **target set** below) from high-throughput cancer biology experiments. The tool queries a number of high-quality resources, in order to assemble useful gene annotations and analyses in an interactive report. The most prominent resources include:

* [Open Targets Platform](https://targetvalidation.org)
* [The Cancer Genome Atlas](https://portal.gdc.cancer.gov/)
* [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp)
* [OmnipathDB](http://omnipathdb.org)
* [STRING](https://string-db.org)
* [Project Score](https://score.depmap.sanger.ac.uk)

 The contents of the final report attempts to shed light on the following questions:

  * Which diseases/tumor types are associated with genes in the target set? Are there druggable proteins in the target set?
  * Which protein complexes are relevant for proteins in the target set?
  * Which subcellular compartments (nucleus, cytosol, plasma membrane etc) are dominant for proteins in the target set?
  * Which protein-protein interactions are known within the target set? Are there interactions between members of the target set and other cancer-relevant proteins (e.g. proto-oncogenes, tumor-suppressors or predicted cancer drivers)? Which proteins constitute hubs in the protein-protein interaction network?
  * Are there specific pathways or functional properties that are enriched within the target set, as compared to a reference/background set?
  * Which members of the target set are frequently mutated in tumor sample cohorts (TCGA, homozygous deletions/amplifications)?
  * Which members of the target set are co-expressed (strong negative or positive correlations) with proto-oncogenes or tumor suppressors in tumor sample cohorts (TCGA, RNAseq)?
  * Which members of the target set are associated with cellular loss-of-fitness in CRISPR/Cas9 whole-genome drop out screens of cancer cell lines (i.e. reduction of cell viability elicited by a gene inactivation)?

**To appear soon**: A more detailed documentation of the underlying routines/analyses performed in _oncoEnrichR_. Example use case and analysis.

### Installation

`devtools::install_github('sigven/oncoEnrichR')`

### Usage

_oncoEnrichR_ works in two basic steps through the following two methods:

   __1.__ `onco_enrich`
  * Takes an input list of human gene/protein identifiers and performs various types of annotations, frequency calculations, and enrichment analyses. Currently acceptable input formats are HGNC symbols, UniProt accessions, Ensembl and Entrez gene identifiers. Returns a *list object* with all contents of the analyses performed

  __2.__ `write`

* takes the contents of the report object retrieved in _1)_ to assemble a structured and interactive _oncoEnrichR_ HTML report

### Example output

A target list of _n = 134_ high-confidence interacting proteins with the c-MYC oncoprotein were previously identified through BioID protein proximity assay in standard cell culture and in tumor xenografts [Dingar et al., J Proteomics, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25452129).

We ran this [list](myc_candidates.xlsx) through the _oncoEnrichR_ analysis workflow using the following configurations for the `onco_enrich` method:

  * `ignore_unknown = TRUE`
  * `query_source = "symbol"`
  * `p_title = "cMYC_BioID_screen"`
  * `p_owner = "Raught et al."`
  * `show_drugs_in_ppi = TRUE`
  * `simplify_go = TRUE`

 and produced the [following HTML report with results](https://folk.uio.no/sigven/cmyc_example.html)
