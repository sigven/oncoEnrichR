## oncoEnrichR - Functional profiling of genesets in the context of cancer

### Contents

- [Overview](#overview)
- [News](#news)
- [Annotation resources](#annotation-resources)
- [Example report](#example-report)
- [Getting started](#getting-started)
- [Contact](#contact)

### Overview

<img align="center" width = "560" height = "700" src="docs/oncoenrichr2.gif">


__oncoEnrichR__ is an R package for functional interrogation of human genesets in the context of cancer.

The package is intended for exploratory analysis and prioritization of a gene list (referred to as **target set** below) from high-throughput cancer biology experiments. The tool queries a number of high-quality data resources in order to assemble useful gene annotations and analyses in an interactive report. The contents of the final report attempts to shed light on the following questions:

  * Which diseases/tumor types are associated with genes in the target set?
  * Which proteins in the target sets are druggable in diffferent cancer conditions (early and late clinical development phases)?
  * Which protein complexes are relevant for proteins in the target set?
  * Which subcellular compartments (nucleus, cytosol, plasma membrane etc) are dominant localizations for proteins in the target set?
  * Which protein-protein interactions are known within the target set? Are there interactions between members of the target set and other cancer-relevant proteins (e.g. proto-oncogenes, tumor-suppressors or predicted cancer drivers)? Which proteins constitute hubs in the protein-protein interaction network?
  * Are there specific pathways, biological processes or molecular functions that are enriched within the target set, as compared to a reference/background set?
  * Which members of the target set are frequently mutated in tumor sample cohorts (TCGA, SNVs/InDels, homozygous deletions, copy number amplifications)?
  * Which members of the target set are co-expressed (strong negative or positive correlations) with cancer-relevant genes (i.e. proto-oncogenes or tumor suppressors) in tumor sample cohorts (TCGA)?
  * Which members of the target set are associated with cellular loss-of-fitness in CRISPR/Cas9 whole-genome drop out screens of cancer cell lines (i.e. reduction of cell viability elicited by a gene inactivation)?

### News

### Annotation resources

Data harvested from the following resources form the backbone of _oncoEnrichR_:

* [Open Targets Platform](https://targetvalidation.org) - drug-target associations and disease-target associations
* [The Cancer Genome Atlas](https://portal.gdc.cancer.gov/) - gene aberration frequencies and co-expression patterns in > 10,000 tumor samples
* [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp) - collection of annotated (e.g. towards pathways) genesets for enrichment/overrepresentation analysis. This includes genesets from [Gene Ontology](http://geneontology.org/), [Reactome](https://reactome.org/), [KEGG](https://www.genome.jp/kegg/pathway.html), [WikiPathways](https://www.wikipathways.org/index.php/WikiPathways), [BIOCARTA](https://maayanlab.cloud/Harmonizome/dataset/Biocarta+Pathways), as well as curated [immunologic](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7) and [cancer-specific](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C6) signatures.
* [OmnipathDB](http://omnipathdb.org) - literature-cuated mammalian signaling pathways
* [STRING](https://string-db.org) - protein-protein interaction database
* [CORUM](https://mips.helmholtz-muenchen.de/corum/) - protein complex database
* [ComPPI](http://comppi.linkgroup.hu/) - subcellular compartment database
* [Project Score](https://score.depmap.sanger.ac.uk) - Database on the effects on cancer cell line viability elicited by CRISPR-Cas9 mediated gene activation



### Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4050882.svg)](https://doi.org/10.5281/zenodo.4050882)


### Getting started

#### Installation
1. `install.packages('devtools')`
2. `devtools::install_github('sigven/oncoEnrichR')`

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)

#### Usage

_oncoEnrichR_ works in two basic steps through the following two methods:

   __1.__ `onco_enrich`
  * Takes an input list of human gene/protein identifiers and performs various types of annotations, frequency calculations, and enrichment analyses. Currently acceptable input formats are HGNC symbols, UniProt accessions, Ensembl and Entrez gene identifiers. Returns a *list object* with all contents of the analyses performed. Arguments and default values:

	  ```r
	  onco_enrich(
	    query,
	    query_source = "symbol",
	    ignore_unknown = FALSE,
	    project_title = "Project title",
	    project_owner = "Project owner",
	    project_description = "Project description",
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
	    show_ppi = TRUE,
	    show_drugs_in_ppi = FALSE,
	    show_disease = TRUE,
	    show_drug = TRUE,
	    show_enrichment = TRUE,
	    show_tcga_aberration = TRUE,
	    show_tcga_coexpression = TRUE,
	    show_subcell_comp = TRUE,
	    show_crispr_lof = TRUE,
	    show_complex = TRUE
	  )
	  ```


	Argument      |Description
	  ------------- |----------------
	  ```query```     |     character vector with gene/query identifiers
	  ```query_source```     |     character indicating source of query (one of 'uniprot_acc', 'symbol', 'entrezgene', or 'ensembl_gene_id')
	  ```ignore_unknown```     |     logical indicating if analysis should continue when uknown query identifiers are encountered
	  ```project_title```     |     project title (report title)
	  ```project_owner```     |     project owner (e.g. lab/PI)
	  ```project_description```     |     brief description of project, how target list was derived
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
	  ```ppi_score_threshold```     |     minimum significance score (0-1000) for protein-protein interactions to be included in network (STRING)
	  ```show_ppi```     |     logical indicating if report should contain protein-protein interaction data (STRING)
	  ```show_drugs_in_ppi```     |     logical indicating if targeted drugs (> phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
	  ```show_disease```     |     logical indicating if report should contain disease associations (Open Targets Platform)
	  ```show_drug``` | logical indicating if report should contain cancer drugs targeted towards proteins in the query list (early and late development phase, from Open Targets Platform)
	  ```show_enrichment```     |     logical indicating if report should contain functional enrichment/over-representation analysis (MSigDB, GO, KEGG, REACTOME etc.)
	  ```show_tcga_aberration```     |     logical indicating if report should contain TCGA aberration plots (amplifications/deletions)
	  ```show_tcga_coexpression```     |     logical indicating if report should contain TCGA co-expression data (RNAseq) of queryset with oncogenes/tumor suppressor genes
	  ```show_subcell_comp```     |     logical indicating if report should list subcellular compartment annotations (ComPPI)
	  ```show_crispr_lof```     |     logical indicating if report should list results from CRISPR/Cas9 loss-of-fitness screens (Project Score)
	  ```show_complex```     |     logical indicating if report should list proteins in known protein complexes (CORUM)



  __2.__ `write`

* takes the contents of the report object retrieved in _1)_ to assemble a structured and interactive _oncoEnrichR_ HTML report

#### Example run

A target list of _n = 134_ high-confidence interacting proteins with the c-MYC oncoprotein were previously identified through BioID protein proximity assay in standard cell culture and in tumor xenografts ([Dingar et al., J Proteomics, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25452129)). We ran this target list through the _oncoEnrichR_ analysis workflow using the following configurations for the `onco_enrich` method:

  * `ignore_unknown = TRUE`
  * `query_source = "symbol"`
  * `project_title = "cMYC_BioID_screen"`
  * `project_owner = "Raught et al."`
  * `show_drugs_in_ppi = TRUE`
  * `simplify_go = TRUE`

 and produced the [following HTML report with results](https://doi.org/10.5281/zenodo.4022174).

 Below are R commands provided to reproduce the example output ("LOCAL_FOLDER") is replaced with a directory on your local computer:

 * `myc_data <- read.csv(system.file("extdata","myc_data.csv", package = "oncoEnrichR"), stringsAsFactors = F)`
 * `myc_report <- oncoEnrichR::onco_enrich(myc_data$symbol, query_source = "symbol", ignore_unknown = T, project_title = "cMYC_BioID_screen", project_owner = "Raught et al.", show_drugs_in_ppi = T, simplify_go = T)`
 * `oncoEnrichR::write(myc_report, "LOCAL_FOLDER", "cmyc_example")`

### Contact

sigven AT ifi.uio.no

### Funding & collaboration

OncoEnrichR is supported by the [Centre for Cancer Cell Reprogramming](https://www.med.uio.no/cancell/english/) at the [University of Oslo](https://www.uio.no)/[Oslo University Hospital](https://radium.no), and [Elixir Norway (Oslo node)](https://elixir.no/organization/organisation/elixir-uio).

<br>
<br>

<p float="left">
  <a href="https://www.med.uio.no/cancell/english/">
     <img src="can-cell.png" width="150" >
  </a>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://elixir.no/organization/organisation/elixir-uio">
     <img src="elixir_norway.png" width="200" />
  </a>
</p>
