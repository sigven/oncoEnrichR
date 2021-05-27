## oncoEnrichR - cancer-dedicated gene set interpretation

### Contents

- [Overview](#overview)
- [News](#news)
- [Annotation resources](#annotation-resources)
- [Example report](#example-report)
- [Getting started](#getting-started)
- [Contact](#contact)

### Overview

__oncoEnrichR__ is an R package for functional interrogation of human genesets in the context of cancer. The package is intended for exploratory analysis and prioritization of a gene list (referred to here as **query set**) from high-throughput cancer biology experiments, e.g. genetic screens (siRNA/CRISPR), protein proximity labeling, or transcriptomics (differential expression). The tool queries a number of high-quality data resources in order to assemble useful gene annotations and analyses in an interactive report (examples from the report shown below).

<br>

<img align="center" border = 2, src="man/figures/oncoenrichr_overview.png"/>

<br><br>

The contents of the final report attempts to shed light on the following questions:

  * Which diseases/tumor types are associated with genes in the query set, and to what extent?
  * Which genes in the query set are attributed with cancer hallmark evidence?
  * Which proteins in the query sets are druggable in diffferent cancer conditions (early and late clinical development phases)? For other proteins in the query set, what is their likelihood of being druggable?
  * Which protein complexes are relevant for proteins in the query set?
  * Which subcellular compartments (nucleus, cytosol, plasma membrane etc) are dominant localizations for proteins in the query set?
  * Are specific tissues or cell types enriched in the query set, considering tissue/cell-type specific expression patterns of target genes?
  * Which protein-protein interactions are known within the query set? Are there interactions between members of the query set and other cancer-relevant proteins (e.g. proto-oncogenes, tumor-suppressors or predicted cancer drivers)? Which proteins constitute hubs in the protein-protein interaction network?
  * Are there specific pathways, biological processes or molecular functions that are enriched within the query set, as compared to a reference/background set?
  * Which members of the query set are frequently mutated in tumor sample cohorts (TCGA, SNVs/InDels, homozygous deletions, copy number amplifications)?
  * Which members of the query set are co-expressed (strong negative or positive correlations) with cancer-relevant genes (i.e. proto-oncogenes or tumor suppressors) in tumor sample cohorts (TCGA)?
  * Which members of the query set are associated with better/worse survival in different cancers, considering high or low gene expression levels in tumors?
  * Which members of the query set are associated with cellular loss-of-fitness in CRISPR/Cas9 whole-genome drop out screens of cancer cell lines (i.e. reduction of cell viability elicited by a gene inactivation)? Which targets are prioritized therapeutic targets, considering fitness effects and genomic biomarkers in combination?

### News
* May 12th 2021: **0.9.6 release**
  * Added cancer hallmarks evidence
  * Fixed minor bugs
  * Drug retrieval from Open Targets Platform v2021.04


### Annotation resources

Data harvested from the following resources form the backbone of _oncoEnrichR_:

* [Open Targets Platform](https://targetvalidation.org) - human drug-target associations and comprehensive disease-target associations
* [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - literature-mined database of drivers, oncogenes and tumor suppressors in cancer
* [The Cancer Genome Atlas](https://portal.gdc.cancer.gov/) - gene aberration frequencies and co-expression patterns in ~9,500 primary tumor samples
* [The Human Protein Atlas](https://www.proteinatlas.org) - expression data for healthy human tissues ([GTex](https://gtexportal.org/home/))/cell types, and prognostic gene expression associations in cancer
* [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp) - collection of annotated (e.g. towards pathways) genesets for enrichment/overrepresentation analysis. This includes genesets from [Gene Ontology](http://geneontology.org/), [Reactome](https://reactome.org/), [KEGG](https://www.genome.jp/kegg/pathway.html), [WikiPathways](https://www.wikipathways.org/index.php/WikiPathways), [BIOCARTA](https://maayanlab.cloud/Harmonizome/dataset/Biocarta+Pathways), as well as curated [immunologic](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7) and [cancer-specific](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C6) signatures.
* [NetPath](http://www.netpath.org) - signaling transduction pathways
* [STRING](https://string-db.org) - protein-protein interaction database
* [CORUM](https://mips.helmholtz-muenchen.de/corum/) - protein complex database
* [ComPPI](http://comppi.linkgroup.hu/) - subcellular compartment database
* [Project Score](https://score.depmap.sanger.ac.uk) - Database on the effects on cancer cell line viability elicited by CRISPR-Cas9 mediated gene activation


### Example report

<a href="https://doi.org/10.5281/zenodo.4751107"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4751107.svg" alt="DOI"></a>

<br>

<img align="center" width = "560" height = "750" border = 2, src="man/figures/oncoenrichr_slideshow2.gif"/>


### Getting started

#### Installation
1. `install.packages('devtools')`
2. `devtools::install_github('sigven/oncoEnrichR')`
3. `library(oncoEnrichR)`

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)

#### Usage

_oncoEnrichR_ performs its operations through the following procedures/methods:

   __1.__ `oncoEnrichR::onco_enrich()`

   * Consists of two main processing steps:

	   1) Takes an input/query list of human gene/protein identifiers (e.g. UniProt accession, RefSeq/Ensembl transcript identifer etc.) as input and conducts uniform identifier conversion

	   2) Performs extensive annotation, enrichment and membership analyses of the query set against underlying data sources on cancer-relevant properties of human genes and their interrelationships.

  * Technically, the method returns a *list object* with all contents of the analyses performed. The specific arguments/options and default values are outlined below:

	  ```r
	  onco_enrich(
	    query,
	    query_id_type = "symbol",
	    ignore_id_err = TRUE,
	    project_title = "Project title",
	    project_owner = "Project owner",
	    project_description = "Project description",
	    bgset = NULL,
	    bgset_id_type = "symbol",
	    bgset_description = "All protein-coding genes",
	    p_value_cutoff_enrichment = 0.05,
	    p_value_adjustment_method = "BH",
	    q_value_cutoff_enrichment = 0.2,
	    min_geneset_size = 10,
	    max_geneset_size = 500,
	    min_subcellcomp_confidence = 1,
	    subcellcomp_show_cytosol = FALSE,
	    simplify_go = TRUE,
	    ppi_add_nodes = 50,
	    ppi_score_threshold = 900,
	    show_ppi = TRUE,
	    show_drugs_in_ppi = TRUE,
	    show_disease = TRUE,
	    show_top_diseases_only = TRUE,
	    show_cancer_hallmarks = TRUE,
	    show_drug = TRUE,
	    show_enrichment = TRUE,
	    show_tcga_aberration = TRUE,
	    show_tcga_coexpression = TRUE,
	    show_subcell_comp = TRUE,
	    show_crispr_lof = TRUE,
	    show_cell_tissue = FALSE,
	    show_prognostic_cancer_assoc = TRUE,
	    show_complex = TRUE)
	  ```


	Argument      |Description
	  ------------- |----------------
	  ```query```     |     character vector with gene/query identifiers
	  ```query_id_type```     |     character indicating type of identifier used for query set (one of "_uniprot_acc_", "_symbol_", "_entrezgene_", "_ensembl_gene_", "_refseq_mrna_", "_refseq_protein_", "_ensembl_protein_", or "_ensembl_mrna_")
	  ```ignore_id_err```     |     logical indicating if analysis should continue when erroneous/unmatched query identifiers are encountered in query or background gene set
	  ```project_title```     |     project title (report title)
	  ```project_owner```     |     project owner (e.g. lab/PI)
	  ```project_description```     |     brief description of project, how target list was derived
	  ```bgset```     |     character vector with gene identifiers, used as reference/background for enrichment/over-representation analysis
	  ```bgset_id_type```     |     character indicating type of identifier used for background set (one of "_uniprot_acc_", "_symbol_", "_entrezgene_", "_ensembl_gene_", "_refseq_mrna_", "_refseq_protein_", "_ensembl_protein_", or "_ensembl_mrna_")
	  ```bgset_description```     |     character with description of background gene set (e.g. 'All lipid-binding proteins (n = 200)')
	  ```p_value_cutoff_enrichment```     |     cutoff p-value for enrichment/over-representation analysis
	  ```p_value_adjustment_method```     |     one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
	  ```q_value_cutoff_enrichment```     |     cutoff q-value for enrichment/over-representation analysis
	  ```min_geneset_size```     |     minimal size of geneset annotated by term for testing in enrichment/over-representation analysis
	  ```max_geneset_size```     |     maximal size of geneset annotated by term for testing in enrichment/over-representation analysis
	  ```min_subcellcomp_confidence```|     minimum confidence level of subcellular compartment annotations (range from 1 to 6, 6 is strongest)
	  ```subcellcomp_show_cytosol``` | logical indicating if subcellular heatmap should highlight cytosol as a subcellular protein location or not
	  ```simplify_go```     |     remove highly similar GO terms in results from GO enrichment/over-representation analysis (recommended)
	  ```ppi_add_nodes```     |     number of neighbouring nodes to add to query set when computing the protein-protein interaction network (STRING)
	  ```ppi_score_threshold```     |     minimum significance score (0-1000) for protein-protein interactions to be included in network (STRING)
	  ```show_ppi```     |     logical indicating if report should contain protein-protein interaction data of query set and their closely interacting partners (STRING), local network communities, and rank of query set based on network centrality
	  ```show_drugs_in_ppi```     |     logical indicating if targeted drugs (>= phase 3) should be displayed in protein-protein interaction network (Open Targets Platform)
	  ```show_disease```     |     logical indicating if report should contain ranked associations to cancer phenotypes (overall), as well as tumor-type specific rankings (association score >= 0.4, minimum number of sources contributing to association >= 2, (Open Targets Platform))
	  ```show_top_diseases_only```| logical indicating if report should only show top (20) cancer phenotypes/disease associations from Open Targets Platform
	  ```show_cancer_hallmarks```| logical indicating if report should show evidence that relates query member genes to cancer hallmarks
	  ```show_drug``` | logical indicating if report should contain cancer drugs targeted towards proteins in the query list (early and late development phase) and tractability/druggability data for all query entries, from Open Targets Platform)
	  ```show_enrichment```     |     logical indicating if report should perform and list functional enrichment/over-representation analysis of query set (MSigDB, GO, KEGG, REACTOME, WikiPathways, NetPath)
	  ```show_tcga_aberration```     |     logical indicating if report should contain TCGA aberration plots (amplifications/deletions, SNVs/InDels (oncoplots))
	  ```show_tcga_coexpression```     |     logical indicating if report should list oncogenes/tumor suppressor genes that significantly correlate with entries in query set in terms of expression (across TCGA cohorts)
	  ```show_subcell_comp```     |     logical indicating if report should list subcellular compartment annotations (ComPPI)
	  ```show_crispr_lof```     |     logical indicating if report should list results from CRISPR/Cas9 loss-of-fitness screens and associated target priority scores (Project Score)
	  ```show_cell_tissue```    | 	logical indicating if report should list results from tissue (GTex)- and cell-type (HPA) specific gene expression patterns in query set
	  ```show_prognostic_cancer_assoc```  |    logical indicating if report should list results from significant associations between gene expression and survival (Human Protein Atlas - Pathology Atlas)
	  ```show_complex```     |     logical indicating if report should show membership of target proteins in known protein complexes (CORUM)



  __2.__ `oncoEnrichR::write()`

  * Consists of two main processing steps:

    1) Transform the contents of the analyses returned by _oncoEnrichR::onco_enrich()_ into various visualizations and interactive tables

    2) Assemble and write the final analysis report through

	   * A) a structured and interactive _oncoEnrichR_ HTML report, OR
	   * B) a multisheet Excel workbook with results from the annotations and analyses provided
     by _oncoEnrichR_

#### Example run

A target list of _n = 134_ high-confidence interacting proteins with the c-MYC oncoprotein were previously identified through BioID protein proximity assay in standard cell culture and in tumor xenografts ([Dingar et al., J Proteomics, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25452129)). We ran this target list through the _oncoEnrichR_ analysis workflow using the following configurations for the `onco_enrich` method:

  * `project_title = "cMYC_BioID_screen"`
  * `project_owner = "Raught et al."`

 and produced the [following HTML report with results](https://doi.org/10.5281/zenodo.4751107).

 Below are R commands provided to reproduce the example output. __NOTE__: Replace "LOCAL_FOLDER" with a directory on your local computer:

 * `library(oncoEnrichR)`
 * `myc_interact_targets <- read.csv(system.file("extdata","myc_data.csv", package = "oncoEnrichR"), stringsAsFactors = F)`
 * `myc_report <- oncoEnrichR::onco_enrich(query = myc_interact_targets$symbol, show_cell_tissue = T, project_title = "cMYC_BioID_screen", project_owner = "Raught et al.")`
 * `oncoEnrichR::write(report = myc_report, file = "LOCAL_FOLDER/myc_report_oncoenrichr.html", format = "html")`
 * `oncoEnrichR::write(report = myc_report, file = "LOCAL_FOLDER/myc_report_oncoenrichr.xlsx", format = "excel")`

### Contact

sigven AT ifi.uio.no

### Funding & collaboration

OncoEnrichR is supported by the [Centre for Cancer Cell Reprogramming](https://www.med.uio.no/cancell/english/) at the [University of Oslo](https://www.uio.no)/[Oslo University Hospital](https://radium.no), and [Elixir Norway (Oslo node)](https://elixir.no/organization/organisation/elixir-uio).

<br>
<br>

<p float="left">
  <a href="https://www.med.uio.no/cancell/english/">
     <img src="man/figures/can-cell.png" width="150" >
  </a>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://elixir.no/organization/organisation/elixir-uio">
     <img src="man/figures/elixir_norway.png" width="200" />
  </a>
</p>
