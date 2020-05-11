## oncoEnrichR

<img align="center" width = "560" height = "700" src="docs/oncoenrichr2.gif">


__oncoEnrichR__ is an R package for functional interrogation of human genesets in the context of cancer.

The package is intended for exploratory analysis and prioritization of a gene list (referred to as **target set** below) from high-throughput cancer biology experiments. The tool queries a number of high-quality data resources in order to assemble useful gene annotations and analyses in an interactive report. The following resources are currently included:

* [Open Targets Platform](https://targetvalidation.org)
* [The Cancer Genome Atlas](https://portal.gdc.cancer.gov/)
* [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp)
* [OmnipathDB](http://omnipathdb.org)
* [STRING](https://string-db.org)
* [CORUM](https://mips.helmholtz-muenchen.de/corum/)
* [ComPPI](http://comppi.linkgroup.hu/)
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

### Installation
3. `install.packages('devtools')`
3. `devtools::install_github('sigven/oncoEnrichR')`

[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)

### Usage

_oncoEnrichR_ works in two basic steps through the following two methods:

   __1.__ `onco_enrich`
  * Takes an input list of human gene/protein identifiers and performs various types of annotations, frequency calculations, and enrichment analyses. Currently acceptable input formats are HGNC symbols, UniProt accessions, Ensembl and Entrez gene identifiers. Returns a *list object* with all contents of the analyses performed

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



  __2.__ `write`

* takes the contents of the report object retrieved in _1)_ to assemble a structured and interactive _oncoEnrichR_ HTML report

### Example output

A target list of _n = 134_ high-confidence interacting proteins with the c-MYC oncoprotein were previously identified through BioID protein proximity assay in standard cell culture and in tumor xenografts ([Dingar et al., J Proteomics, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25452129)).

We ran this target list, as found in [this Excel file](myc_candidates.xlsx), through the _oncoEnrichR_ analysis workflow using the following configurations for the `onco_enrich` method:

  * `ignore_unknown = TRUE`
  * `query_source = "symbol"`
  * `p_title = "cMYC_BioID_screen"`
  * `p_owner = "Raught et al."`
  * `show_drugs_in_ppi = TRUE`
  * `simplify_go = TRUE`

 and produced the [following HTML report with results](https://folk.uio.no/sigven/cmyc_example.html).

 R commands to reproduce the example output (requires that the [openxlsx package](https://github.com/ycphs/openxlsx) is installed on your local computer), and that "_LOCAL_FOLDER_" is replaced with a directory on your local computer:

 * `myc_data_raw <- openxlsx::read.xlsx("myc_candidates.xlsx",colNames = F)`
 * `myc_report <- oncoEnrichR::onco_enrich(myc_data_raw$X1, query_source = "symbol", ignore_unknown = T, p_title = "cMYC_BioID_screen", p_owner = "Raught et al.", show_drugs_in_ppi = T, simplify_go = T)`
 * `oncoEnrichR::write(myc_report, "LOCAL_FOLDER", "cmyc_example")`
