
# Version 1.0.7

## Added

* Data updates
  * WikiPathways (20211110)
  * Open Targets Platform (2021.11)
  * EFO (v3.36.0)
  * TCGA (v31, 2021-10-29)
  * CancerMine (v40)
  * Human Protein Atlas (v21)
  * UniProt KB (2021_04)
  * PFAM domains (release 35, November 2021)
  
## Changed

- Management of annotation databases have been re-designed, reducing the size 
and installation of *oncoEnrichR* significantly. The re-design results in
the addition of an initial step in the analysis; **database loading**.


# Version 1.0.6

## Added

* Data updates 
  * KEGG (20211013)
  * WikiPathways (20211010)
  * Open Targets Platform (2021.09)
  * EFO (v3.35.0), DiseaseOntology (v2021-10-11)
  * TCGA (v30, 2021-09-23)
  * Protein complex annotations - additions from [Compleat](https://fgr.hms.harvard.edu/compleat), [ComplexPortal](https://www.ebi.ac.uk/complexportal/home), [PDB](https://www.rcsb.org/), [hu.MAP2](http://humap2.proteincomplexes.org/)

* New functionality
  * Report modules
     * **Regulatory interactions** - ([DoRothEA](https://saezlab.github.io/dorothea/))
        * View overlap between members of queryset with previously established regulatory interactions (TF-target relationships) from the DoRothEA resource
     * **Ligand-receptor interactions** -  ([CellChatDB](http://www.cellchat.org/))
        * Explore potential ligand-receptor interactions in the query set, as found in the CellChat database

  * Report styling options
     * Choose your own Bootswatch theme for HTML report - option `html_report_theme`
     * Choose where the table of contents are placed in the HTML report (floating-left or static at the top) - option `html_floating_toc`

  * Arguments to `oncoEnichR::onco_enrich()`:
     * `html_floating_toc` - logical, if set to FALSE, table of contents in HTML report will be placed on top of the main document
     * `html_report_theme` - character, choose between different Bootswatch themes for style in HTML report
     * `show_regulatory` - logical, show regulatory interactions module (DoRothEA)
     * `min_confidence_reg_interaction` - character, minimum confidence of regulatory interactions included from DoRothEA ('A','B','C', or 'D')
     * `show_ligand_receptor` - logical, show ligand-receptor interaction module (CellChatDB)
     * `num_terms_enrichment_plot` - integer, number of enriched Gene Ontology terms (max) to show in enrichment barplots (module *Functional Enrichment*)

  * Default gene rankings
     * In the modules **Regulatory Interactions** and **Tissue and cell type enrichment**, interactions/targets are ranked according to a quantitative cancer association score *cancer_max_rank* (maximum gene-cancer association rank (Open Targets Platform) across primary sites)

## Fixed
  * Report *disclaimer* occurring at the bottom of the report is no longer missing from web-based analysis (oncotools.elixir.no) - file `_site.yml` in `inst/templates`

# Version 1.0.5

## Added

* Data updates: KEGG (20210916), WikiPathways (20210910), CancerMine (v39),
EFO (v3.34.0)

# Version 1.0.4

## Added

* Data updates: Open Targets Platform (21.06), CancerMine (v37)
* Software upgrades: multiple R packages, including clusterProfiler (4.0.2)

# Version 1.0.3

## Added

* Functionality for Galaxy-specific HTML output (non-selfcontained HTML)
* Data updates: WikiPathways (20210610), KEGG (20210610), CancerMine (v36, 20210611)
* Web-based interface to oncoEnrichR available: [oncotools.elixir.no](https://oncotools.elixir.no)

## Fixed

* Minor typos in HTML output
* Ambiguous Entrez gene identifiers (towards gene symbols)

# Version 1.0.0

## New

* Added prognostic Z-scores also for methylation features ([Smith et al., _bioRxiv_, 2021](https://www.biorxiv.org/content/10.1101/2021.06.01.446243v1))
* To keep the file size of output Excel/HTML at managable levels, maximum size of queryset is reduced to **n = 600 genes**

## Fixed

* Fixed assymetrical scales of prognostic Z-score plots
* To decrease the size of the HTML report, limited datatable output of recurrent somatic SNVs/InDels and sCNA frequencies to top 2,500 (full set is kept in Excel output)

# Version 0.9.7 (June 4th 2021) 

## New

* More prognostic associations - harvested from [Smith et al., _bioRxiv_, 2021](https://www.biorxiv.org/content/10.1101/2021.06.01.446243v1)

## Fixed

* DESCRIPTION: Added `TissueEnrich` and `assertable` dependency

## Data upgrade

* Target-disease associations from Open Targets Platform (v2021.04)
* WikiPathways (20210510)
* Disease Ontology (20210603)
* KEGG (20210531)

