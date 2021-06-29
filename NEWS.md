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

## Software upgrade

* Bioconductor 3.13
* Logging now performed with [log4r](https://github.com/johnmyleswhite/log4r)

# Version 0.9.6 (May 12th 2021) 

## Added

* Cancer hallmarks evidence

## Fixed

* Minor bugs (background query specificiation)

## Data upgrade

* Drug retrieval from Open Targets Platform (v2021.04)

# Version 0.9.4 (April 25th 2021)

##  Fixed 

* Cleaned Excel output

# Version 0.9.3 (April 24th 2021)

## Added

* Option to _onco_enrich_ that may ignore cytosol as a structure in subcelluar heatmap
* Interactive table with recurrent somatic variants (TCGA) overlapping query set in the _SNVs/InDels_ section
* Barplot for subcellular locations
* Barplot for most significant GO terms enriched in query set

## Fixed

* Bug in PPI's for query sets with few/none interactions

## Data upgrade

* MSigDB/GO/Reactome (MSigDB v7.4), UniProt (2021_02)
