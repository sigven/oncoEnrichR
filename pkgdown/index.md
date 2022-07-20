---
editor_options: 
  markdown: 
    wrap: 72
---

<br>

# oncoEnrichR

**oncoEnrichR** is an R package for functional interrogation of human
genesets in the context of cancer. It is primarily intended for
exploratory analysis, interpretation, and prioritization of long gene
lists, which represent a common output from many types of
high-throughput cancer biology screens.

**oncoEnrichR** can be used to interrogate results from e.g. genetic
screens (siRNA/CRISPR), protein proximity labeling, or transcriptomics
(differential expression). The tool queries a variety of high-quality
data resources in order to assemble useful gene annotations and analyses
in an interactive report (examples from the report shown below).

Web-based access to **oncoEnrichR** is available at <font size="4">
[**https://oncotools.elixir.no**](https://oncotools.elixir.no/tool_runner?tool_id=galaxy-ntnu.bioinfo.no%2Ftoolshed_nels%2Frepos%2Fradmilko%2Foncoenrichr%2Foncoenrichr_wrapper%2F1.2.1)</font>

<br>

<img src="img/oncoenrichr_overview.png" align="center" border="0"/>

<br><br>

## Questions adressed by oncoEnrichR

The contents of the analysis report provided by oncoEnrichR address the
following scientific questions for a given gene list:

-   Which diseases/tumor types are known to be associated with genes in
    the query set, and to what extent?
-   Which genes in the query set are attributed with cancer hallmark
    evidence?
-   Which proteins in the query sets are druggable in diffferent cancer
    conditions (early and late clinical development phases)? For other
    proteins in the query set, what is their likelihood of being
    druggable?
-   Which protein complexes involve proteins in the query set?
-   Which subcellular compartments (nucleus, cytosol, plasma membrane
    etc) are dominant localizations for proteins in the query set?
-   Are specific tissues or cell types enriched in the query set,
    considering tissue/cell-type specific expression patterns of target
    genes?
-   Which protein-protein interactions are known within the query set?
    Are there interactions between members of the query set and other
    cancer-relevant proteins (e.g. proto-oncogenes, tumor-suppressors or
    predicted cancer drivers)? Which proteins constitute hubs in the
    protein-protein interaction network?
-   Which known regulatory interactions (TF-target) are found within the
    query set, and what is their mode of regulation (repressive vs.
    stimulating)?
-   Are there occurrences of known ligand-receptor interactions within
    the query set?
-   Are there specific pathways, biological processes, or pre-defined
    molecular signatures that are enriched within the query set, as
    compared to a reference/background set?
-   Which members of the query set are frequently mutated in tumor
    sample cohorts (TCGA, SNVs/InDels, homozygous deletions, copy number
    amplifications)?
-   Which members of the query set are co-expressed (strong negative or
    positive correlations) with cancer-relevant genes (i.e.
    proto-oncogenes or tumor suppressors) in tumor sample cohorts
    (TCGA)?
-   Which members of the query set are associated with better/worse
    survival in different cancers, considering high or low gene
    expression levels, mutation, or copy number status in tumors?
-   Which members of the query set are predicted as partners of
    synthetic lethality interactions?
-   Which members of the query set are associated with cellular
    loss-of-fitness in CRISPR/Cas9 whole-genome drop out screens of
    cancer cell lines (i.e. reduction of cell viability elicited by a
    gene inactivation)? Which targets are prioritized therapeutic
    targets, considering fitness effects and genomic biomarkers in
    combination?


See also the [output views](articles/output.html) that addresses each of the questions above.

## News

-   July 13th 2022: [**1.2.1
    release**](articles/CHANGELOG.html#version-1-2-1)
-   June 24th 2022: [**1.2.0
    release**](articles/CHANGELOG.html#version-1-2-0)
-   March 31st 2022: [**1.0.9
    release**](articles/CHANGELOG.html#version-1-0-9)
-   February 20th 2022: [**1.0.8
    release**](articles/CHANGELOG.html#version-1-0-8)

## Example report

<a href="https://doi.org/10.5281/zenodo.6828941"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6828941.svg" alt="DOI"/></a>

<br>

<img src="img/oncoenrichr_slideshow.gif" align="center" width="560" height="750" border="2,"/>

## Getting started

-   [Installation instructions](articles/installation.html)
-   [How to run oncoEnrichR](articles/running.html)
-   [Explore the various analysis outputs produced with
    oncoEnrichR](articles/output.html)
-   [Annotation resources available in
    oncoEnrichR](articles/annotation_resources.html)

## Contact

sigven AT ifi.uio.no

## Funding and Collaboration

OncoEnrichR is supported by the [Centre for Cancer Cell
Reprogramming](https://www.med.uio.no/cancell/english/) at the
[University of Oslo](https://www.uio.no)/[Oslo University
Hospital](https://radium.no), and [Elixir Norway (Oslo
node)](https://elixir.no/organization/organisation/elixir-uio).

<br> <br>

<p float="left">

<a href="https://www.med.uio.no/cancell/english/">
<img src="img/can-cell.png" width="150"/> </a>        
<a href="https://elixir.no/organization/organisation/elixir-uio">
<img src="img/elixir_norway.png" width="200"/> </a>

</p>
