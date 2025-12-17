# Rank a gene list according to cancer relevance to a tumor type/site

Function that ranks a list of human gene identifiers with respect to
strength of association to a particular tumor type/site. Underlying
association evidence for the cancer rank per site is pulled from the
Open Targets Platform, see more details
[here](https://sigven.github.io/oncoEnrichR/articles/cancer_gene_rank.html).

## Usage

``` r
cancer_association_rank(
  query = NULL,
  tumor_site = "Breast",
  query_id_type = "symbol",
  ignore_id_err = TRUE,
  include_gene_summary = FALSE,
  oeDB = NULL
)
```

## Arguments

- query:

  character vector with gene/query identifiers (minimum 1, maximum
  2,500)

- tumor_site:

  character indicating primary tumor site of interest

- query_id_type:

  character indicating source of query (one of "uniprot_acc",
  "symbol","entrezgene", or "ensembl_gene", "ensembl_mrna",
  "refseq_transcript_id", "ensembl_protein", "refseq_protein")

- ignore_id_err:

  logical indicating if analysis should continue when uknown query
  identifiers are encountered

- include_gene_summary:

  logical indicating if gene summary (NCBI/UniProt) should be included
  in output data tables

- oeDB:

  oncoEnrichR data repository object - as returned from
  [`load_db()`](https://sigven.github.io/oncoEnrichR/reference/load_db.md)

## Value

A list with two data frames: query genes ranked according to cancer
relevance for the specified tumor site, and unranked query genes (no
evidence found)
