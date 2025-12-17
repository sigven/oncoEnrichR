# Output

## Output files

The output of the analyses performed with oncoEnrichR are provided in
two different formats:

- An interactive [quarto](https://quarto.org)-generated HTML report
- A multi-sheet Excel workbook

### HTML report

Here, we showcase screenshots from the oncoEnrichR interactive HTML
report, following the various types of questions that can be answered
with the tool.

#### Query verification

- Are all my gene/query identifiers correctly identified? Are some of my
  query identifiers outdated, and no longer considered the primary
  identifier?

  

[![](module_screenshots/query_verification.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/query_verification.jpg)

  

#### Poorly characterized genes

- Which members of the query set have a poorly characterized/unknown
  function?

  

[![](module_screenshots/unknown_function.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/unknown_function.jpg)

  

#### Cancer associations

- Which tumor types are associated with genes in the queryset? To what
  extents? Which genes are classified as tumor suppressors,
  proto-oncogenes, or potential cancer driver genes?

  

[![](module_screenshots/cancer_association_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/cancer_association_I.jpg)  
  
  
[![](module_screenshots/cancer_association_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/cancer_association_II.jpg)

  

#### Cancer hallmark evidence

- Which genes have been attributed to the hallmarks of cancer?

  

[![](module_screenshots/cancer_hallmarks.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/cancer_hallmarks.jpg)

  

#### Drug associations

- Are there available targeted cancer drugs for query set members? For
  which cancer indications are there approved drugs targeted towards
  query set members? What are the target tractabilities of all query set
  members?

  

[![](module_screenshots/drug_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/drug_I.jpg)  
  
  
[![](module_screenshots/drug_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/drug_II.jpg)

  

#### Synthetic lethality

- Which genes are predicted to take part in synthetic lethality
  interactions?

  

[![](module_screenshots/synthetic_lethality.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/synthetic_lethality.jpg)

  

#### Gene fitness scores

- Which members in the query set have a significant effect on cell
  fitness, considering large-scale genome-wide CRISPR–Cas9 dropout
  screening data of cancer cell lines?

  

[![](module_screenshots/fitness_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/fitness_I.jpg)  
  
  
[![](module_screenshots/fitness_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/fitness_II.jpg)

  

#### Protein complexes

- Are members of the query set involved in cancer-relevant protein
  complexes?

  

[![](module_screenshots/protein_complex.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/protein_complex.jpg)

  

#### Protein domains

- Which protein domains are most frequent among members of the query
  set?

  

[![](module_screenshots/protein_domains.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/protein_domains.jpg)

  

#### Protein-protein interaction network

- How are protein-protein interactions represented in the query set?
  Which proteins are hubs in the interaction network? Which members of
  the query set form community structures in the network?

  

[![](module_screenshots/ppi_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/ppi_I.jpg)  
  
  
[![](module_screenshots/ppi_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/ppi_II.jpg)

  

#### Regulatory interactions

- Are there transcriptional regulatory interactions in the query set?
  Both TF/regulator and target present? What are the regulatory targets
  for other TFs present? What are the regulators of other targets?

  

[![](module_screenshots/regulatory_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/regulatory_I.jpg)  
  
  
[![](module_screenshots/regulatory_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/regulatory_II.jpg)

  

#### Subcellular compartments

- How are subcellular compartments represented in the query set?

  

[![](module_screenshots/subcell_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/subcell_I.jpg)  
  
  
[![](module_screenshots/subcell_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/subcell_II.jpg)

  

#### Function and pathway enrichment

- Which biological pathways/processes, molecular functions, or cancer
  gene signatures are enriched in the query set?

  

[![](module_screenshots/enrichment_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/enrichment_I.jpg)  
  
  
[![](module_screenshots/enrichment_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/enrichment_II.jpg)

  

#### Tumor aberration frequencies

- To what extent are genes in the query set mutated in tumor samples,
  through copy number alterations or point mutations/indels?

  

[![](module_screenshots/aberration_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/aberration_I.jpg)  
  
  
[![](module_screenshots/aberration_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/aberration_II.jpg)  
  
  
[![](module_screenshots/aberration_III.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/aberration_III.jpg)  
  
  
[![](module_screenshots/aberration_IV.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/aberration_IV.jpg)

  

#### Tumor co-expression

- Are members of the query set co-expressed with other cancer genes in
  tumor samples?

  

[![](module_screenshots/tumor_coexpression.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/tumor_coexpression.jpg)

  

#### Prognostic associations

- Which genes are linked to survival in cancer patients, considering
  expression or methylation levels, or mutation status of query genes in
  tumor samples?

  

[![](module_screenshots/prognostic_I.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/prognostic_I.jpg)  
  
  
[![](module_screenshots/prognostic_II.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/prognostic_II.jpg)

### Excel workbook

- The Excel workbook contains all results from the various analyses
  conducted, organized into multiple sheets (as indicated through red
  circles in the screenshot below)

  

[![](module_screenshots/excel.jpg)](https://sigven.github.io/oncoEnrichR/articles/module_screenshots/excel.jpg)

  
