validate_query_genes <- function(qgenes,
                               q_id_type = "symbol",
                               qtype = "target",
                               ignore_id_err = F,
                               genedb = NULL,
                               transcript_xref = NULL) {

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  stopifnot(!is.null(q_id_type))
  stopifnot(is.character(qgenes))
  stopifnot(q_id_type == "symbol" |
              q_id_type == "entrezgene" |
              q_id_type == "ensembl_mrna" |
              q_id_type == "ensembl_protein" |
              q_id_type == "refseq_protein" |
              q_id_type == "refseq_transcript_id" |
              q_id_type == "uniprot_acc" |
              q_id_type == "ensembl_gene")
  stopifnot(is.character(qgenes))
  qgenes <- qgenes[!is.na(qgenes)]
  stopifnot(!is.null(genedb))
  stopifnot(!is.null(transcript_xref))
  validate_db_df(genedb, dbtype = "genedb")
  validate_db_df(transcript_xref, dbtype = "transcript_xref")

  alias2entrez <-
    transcript_xref |>
    dplyr::filter(.data$property == "alias") |>
    dplyr::rename(alias = "value") |>
    dplyr::select(-c("property"))

  target_genes <- data.frame(
    'qid' = unique(qgenes),
    stringsAsFactors = F)

  target_genes$qid <- stringr::str_trim(
    target_genes$qid
  )

  if (q_id_type == "entrezgene") {
    target_genes$qid <- as.integer(target_genes$qid)
  }

  if (q_id_type == "ensembl_gene" |
      q_id_type == "ensembl_protein" |
      q_id_type == "ensembl_mrna" |
      q_id_type == "refseq_transcript_id" |
      q_id_type == "refseq_protein_id") {
    target_genes$qid <- stringr::str_replace(
      target_genes$qid, "\\.[0-9]{1,}","")
  }


  gdb <- genedb |>
    dplyr::select(c("symbol",
                  "entrezgene",
                  "name")) |>
    dplyr::distinct()

  queryset <- list()
  queryset[['found']] <- data.frame()
  queryset[['not_found']] <- data.frame()
  queryset[['all']] <- data.frame()
  queryset[['match_status']] <- "perfect_go"

  qtype_id <- 'symbol'
  if (q_id_type == 'entrezgene' |
     q_id_type == 'symbol'){

    qtype_id <- q_id_type

    target_genes <- target_genes |>
      dplyr::left_join(gdb,
                       by = c("qid" = qtype_id), relationship = "many-to-many") |>
      dplyr::mutate(!!rlang::sym(qtype_id) := qid) |>
      dplyr::distinct()

  } else {

    qtype_id <- q_id_type
    if (q_id_type == "refseq_protein") {
      qtype_id <- "refseq_protein_id"
    }
    if (q_id_type == "ensembl_mrna") {
      qtype_id <- "ensembl_transcript_id"
    }
    if (q_id_type == "ensembl_protein") {
      qtype_id <- "ensembl_protein_id"
    }
    if (q_id_type == "ensembl_gene") {
      qtype_id <- "ensembl_gene_id"
    }

    gene_xref_map <-
      transcript_xref |>
      dplyr::filter(.data$property == qtype_id) |>
      dplyr::rename(!!rlang::sym(qtype_id) := value) |>
      dplyr::select(c("entrezgene"), rlang::sym(qtype_id))

    target_genes <- as.data.frame(
      target_genes |>

        ## map query to entrezgene
        dplyr::left_join(
          gene_xref_map,
          by = c("qid" = qtype_id), relationship = "many-to-many") |>
        dplyr::mutate(!!rlang::sym(qtype_id) := qid) |>
        ## append other gene annotations
        dplyr::left_join(gdb, by = c("entrezgene"), relationship = "many-to-many") |>
        dplyr::distinct() |>
        dplyr::group_by(
          .data$symbol,
          .data$entrezgene,
          .data$name) |>
        dplyr::summarise(
          !!rlang::sym(qtype_id) := paste(
            !!rlang::sym(qtype_id), collapse = ","),
          qid = paste(.data$qid, collapse = ","),
          .groups = "drop")
    )

  }

  queryset[['found']] <- target_genes |>
    dplyr::filter(!is.na(.data$symbol) &
                    !is.na(.data$entrezgene)) |>
    dplyr::mutate(alias = F)

  queryset[['not_found']] <- target_genes |>
    dplyr::filter(is.na(.data$symbol) |
                    is.na(.data$entrezgene))

  if (nrow(queryset[['not_found']]) > 0) {
    if (ignore_id_err == T) {
      queryset[['match_status']] <- "imperfect_go"

      if (q_id_type == 'symbol') {
        lgr::lgr$info( paste0("WARNING: ", qtype, " gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid, collapse = ", ")))
        lgr::lgr$info( paste0("Trying to map ", qtype, " gene identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid, collapse = ", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], c("qid")),
            alias2entrez,
            by = c("qid" = "alias"), relationship = "many-to-many")

        ## Check that alias is not an alias for existing query entries (found)
        ## anti_join against found entries

        if (nrow(query_as_alias) > 0) {

          if (nrow(queryset[['found']]) > 0) {
            query_as_alias <- query_as_alias |>
              dplyr::anti_join(queryset[['found']],
                               by = "entrezgene") |>
              dplyr::distinct()
          }
          if (nrow(query_as_alias) > 0) {

            query_as_alias <- query_as_alias |>
              dplyr::left_join(gdb,
                               by = "entrezgene", relationship = "many-to-many") |>
              dplyr::distinct() |>
              dplyr::mutate(alias = T)

            lgr::lgr$info(
              paste0("Mapped query identifiers as gene aliases ",
                     paste0(query_as_alias$qid, collapse = ", ")," ---> ",
                     paste0(query_as_alias$symbol, collapse = ", ")))

            queryset[['found']] <-
              dplyr::bind_rows(queryset[['found']], query_as_alias)
            queryset[['not_found']] <- queryset[['not_found']] |>
              dplyr::anti_join(query_as_alias, by = "qid")

            if (nrow(queryset[['not_found']]) > 0) {
              lgr::lgr$warn(
                paste0(stringr::str_to_title(qtype)," gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid, collapse = ", "),
                       " (make sure that unambiguous primary identifiers/symbols are used)"))
            } else {
              queryset[['match_status']] <- "perfect_go"
            }
          }
        } else {
          lgr::lgr$warn(
            paste0("Query gene identifiers NOT found: ",
                   paste0(queryset[['not_found']]$qid, collapse = ", "),
                   " (make sure that unambiguous primary identifiers/symbols are used)"))
        }

      } else {
        lgr::lgr$warn(paste0(
          stringr::str_to_title(qtype), " gene identifiers NOT found: ",
          paste0(queryset[['not_found']]$qid, collapse = ", ")))
      }

      ## Indicate that processing should stop when encountering invalid query identifiers
    } else {
      queryset[['match_status']] <- "imperfect_stop"

      if (q_id_type == 'symbol') {
        lgr::lgr$warn( paste0("WARNING: ", qtype, " gene identifiers NOT found as primary symbols: ",paste0(queryset[['not_found']]$qid, collapse = ", ")))
        lgr::lgr$info( paste0("Trying to map ", qtype, " gene identifiers as gene aliases/synonyms: ",paste0(queryset[['not_found']]$qid, collapse = ", ")))

        query_as_alias <-
          dplyr::inner_join(
            dplyr::select(queryset[['not_found']], c("qid")),
            alias2entrez,
            by = c("qid" = "alias"), relationship = "many-to-many")

        if (nrow(query_as_alias) > 0) {
          query_as_alias <- query_as_alias |>
            dplyr::left_join(gdb,
                             by = "entrezgene", relationship = "many-to-many") |>
            dplyr::distinct() |>
          dplyr::mutate(alias = T)


          lgr::lgr$info(
            paste0("Mapped ", qtype, " gene identifiers as gene aliases ",
                   paste0(query_as_alias$qid, collapse = ", ")," ---> ",
                   paste0(query_as_alias$symbol, collapse = ", ")))

          queryset[['found']] <-
            dplyr::bind_rows(queryset[['found']], query_as_alias)
          queryset[['not_found']] <- queryset[['not_found']] |>
            dplyr::anti_join(query_as_alias, by = "qid")

          if (nrow(queryset[['not_found']]) > 0) {
            lgr::lgr$info(
              paste0("ERROR: ", qtype, " gene identifiers NOT found: ",
                     paste0(queryset[['not_found']]$qid, collapse = ", "),
                     " (make sure that unambiguous primary identifiers/symbols are used)"))
          } else {
            queryset[['match_status']] <- "perfect_go"

          }
        } else {
          lgr::lgr$info( paste0("ERROR: ", qtype, " gene identifiers NOT found: ",
                         paste0(queryset[['not_found']]$qid, collapse = ", "),
                         " (make sure that unambiguous primary identifiers/symbols are used)"))

        }
      }
      else{
        lgr::lgr$info( paste0("ERROR: ", qtype, " gene identifiers NOT found: ",
                       paste0(queryset[['not_found']]$qid, collapse = ", "),
                       " (make sure that unambiguous primary identifiers/symbols are used)"))
      }
    }
  }
  if (nrow(queryset[['found']]) == length(qgenes)) {
    lgr::lgr$info( paste0('SUCCESS: Identified all genes (n = ',
                             nrow(queryset[['found']]),') in ',qtype,' set'))
  } else {
    if (nrow(queryset[['found']]) == 0) {
      lgr::lgr$info( paste0(
        "ERROR: NO ", qtype, " gene identifiers found: ",
        paste0(target_genes$qid, collapse = ", "),
        " - wrong query_id_type (",q_id_type,")?","\n"))
      queryset[['match_status']] <- "imperfect_stop"
    } else {
      lgr::lgr$info(
        paste0('Identified n = ',
               nrow(queryset[['found']]),' entries in ', qtype,
               ' set (n = ',
               nrow(queryset[['not_found']]),' invalid entries)'))
    }

  }

  if (nrow(queryset[['found']]) > 0) {
    queryset[['found']] <- queryset[['found']] |>
      dplyr::mutate(status = 'found') |>
      dplyr::mutate(
        status = dplyr::if_else(
          .data$status == "found" & .data$alias == T,
          as.character("found_as_alias"),
          as.character(.data$status))
      ) |>
      dplyr::arrange(dplyr::desc(.data$status), .data$symbol) |>
      dplyr::mutate(
        genename = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",
                          .data$entrezgene,"' target='_blank'>",.data$name,"</a>")
      )
  }
  if (nrow(queryset[['not_found']]) > 0) {
    queryset[['not_found']] <- queryset[['not_found']] |>
      dplyr::mutate(status = 'not_found') |>
      dplyr::mutate(genename = NA) |>
      dplyr::arrange(.data$qid)
  }

  if (nrow(queryset[['found']]) > 0 | nrow(queryset[['not_found']]) > 0) {

    queryset[['all']] <- as.data.frame(
      queryset[['not_found']] |>
      dplyr::bind_rows(queryset[['found']]) |>
      dplyr::rename(query_id = "qid") |>
      dplyr::select(c("query_id",
                    "status",
                    "symbol",
                    "genename")) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          symbol = dplyr::if_else(
            .data$status == "not_found",
            as.character(NA),
            as.character(.data$symbol)))
    )
  }

  queryset[['found']]$qid <- NULL
  queryset[['not_found']]$qid <- NULL


  return(queryset)

}

#' Function that validates the oncoEnrichR database object
#'
#' @param oe_db list object with annotation data for oncoEnrichR
#'
#' @keywords internal
#'
validate_db <- function(oe_db) {

  ## check that db is of list type

  db_entries <-
    c("cancerdrugdb",
      "genedb",
      "hpa",
      "ligandreceptordb",
      "otdb",
      "pathwaydb",
      "pfamdb",
      "depmapdb",
      "survivaldb",
      "release_notes",
      "slparalogdb",
      "subcelldb",
      "tcgadb",
      "tftargetdb")
      #"tissuecelldb")

  for (db in db_entries) {
    if (!(db %in% names(oe_db))) {
      lgr::lgr$info(paste0("ERROR: '",db,"' NOT found in oncoEnrichR db object"))
      return(-1)
    }

  }
  return(0)
}

#' Function that validates a particular db object (data.frame) in oncoEnrichR
#'
#' @param df data.frame with annotation data for oncoEnrichR
#' @param dbtype type of oncoEnrichR datasource
#'
#' @keywords internal
#'

validate_db_df <- function(df, dbtype = "genedb") {

  val <- assertthat::validate_that(
    is.data.frame(df)
  )
  if (!is.logical(val)) {
    stop(val)
  }

  dbtypes <- c("genedb",
               "hpadb",
               "tcga_aberration",
               "tcga_diagnosis_code",
               "tcga_site_code",
               "tcga_clinical_strata_code",
               "tcga_coexpression",
               "tcga_recurrent_variants",
               "protein_complex",
               "protein_domain",
               "dorothea",
               "ligand_receptor_db",
               "ligand_receptor_xref",
               "transcript_xref",
               "biogrid",
               "compartments",
               "oeDB",
               "tf_target_interactions",
               "go_gganatogram",
               "opentarget_disease_assoc",
               "opentarget_disease_site_rank",
               #"enrichment_db_hpa_singlecell",
               #"enrichment_db_hpa_tissue",
               "ppi_nodes",
               "slparalog",
               "fitness_scores",
               "target_priority_scores",
               "survival_km_cshl",
               "ppi_edges")
  if (!(dbtype %in% dbtypes)) {
    stop(
      paste0("dbtype '",dbtype,
             "' not recognized, possible s are: ",
             paste(sort(dbtypes), collapse = ", ")))
  }

  if (dbtype == "genedb") {
    cols <- c('symbol',
              'entrezgene',
              'oncogene',
              'oncogene_confidence_level',
              'tumor_suppressor',
              'tsg_confidence_level',
              'cancer_driver',
              'hgnc_id',
              'tsg_support',
              'oncogene_support',
              'driver_support',
              'cancergene_evidence',
              'go_term_link',
              'name',
              'gene_summary',
              'gene_biotype',
              'AB_tractability_category',
              'AB_tractability_support',
              'SM_tractability_category',
              'SM_tractability_support',
              'genename',
              'targeted_cancer_drugs_lp',
              'targeted_cancer_drugs_ep',
              'approved_drugs',
              'num_go_terms',
              'cancer_max_rank',
              'unknown_function_rank',
              'has_gene_summary')
  }

  if (dbtype == "opentarget_disease_assoc") {
    cols <- c("disease_efo_id",
              "direct_ot_association",
              "ot_association_score",
              "ot_datatype_support",
              "ensembl_gene_id",
              "efo_name",
              "primary_site",
              "cancer_phenotype",
              "ot_link")
  }

  if (dbtype == "biogrid") {
    cols <- c("entrezgene_A",
              "entrezgene_B",
              "method",
              "pmid",
              "throughput")
  }

  if (dbtype == "tcga_aberration") {
    cols <- c("symbol",
              "variant_type",
              "samples_mutated",
              "tot_samples",
              "percent_mutated",
              "percentile",
              "decile",
              "clinical_strata_code",
              "diagnosis_code",
              "site_code")
  }

  if (dbtype == "slparalog") {
    cols <- c("entrezgene_A1",
              "symbol_A1",
              "entrezgene_A2",
              "symbol_A2",
              "prediction_score",
              "prediction_percentile",
              "sequence_identity_pct",
              "family_size",
              "conservation_score",
              "ppi_overlap")
  }

  if (dbtype == "tcga_coexpression") {
    cols <- c("symbol",
              "symbol_partner",
              "r",
              "p_value",
              "tumor",
              "tumor_suppressor",
              "oncogene",
              "cancer_driver")
  }

  if (dbtype == "tcga_diagnosis_code") {
    cols <- c("primary_diagnosis",
              "diagnosis_code")
  }

  if (dbtype == "tcga_site_code") {
    cols <- c("primary_site",
              "site_code")
  }

  if (dbtype == "tcga_clinical_strata_code") {
    cols <- c("clinical_strata",
              "clinical_strata_code")
  }


  if (dbtype == "opentarget_disease_site_rank") {
    cols <- c("primary_site",
              "ensembl_gene_id",
              "tissue_assoc_score",
              "tissue_assoc_rank",
              "global_assoc_score",
              "global_assoc_rank")
  }

  if (dbtype == "enrichment_db_hpa_singlecell") {
    cols <- c("ensembl_gene_id",
              "category",
              "cell_type")
  }

  if (dbtype == "enrichment_db_hpa_tissue") {
    cols <- c("ensembl_gene_id",
              "category",
              "tissue")
  }

  if (dbtype == "hpadb") {
    cols <- c("ensembl_gene_id",
              "property",
              "value")
  }


  if (dbtype == "tf_target_interactions") {
    cols = c("queryset_overlap",
                 "literature_support",
                 "interaction_sources",
                 "regulator_cancer_max_rank",
                 "target_cancer_max_rank",
                 "regulator",
                 "regulator_name",
                 "target",
                 "target_name",
                 "confidence_level",
                 "mode_of_regulation")
  }


  if (dbtype == "go_gganatogram") {
    cols <- c("ggcompartment",
                  "go_id")
  }
  if (dbtype == "protein_complex") {
    cols <- c('complex_id',
              'complex_name',
              'purification_method',
              'complex_comment',
              'disease_comment',
              'sources',
              'confidence',
              'complex_literature',
              'complex_literature_support')
  }
  if (dbtype == "protein_domain") {
    cols <- c('uniprot_acc',
              'pfam_id',
              'pfam_short_name',
              'pfam_long_name',
              'domain_freq')
  }

  if (dbtype == "ligand_receptor_xref") {
    cols <- c('interaction_id',
              'symbol',
              'class')
  }

  if (dbtype == "ligand_receptor_db") {
    cols <- c('interaction_id',
              'interaction_name',
              'annotation',
              'pathway_name',
              'interaction_members',
              'ligand',
              'receptor',
              'agonist',
              'antagonist',
              'co_A_receptor',
              'co_I_receptor',
              'literature_support')
  }
  if (dbtype == "dorothea") {
    cols <- c('regulator',
              'target',
              'interaction_sources',
              'confidence_level',
              'mode_of_regulation',
              'tf_target_literature_support',
              'tf_target_literature')
  }
  if (dbtype == "transcript_xref") {
    cols <- c('entrezgene',
              'property',
              'value')
  }
  if (dbtype == "survival_km_cshl") {
    cols <- c('symbol',
              'tcga_cohort',
              'z_score')
  }

  if (dbtype == "target_priority_scores") {
    cols <- c('symbol',
              'gene_id',
              'priority_score',
              'tumor_type')
  }

  if (dbtype == "fitness_scores") {
    cols <- c('symbol',
              'model_name',
              'model_id',
              'loss_of_fitness',
              'scaled_BF',
              'tissue',
              'cancer_type',
              'tissue_status',
              'entrezgene',
              'sample_site',
              'gene_id_project_score')
  }

  if (dbtype == "comppidb") {
    cols <- c('uniprot_acc',
              'go_id',
              'go_term',
              'confidence',
              'annotation_source',
              'annotation_type')
  }

  if (dbtype == "compartments") {
    cols <- c('entrezgene',
              'go_id',
              'go_term',
              'confidence',
              'annotation_source',
              'annotation_channel')
  }


  if (dbtype == "ppi_nodes") {
    cols <- c('symbol',
              'entrezgene',
              'genename',
              'query_node',
              'cancer_driver',
              'id',
              'tumor_suppressor',
              'oncogene')
  }

  if (dbtype == "ppi_edges") {
    cols <- c('symbol_A',
              'symbol_B',
              'entrezgene_A',
              'entrezgene_B',
              'oncogene_A',
              'oncogene_B',
              'tsgene_A',
              'tsgene_B',
              'cdriver_A',
              'cdriver_B',
              'querynode_A',
              'querynode_B',
              'weight',
              'from',
              'to',
              'interaction_symbol')

  }
  assertable::assert_colnames(df,
                              colnames = cols,
                              only_colnames = F,
                              quiet = T)

  if (dbtype == 'transcript_xref') {
    identifiers_expected <- c('alias',
                     'ensembl_gene_id',
                     'ensembl_protein_id',
                     'ensembl_transcript_id',
                     'refseq_protein_id',
                     'refseq_transcript_id',
                     'symbol',
                     'uniprot_acc')

    identifiers_found <- sort(unique(df$property))

    assertthat::assert_that(
      identical(identifiers_expected, identifiers_found),
      msg = paste0("Identifier types present in 'transcript_xref' data frame",
                   " ('property' column) does not have all necessary s")
    )

    assertthat::assert_that(
      typeof(df$entrezgene) == "integer",
      msg = paste0("Type of 'entrezgene' column in 'transcript_xref' ",
                   "is not of type 'integer'")
    )


  }

  return(0)

}

add_excel_sheet <- function(
  report = NULL,
  workbook = NULL,
  analysis_output = "disease_association",
  tableStyle = "TableStyleMedium15") {

  invisible(assertthat::assert_that(!is.null(report)))
  invisible(assertthat::assert_that(!is.null(report$data)))
  invisible(assertthat::assert_that(!is.null(workbook)))

  target_df <- data.frame()

  lgr::lgr$info(
    paste0("Adding Excel sheet to workbook - ", analysis_output))

  if (analysis_output == "settings") {
    target_df <- data.frame(
      category = 'CANCER_ASSOCIATION',
      configuration = 'show_disease',
      value = as.character(report$config$show$disease),
      stringsAsFactors = F
    )
    target_df <- target_df |>
      dplyr::bind_rows(
        data.frame(
          category = 'CANCER_ASSOCIATION',
          configuration = 'show_top_diseases',
          value = as.character(report$config$disease$show_top_diseases),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'CANCER_HALLMARK',
          configuration = 'show_hallmark',
          value = as.character(report$config$show$cancer_hallmark),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'UNKNOWN_FUNCTION',
          configuration = 'show_unknown_function',
          value = as.character(report$config$show$unknown_function),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'DRUG_KNOWN',
          configuration = 'show_drug',
          value = as.character(report$config$show$drug),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SYNTHETIC_LETHALITY',
          configuration = 'show_synleth',
          value = as.character(report$config$show$synleth),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PROTEIN_COMPLEX',
          configuration = 'show_protein_complex',
          value = as.character(report$config$show$protein_complex),
          stringsAsFactors = F
        ),

        data.frame(
          category = 'ENRICHMENT',
          configuration = 'show_enrichment',
          value = as.character(report$config$show$enrichment),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_p_value_adj',
          value = as.character(report$config$enrichment$p_adjust_method),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_p_value_cutoff',
          value = as.character(report$config$enrichment$p_value_cutoff),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_q_value_cutoff',
          value = as.character(report$config$enrichment$q_value_cutoff),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_min_geneset_size',
          value = as.character(report$config$enrichment$min_gs_size),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_max_geneset_size',
          value = as.character(report$config$enrichment$max_gs_size),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'enrichment_simplify_go',
          value = as.character(report$config$enrichment$simplify_go),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'bgset_description',
          value = as.character(report$config$enrichment$bgset_description),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ENRICHMENT',
          configuration = 'bgset_size',
          value = as.character(report$config$enrichment$bgset_size),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'REGULATORY',
          configuration = 'show_regulatory',
          value = as.character(report$config$show$regulatory),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'REGULATORY',
          configuration = 'regulatory_min_confidence',
          value = as.character(report$config$regulatory$min_confidence),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SUBCELL',
          configuration = 'show_subcell',
          value = as.character(report$config$show$subcellcomp),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SUBCELL',
          configuration = 'subcellcomp_min_confidence',
          value = as.character(report$config$subcellcomp$minimum_confidence),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SUBCELL',
          configuration = 'subcellcomp_min_channels',
          value = as.character(report$config$subcellcomp$minimum_channels),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SUBCELL',
          configuration = 'subcellcomp_channels',
          value = paste(
            as.character(
              report$config$subcellcomp$channels),
            collapse = ", "),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'SUBCELL',
          configuration = 'show_cytosol',
          value = as.character(report$config$subcellcomp$show_cytosol),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'CELL_TISSUE',
          configuration = 'show_cell_tissue',
          value = as.character(report$config$show$cell_tissue),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'show_ppi',
          value = as.character(report$config$show$ppi),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_string_min_score',
          value = as.character(report$config$ppi$string$minimum_score),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_string_network_type',
          value = as.character(report$config$ppi$string$network_type),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_biogrid_min_evidence',
          value = as.character(report$config$ppi$biogrid$minimum_evidence),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_show_drugs',
          value = as.character(report$config$ppi$string$show_drugs),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_add_nodes',
          value = as.character(report$config$ppi$string$add_nodes),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PPI',
          configuration = 'ppi_show_isolated_nodes',
          value = as.character(report$config$ppi$string$show_isolated_nodes),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'FITNESS',
          configuration = 'show_fitness',
          value = as.character(report$config$show$fitness),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'FITNESS',
          configuration = 'fitness_max_score',
          value = as.character(report$config$fitness$max_BF_score),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'ABERRATION',
          configuration = 'show_aberration',
          value = as.character(report$config$show$aberration),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'COEXPRESSION',
          configuration = 'show_coexpression',
          value = as.character(report$config$show$coexpression),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'LIGAND_RECEPTOR',
          configuration = 'show_ligand_receptor',
          value = as.character(report$config$show$ligand_receptor),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'PROGNOSTIC',
          configuration = 'show_prognostic',
          value = as.character(report$config$show$cancer_prognosis),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'project_title',
          value = as.character(report$config$project_title),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'project_description',
          value = as.character(report$config$project_description),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'project_owner',
          value = as.character(report$config$project_owner),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'query_ignore_err',
          value = as.character(report$config$query$ignore_err),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'query_id_type',
          value = as.character(report$config$query$id_type),
          stringsAsFactors = F
        ),
        data.frame(
          category = 'OTHER',
          configuration = 'bgset_id_type',
          value = as.character(report$config$bgset$id_type),
          stringsAsFactors = F
        )
      )

  }

  if (analysis_output == "query") {
    if (is.data.frame(report$data$query$target)) {
      if (NROW(report$data$query$target) > 0) {
        target_df <- report$data$query$target |>
          dplyr::mutate(
            genename =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$genename)
                )
              )
          )
      }
    }
  }

  if (analysis_output == "unknown_function") {
    if (is.data.frame(report$data$unknown_function$hits_df)) {
      if (NROW(report$data$unknown_function$hits_df) > 0) {
        target_df <- report$data$unknown_function$hits_df |>
          dplyr::rename(go_terms = .data$go_term_link) |>
          dplyr::mutate(
            annotation_source = "Gene Ontology (2022-12)/NCBI Gene/UniProt (2022_05)",
            version = NA) |>
          dplyr::select(c("annotation_source",
                        "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          ) |>
          dplyr::mutate(
            go_terms =
              stringr::str_trim(
                textclean::replace_html(.data$go_terms)
              )
          ) |>
          dplyr::mutate(
            gene_summary =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$gene_summary)
                )
              )
          )
      }
    }
  }

  if (analysis_output == "cancer_association") {
    if (is.data.frame(report$data$disease$target)) {
      if (NROW(report$data$disease$target) > 0) {
        target_df <- report$data$disease$target |>
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) |>
          dplyr::select(-c("cancer_association_links",
                           "disease_association_links",
                           "cancergene_evidence",
                           "gene_summary")) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          )
      }
    }
  }

  if (analysis_output == "cancer_hallmark") {
    if (is.data.frame(report$data$cancer_hallmark$target)) {
      if (NROW(report$data$cancer_hallmark$target) > 0) {
        target_df <- report$data$cancer_hallmark$target |>
          dplyr::mutate(
            symbol =
              stringr::str_squish(
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
              )
          ) |>
          dplyr::select(-c("literature_support"))
      }
    }
  }


  if (analysis_output == "drug_known") {
    if (is.data.frame(report$data$drug$target_drugs)) {
      if (NROW(report$data$drug$target_drugs)) {
        target_df <- report$data$drug$target_drugs |>
          dplyr::mutate(
            annotation_source = report$config$resources$opentargets$name,
            version = report$config$resources$opentargets$version) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          ) |>
          dplyr::mutate(
            approved_drugs =
              stringr::str_trim(
                textclean::replace_html(.data$approved_drugs)
              )
          ) |>
          dplyr::mutate(
            drugs_late_phase =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(.data$drugs_late_phase)
                  )
                ),
                " , ",
                ", "
              )
          ) |>
          dplyr::mutate(
            drugs_early_phase =
              stringr::str_replace_all(
                stringr::str_squish(
                  stringr::str_trim(
                    textclean::replace_html(.data$drugs_early_phase)
                  )
                ),
                " , ",
                ", "
              )
          )
      }
    }
  }

  if (analysis_output == "drug_tractability") {
    if (is.data.frame(report$data$drug$tractability_ab) &
        is.data.frame(report$data$drug$tractability_sm)) {
      if (NROW(report$data$drug$tractability_ab) |
          NROW(report$data$drug$tractability_sm)) {

        target_df <- data.frame()

        if (NROW(report$data$drug$tractability_sm)) {
          df <- report$data$drug$tractability_sm |>
            dplyr::mutate(
              annotation_source = report$config$resources$opentargets$name,
              version = report$config$resources$opentargets$version) |>
            dplyr::rename(tractability_category = "SM_tractability_category",
                          tractability_support = "SM_tractability_support") |>
            dplyr::mutate(tractability_drugtype = "Small molecules/compounds") |>
            dplyr::select(c("annotation_source", "version",
                          "tractability_drugtype"),
                          dplyr::everything()) |>
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
            ) |>
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(.data$tractability_support)
                )
            )

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }

        if (NROW(report$data$drug$tractability_ab)) {
          df <- report$data$drug$tractability_ab |>
            dplyr::mutate(
              annotation_source = report$config$resources$opentargets$name,
              version = report$config$resources$opentargets$version) |>
            dplyr::rename(tractability_category = .data$AB_tractability_category,
                          tractability_support = .data$AB_tractability_support) |>
            dplyr::mutate(tractability_drugtype = "Antibody") |>
            dplyr::select(c("annotation_source", "version",
                          "tractability_drugtype"),
                          dplyr::everything()) |>
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$symbol)
                )
            ) |>
            dplyr::mutate(
              tractability_support =
                stringr::str_trim(
                  textclean::replace_html(.data$tractability_support)
                )
            )

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }
      }
    }
  }

  if (analysis_output == "ligand_receptor") {
    for (c in c('secreted_signaling','cell_cell_contact',
               'ecm_receptor')) {

      if (is.data.frame(report$data$ligand_receptor[[c]])) {
        if (NROW(report$data$ligand_receptor[[c]]) > 0) {

          df <- report$data$ligand_receptor[[c]] |>
            dplyr::mutate(
              annotation_source = report$config$resources[['cellchatdb']]$name,
              version = report$config$resources[['cellchatdb']]$version) |>
            dplyr::select(c("annotation_source", "version"),
                          dplyr::everything()) |>
            dplyr::mutate(
              literature_support =
                stringr::str_trim(
                  textclean::replace_html(.data$literature_support)
                )
            )

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }
      }

    }
  }

  if (analysis_output == "protein_domain") {

    if (is.data.frame(report$data$protein_domain$target)) {
      if (NROW(report$data$protein_domain$target) > 0) {

        df <- report$data$protein_domain$target |>
          dplyr::mutate(
            annotation_source = report$config$resources[['pfam']]$name,
            version = report$config$resources[['pfam']]$version) |>
          dplyr::select(c("annotation_source",
                        "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            protein_domain =
              stringr::str_trim(
                textclean::replace_html(.data$protein_domain)
              ),
            target_genes =
              stringr::str_trim(
                textclean::replace_html(.data$target_genes)
              )
          )

        target_df <- target_df |>
          dplyr::bind_rows(df)
      }
    }
  }

  if (analysis_output == "protein_complex") {

    for (c in c('omnipath','humap2')) {

      if (is.data.frame(report$data$protein_complex[[c]])) {
        if (NROW(report$data$protein_complex[[c]]) > 0) {

          res_name <- c
          if(c == 'omnipath'){
            res_name <- 'omnipathr'
          }

          df <- report$data$protein_complex[[c]] |>
            dplyr::mutate(
              annotation_source = report$config$resources[[res_name]]$name,
              version = report$config$resources[[res_name]]$version) |>
            dplyr::select(c("annotation_source",
                            "version"),
                          dplyr::everything()) |>
            dplyr::mutate(
              complex_genes =
                stringr::str_replace_all(
                  stringr::str_squish(
                    stringr::str_trim(
                      textclean::replace_html(.data$complex_genes)
                    )
                  ),
                  " , ",
                  ", "
                )
            )

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }
      }
    }

    if (NROW(target_df) > 0 &
       "literature" %in% colnames(target_df)) {
      target_df <- target_df |>
        dplyr::mutate(
          literature =
            stringr::str_trim(
              textclean::replace_html(.data$literature)
            )
        ) |>
        dplyr::mutate(
          complex_name =
            stringr::str_trim(
              textclean::replace_html(.data$complex_name)
            )
        )
    }

  }


  if (analysis_output == "prognostic_association_I") {
    if (is.data.frame(report$data$cancer_prognosis$hpa$assocs)) {
      if (NROW(report$data$cancer_prognosis$hpa$assocs) > 0) {
        target_df <- report$data$cancer_prognosis$hpa$assocs |>
          dplyr::mutate(
            annotation_source = report$config$resources$hpa$name,
            version = report$config$resources$hpa$version) |>
          dplyr::select(c("annotation_source",
                        "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            tumor_types =
              stringr::str_trim(
                textclean::replace_html(.data$tumor_types)
              )
          )
      }
    }
  }

  if (analysis_output == "synthetic_lethality") {
    for (t in c('single_pair_member','both_in_pair')) {

      if (is.data.frame(report$data$synleth[[t]])) {
        if (NROW(report$data$synleth[[t]]) > 0) {
          df <- report$data$synleth[[t]] |>
            dplyr::mutate(
              annotation_source = "De Kegel et al., Cell Systems, 2021",
              version = "v1") |>
            dplyr::mutate(feature_type = t) |>
            dplyr::mutate(
              genename_A =
                stringr::str_trim(
                  textclean::replace_html(.data$genename_A)
                )
            ) |>
            dplyr::mutate(
              genename_B =
                stringr::str_trim(
                  textclean::replace_html(.data$genename_B)
                )
            ) |>
            dplyr::arrange(.data$feature_type,
                           dplyr::desc(.data$prediction_score)) |>
            dplyr::select(c("annotation_source",
                            "version",
                          "feature_type"), dplyr::everything())

          target_df <- target_df |>
            dplyr::bind_rows(df)

        }
      }
    }
  }


  if (analysis_output == "prognostic_association_II") {
    for (t in c('cna','mut','exp','meth')) {
      if (is.data.frame(report$data$cancer_prognosis$km_cshl$assocs[[t]])) {
        if (NROW(report$data$cancer_prognosis$km_cshl$assocs[[t]]) > 0) {
          df <- report$data$cancer_prognosis$km_cshl$assocs[[t]] |>
            dplyr::mutate(
              annotation_source = report$config$resources$tcga_survival$name,
              version = report$config$resources$tcga_survival$version) |>
            dplyr::mutate(feature_type = t) |>
            dplyr::arrange(.data$feature_type, .data$z_score) |>
            dplyr::select(c("annotation_source", "version"),
                          dplyr::everything())

          target_df <- target_df |>
            dplyr::bind_rows(df)

        }
      }
    }
  }

  if (analysis_output == "coexpression") {

    ## co-expression
    if (is.data.frame(report$data$tcga$coexpression)) {
      if (NROW(report$data$tcga$coexpression) > 0) {
        target_df <- report$data$tcga$coexpression |>
          dplyr::mutate(
            annotation_source = report$config$resources$tcga$name,
            version = report$config$resources$tcga$version) |>
          dplyr::rename(tcga_cohort = "tumor") |>
          dplyr::select(-c("entrezgene")) |>
          dplyr::rename(partner_oncogene = "oncogene",
                        partner_tumor_suppressor = "tumor_suppressor",
                        partner_cancer_driver = "cancer_driver") |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything())
      }
    }
  }

  if (analysis_output == "regulatory") {

    ## regulatory interactions
    for (c in c('pancancer','global')) {

      if (is.data.frame(report$data$regulatory$interactions[[c]])) {
        if (NROW(report$data$regulatory$interactions[[c]]) > 0) {
          df <- report$data$regulatory$interactions[[c]] |>
            dplyr::mutate(
              dorothea_collection = c,
              annotation_source = report$config$resources$dorothea$name,
              version = report$config$resources$dorothea$version) |>
            dplyr::select(c("annotation_source",
                            "version",
                          "dorothea_collection"),
                          dplyr::everything()) |>
            dplyr::mutate(
              target_name =
                stringr::str_trim(
                  textclean::replace_html(.data$target_name)
                )
            ) |>
            dplyr::mutate(
              regulator_name =
                stringr::str_trim(
                  textclean::replace_html(.data$regulator_name)
                )
            ) |>
            dplyr::mutate(
              literature_support =
                stringr::str_trim(
                  textclean::replace_html(.data$literature_support)
                )
            )

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }
      }
    }

  }

  if (analysis_output == "recurrent_variants") {

    if (is.data.frame(report$data$tcga$recurrent_variants)) {
      if (NROW(report$data$tcga$recurrent_variants) > 0) {
        df <-
          report$data$tcga$recurrent_variants

        colnames(df) <- tolower(colnames(df))
        df <- as.data.frame(
          df |>
            dplyr::mutate(
              site_recurrence = as.numeric(.data$site_recurrence)
            ) |>
            dplyr::arrange(
              dplyr::desc(.data$total_recurrence),
              dplyr::desc(.data$site_recurrence)) |>
            dplyr::mutate(
              annotation_source = report$config$resources$tcga$name,
              version = report$config$resources$tcga$version) |>
            dplyr::mutate(
              ensembl_gene_id =
                stringr::str_trim(
                  textclean::replace_html(.data$ensembl_gene_id)
                )
            ) |>
            dplyr::mutate(
              ensembl_transcript_id =
                stringr::str_trim(
                  textclean::replace_html(.data$ensembl_transcript_id)
                )
            ) |>
            dplyr::mutate(
              protein_domain =
                stringr::str_trim(
                  textclean::replace_html(.data$protein_domain)
                )
            ) |>
            dplyr::mutate(
              cosmic_mutation_id =
                stringr::str_trim(
                  textclean::replace_html(.data$cosmic_mutation_id)
                )
            ) |>
            dplyr::mutate(
              site_recurrence = paste(.data$primary_site,
                                      .data$site_recurrence, sep =":")
            ) |>
            dplyr::group_by(
              .data$symbol, .data$consequence,
              .data$protein_change, .data$protein_domain,
              .data$mutation_hotspot,
              .data$loss_of_function,
              .data$ensembl_gene_id,
              .data$ensembl_transcript_id,
              .data$total_recurrence,
              .data$cosmic_mutation_id
            ) |>
            dplyr::summarise(site_recurrence = paste(
              .data$site_recurrence, collapse = ", "
            ), .groups = "drop") |>
            dplyr::arrange(
              dplyr::desc(.data$total_recurrence)
            )

        )

        target_df <- target_df |>
          dplyr::bind_rows(df)
      }
    }
  }

  if (analysis_output == "aberration") {

    ## cna aberrations
    for (t in c('cna_ampl','cna_homdel')) {
      if (is.data.frame(report$data$tcga$aberration$table[[t]])) {
        if (NROW(report$data$tcga$aberration$table[[t]]) > 0) {
          df <-
            report$data$tcga$aberration$table[[t]] |>
            dplyr::mutate(
              annotation_source = report$config$resources$tcga$name,
              version = report$config$resources$tcga$version) |>
            dplyr::mutate(
              symbol =
                stringr::str_trim(
                  textclean::replace_html(.data$gene)
                )
            ) |>
            dplyr::select(-c("gene")) |>
            dplyr::select(c("annotation_source", "version",
                          "symbol"), dplyr::everything())

          target_df <- target_df |>
            dplyr::bind_rows(df)
        }
      }
    }

  }

  if (analysis_output == "subcellcomp") {
    if (is.data.frame(report$data$subcellcomp$all)) {
      if (NROW(report$data$subcellcomp$all) > 0) {
        target_df <- report$data$subcellcomp$all |>
          dplyr::mutate(
            annotation_source = report$config$resources$compartments$name,
            version = report$config$resources$compartments$version) |>
          dplyr::select(-c("ggcompartment")) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything()) |>
          dplyr::mutate(
            compartment =
              stringr::str_trim(
                textclean::replace_html(.data$compartment)
              )
          ) |>
          dplyr::mutate(
            genename =
              stringr::str_trim(
                textclean::replace_html(.data$genename)
              )
          )
      }
    }
  }

  if (analysis_output == "enrichment") {
    enrichment_df <- data.frame()
    for (e in c('go','wikipathway','kegg','msigdb','netpath')) {
      if (NROW(report$data$enrichment[[e]]) > 0) {
        enrichment_df <- enrichment_df |>
          dplyr::bind_rows(
            report$data$enrichment[[e]] |>
              dplyr::mutate(
                annotation_source = report$config$resources[[e]]$name,
                version = report$config$resources[[e]]$version) |>
              dplyr::rename(category = "db",
                            entrezgene = "gene_id")
          ) |>
          dplyr::select(-c("gene_symbol_link",
                           "description_link")) |>
          dplyr::select(c("annotation_source", "version",
                        "category", "description"),
                        dplyr::everything())
      }
    }
    target_df <- enrichment_df
  }



  if (analysis_output == "fitness_scores") {
    if (is.data.frame(report$data$fitness$fitness_scores$targets)) {
      if (NROW(report$data$fitness$fitness_scores$targets) > 0) {
        target_df <- report$data$fitness$fitness_scores$targets |>
          dplyr::mutate(
            annotation_source = report$config$resources$depmap$name,
            version = report$config$resources$depmap$version) |>
          dplyr::select(-c("symbol_link_ps", "model_link_ps",
                           "n_gene")) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything())
      }
    }
  }

  if (analysis_output == "fitness_prioritized") {
    if (is.data.frame(report$data$fitness$target_priority_scores$targets)) {
      if (NROW(report$data$fitness$target_priority_scores$targets) > 0) {
        target_df <- report$data$fitness$target_priority_scores$targets |>
          dplyr::mutate(
            annotation_source = report$config$resources$depmap$name,
            version = report$config$resources$depmap$version) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything())
      }
    }
  }

  if (analysis_output == "ppi_string") {

    target_df <- data.frame()
    if (is.data.frame(report$data$ppi$string$complete_network$edges)) {
      if (NROW(report$data$ppi$string$complete_network$edges) > 0) {
        target_df <- report$data$ppi$string$complete_network$edges |>
          dplyr::filter(!is.na(.data$entrezgene_A) &
                          !is.na(.data$entrezgene_B)) |>
          dplyr::arrange(dplyr::desc(.data$score)) |>
          dplyr::rename(is_target_A = "querynode_A",
                        is_target_B = "querynode_B",
                        string_score = "score",
                        string_nscore = "nscore",
                        string_fscore = "fscore",
                        string_pscore = "pscore",
                        string_ascore = "ascore",
                        string_escore = "escore",
                        string_dscore = "dscore",
                        string_tscore = "tscore"
                        ) |>
          dplyr::select(c("symbol_A",
                        "symbol_B",
                        "is_target_A",
                        "is_target_B",
                        "string_score",
                        "string_nscore",
                        "string_fscore",
                        "string_pscore",
                        "string_ascore",
                        "string_escore",
                        "string_dscore",
                        "string_tscore")) |>
          dplyr::mutate(
            annotation_source = report$config$resources$string$name,
            version = report$config$resources$string$version) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything())
      }
    }

  }

  if (analysis_output == "ppi_biogrid") {

    target_df <- data.frame()
    if (is.data.frame(report$data$ppi$biogrid$complete_network$edges)) {
      if (NROW(report$data$ppi$biogrid$complete_network$edges) > 0) {
        target_df <- report$data$ppi$biogrid$complete_network$edges |>
          dplyr::filter(!is.na(.data$entrezgene_A) &
                          !is.na(.data$entrezgene_B)) |>
          dplyr::rename(is_target_A = "querynode_A",
                        is_target_B = "querynode_B",
                        biogrid_evidence = "evidence",
                        biogrid_num_evidence = "num_evidence_items"
          ) |>
          dplyr::select(c("symbol_A",
                          "symbol_B",
                          "is_target_A",
                          "is_target_B",
                          "biogrid_evidence",
                          "biogrid_num_evidence")) |>
          dplyr::arrange(
            dplyr::desc(.data$biogrid_num_evidence)) |>
          dplyr::mutate(
            annotation_source = report$config$resources$biogrid$name,
            version = report$config$resources$biogrid$version) |>
          dplyr::select(c("annotation_source", "version"),
                        dplyr::everything())
      }
    }

  }

  if (analysis_output == "cell_tissue") {

    target_df <- data.frame()
    for (e in c("tissue_enrichment", "scRNA_enrichment")) {
      if (is.data.frame(report$data$cell_tissue[[e]]$per_gene)) {
        if (NROW(report$data$cell_tissue[[e]]$per_gene) > 0) {
          if (e == "tissue_enrichment") {
            df <-
              report$data$cell_tissue[[e]]$per_gene |>
              dplyr::mutate(
                annotation_source = report$config$resources$gtex$name,
                version = report$config$resources$gtex$version,
                category = stringr::str_replace(e,"_enrichment","")) |>
              dplyr::select(c("annotation_source",
                            "version",
                            "category"),
                            dplyr::everything()) |>
              dplyr::rename(tissue_or_celltype = "tissue")

          } else {
            df <-
              report$data$cell_tissue[[e]]$per_gene |>
              dplyr::mutate(
                annotation_source = report$config$resources$hpa$name,
                version = report$config$resources$hpa$version,
                category = stringr::str_replace(e,"_enrichment","")) |>
              dplyr::select(c("annotation_source",
                              "version",
                            "category"), dplyr::everything()) |>
              dplyr::rename(tissue_or_celltype = "cell_type")
          }
          target_df <- target_df |>
            dplyr::bind_rows(df) |>
            dplyr::mutate(
              genename =
                stringr::str_trim(
                  textclean::replace_html(.data$genename)
                )
            )
        }
      }
    }
  }

  if (nrow(target_df) > 0) {

    analysis_output = stringr::str_replace(
      analysis_output,"_ps_","_")

    openxlsx::addWorksheet(workbook,
                           sheetName = toupper(analysis_output))

    ## set automatic column widths
    openxlsx::setColWidths(workbook,
                           sheet = toupper(analysis_output),
                           cols = 1:ncol(target_df),
                           widths = "auto")

    ## write with default Excel Table style
    openxlsx::writeDataTable(workbook,
                             sheet = toupper(analysis_output),
                             x = target_df,
                             startRow = 1,
                             startCol = 1,
                             colNames = TRUE,
                             tableStyle = tableStyle)
  }

  return(workbook)
}


file_is_writable <- function(path) {
  assertthat::is.string(path) &&
    file.exists(path) &&
    assertthat::is.writeable(path)
}

#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

utils::globalVariables(c("."))
