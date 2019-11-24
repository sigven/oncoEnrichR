# library(magrittr)
# autophagy_candidates <- openxlsx::read.xlsx("/Users/sigven/research/cancell/moncho/data/Autophagy_screens_candidate_list.20190913.xlsx", sheet = 1) %>%
#   janitor::clean_names() %>%
#   dplyr::mutate_if(is.factor, as.character) %>%
#   dplyr::mutate(x_50_increase_starvation = stringr::str_trim(x_50_increase_starvation)) %>%
#   dplyr::mutate(x_50_decrease_starvation = stringr::str_trim(x_50_decrease_starvation))
#
# mitophagy_candidates <- c(autophagy_candidates[autophagy_candidates$x_50_increase_mitophagy != "",]$x_50_increase_mitophagy,
#                           autophagy_candidates[autophagy_candidates$x_50_decrease_mitophagy != "",]$x_50_decrease_mitophagy)
#
# starvation_candidates <- c(autophagy_candidates[autophagy_candidates$x_50_increase_starvation != "",]$x_50_increase_starvation,
#                           autophagy_candidates[autophagy_candidates$x_50_decrease_starvation != "",]$x_50_decrease_starvation)
#
# hits <- list()
# hits[['mitophagy']] <- oncoEnrichR::verify_gene_entries(mitophagy_candidates)
# hits[['starvation']] <- oncoEnrichR::verify_gene_entries(starvation_candidates)
# hits[['universe']] <- oncoEnrichR::verify_gene_entries(autophagy_candidates$si_rna_list)
# hits[['universe']] <- oncoEnrichR::genedb %>% dplyr::filter(!is.na(entrezgene))
# # #
# # #
# gene_enrich_report <- oncoEnrichR::generate_report(hits$starvation$entrezgene,
#                                                           p_title = "StarvationMoncho",
#                                                           p_owner = "Laura Trachsel Moncho / Simonsen lab",
#                                                           background_fname = "/Users/sigven/research/cancell/moncho/moncho_background_starvation.txt",
#                                                           background_enrichment_entrez = hits$universe$entrezgene,
#                                                           background_enrichment_description = "All lipid-binding proteins (n = 198)",
#                                                          gtex_atlasassay_groups = c("g9"))
# oncoEnrichR::write_report(project_directory = "/Users/sigven/research/software/oncoEnrichR", report_name = "StarvationMoncho")
# #
# gene_enrich_report <- oncoEnrichR::generate_report(hits$mitophagy$entrezgene,
#                                                          p_title = "MitophagyMoncho",
#                                                          p_owner = "Laura Trachsel Moncho / Simonsen lab",
#                                                          background_fname = "/Users/sigven/research/cancell/moncho/moncho_background_mitophagy.txt",
#                                                          background_enrichment_entrez = hits$universe$entrezgene,
#                                                          background_enrichment_description = "All lipid-binding proteins (n = 198)")
# oncoEnrichR::write_report(project_directory = "/Users/sigven/research/software/oncoEnrichR", report_name = "MitophagyMoncho")
#
# fgfr4_hits <- read.table(file="/Users/sigven/research/cancell/wesche/fgfr4.txt",stringsAsFactors = F,sep="\n")
# colnames(fgfr4_hits) <- c('symbol')
# hits[['fgfr4']] <- oncoEnrichR::verify_gene_entries(fgfr4_hits$symbol)
#
# gene_enrich_report <- oncoEnrichR::generate_report(hits$fgfr4$entrezgene, p_title = "Fgfr4Wesche",
#                                                          p_owner = "Jørgen Wesche / Wesche lab",
#                                                           background_fname = "/Users/sigven/research/cancell/wesche/project_background_fgfr4.txt",
#                                                           background_enrichment_entrez = NULL,
#                                                           ppi_add_nodes = 1,
#                                                          gtex_atlasassay_groups = c("g9"),
#                                                           background_enrichment_description = "All protein coding genes")
# oncoEnrichR::write_report(project_directory = "/Users/sigven/research/software/oncoEnrichR", report_name = "Fgfr4Wesche")



#
# qsource <- 'symbol'
# target_genes <- list()
# target_genes[['increase']] <- oncoEnrichR::target_disease_associations(qgenes_increase, min_association_score = 0.75)
# target_genes[['decrease']] <- oncoEnrichR::target_disease_associations(qgenes_decrease, min_association_score = 0.75)
# # #
# all_enrichment_results <- list()
# all_enrichment_results[['decrease']] <- data.frame()
# all_enrichment_results[['increase']] <- data.frame()
# # #
# for(tset in c('decrease','increase')){
#   for(c in names(oncoEnrichR::msigdb[['COLLECTION']])){
#     for(subcat in names(oncoEnrichR::msigdb[['COLLECTION']][[c]])){
#       cat(c,subcat,'\n')
#       if(c == "C5"){
#         #enr <- oncoEnrichR::get_go_enrichment(hits[[tset]]$entrezgene, hits[['universe']]$entrezgene,
#         #                                     minGSSize = 5, ontology = subcat, genedb = oncoEnrichR::genedb)
#         enr <- oncoEnrichR::get_go_enrichment(hits[[tset]]$entrezgene,
#                                              minGSSize = 5,
#                                              ontology = subcat,
#                                              genedb = oncoEnrichR::genedb)
#         if(!is.null(enr)){
#           all_enrichment_results[[tset]] <- dplyr::bind_rows(all_enrichment_results[[tset]], enr)
#         }
#       }else{
#         enr <- oncoEnrichR::get_universal_enrichment(hits[[tset]]$entrezgene, genedb = oncoEnrichR::genedb,
#                                                     minGSSize = 4,
#                                                     TERM2GENE = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2GENE,
#                                                     TERM2NAME = oncoEnrichR::msigdb$COLLECTION[[c]][[subcat]]$TERM2NAME,
#                                                     TERM2SOURCE = oncoEnrichR::msigdb$TERM2SOURCE)
#         if(!is.null(enr)){
#           all_enrichment_results[[tset]] <- dplyr::bind_rows(all_enrichment_results[[tset]], enr)
#         }
#       }
#     }
#   }
# }
# # #
# enr <- oncoEnrichR::get_universal_enrichment(hits[['increase']]$entrezgene, genedb = oncoEnrichR::genedb,
#                                             minGSSize = 4,
#                                             TERM2GENE = oncoEnrichR::wikipathwaydb$TERM2GENE,
#                                             TERM2NAME = oncoEnrichR::wikipathwaydb$TERM2NAME,
#                                             TERM2SOURCE = oncoEnrichR::wikipathwaydb$TERM2SOURCE)
# #
# co_exp <- oncoEnrichR::gtex_co_expression(hits$increase$entrezgene, genedb = oncoEnrichR::genedb, gids = c("g32","g9","g45","g10"))

# axl_hits <- openxlsx::read.xlsx("/Users/sigven/research/cancell/stenmark/hits_AXL_10_and_LBD.xlsx",sheet = 1,startRow =1, colNames = T) %>%
#   janitor::clean_names() %>%
#   dplyr::select(accession) %>%
#   dplyr::inner_join(oncoEnrichR::uniprot_xref,by=c("accession" = "uniprot_acc"))
#
# gene_enrich_report <- oncoEnrichR::generate_report_data(axl_hits$symbol, p_title="AxlStenmark", p_owner = "Harald Stenmark/Stenmark lab (Polen collaboration)",background_fname = "/Users/sigven/research/cancell/stenmark/project_background_axl.txt")
