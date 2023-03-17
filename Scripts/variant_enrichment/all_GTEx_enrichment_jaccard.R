#'---
#' title: Get all jaccard variant enrichment results
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 2000
#'   input:
#'     - gg_rds_opt_q:  '`sm config["htmlOutputPath"] + "/variant_enrichment/GTEx_tissues/optQ_enrichment_jaccard.png"`'
#'     - gg_rds_opt_q_newFilt:  '`sm config["htmlOutputPath"] + "/variant_enrichment/GTEx_tissues/optQ_newFilt_enrichment_jaccard.png"`'
#'   output:
#'     - all_done: '`sm config["DATADIR"] + "/variant_enrichment/GTEx_tissues/jaccard_var_enrich_all_done.txt"`'
#'   type: script
#'---


# #'     - gg_rds_all_q:  '`sm expand(config["htmlOutputPath"] + "/variant_enrichment/GTEx_tissues/q{q}_enrichment_jaccard.png", q=config["q"])`'

#+ create dummy file
file.create(snakemake@output$all_done)


