#'---
#' title: Get all results of FRASER analysis
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 10000
#'   input:
#'    - resultTableGene: '`sm expand("/s/project/gtex_genetic_diagnosis/v8/processed_results" +
#'                          "/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'    - all_res_tables_done: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results" +
#'                          "/aberrant_splicing/results/gencode34/all_fraser1_res.done"`'
#'   type: script
#'---

file.create(snakemake@output$all_res_tables_done)
