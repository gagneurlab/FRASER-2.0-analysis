#'---
#' title: Get result tables for all tissues
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 2000
#'   input:
#'     - res_gtex: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_25_minN10/PCA__pc0.1/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv", dataset=config["tissues_for_reproducibility"], k=[20], delta=[0.1, 0.2, 0.3])`'
#'     - res_udn: '`sm expand(config["DATADIR"] + "/udn/FRASER2_results/minK{k}_25_minN10/PCA__pc0.1/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv", k=[20], dataset=["Fibroblasts", "blood_polya", "blood_total_rna"], delta=[0.1, 0.2, 0.3])`'  
#'   output:
#'     - all_done: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/all_optQ_results_done.txt"`'
#'   type: script
#'---

#+ create dummy file
out <- file.create(snakemake@output$all_done)
