#'---
#' title: Venn diagram across GTEx (FRASER2 vs FRASER)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/venn_minK{k}_{q}_minN{n}_{implementation}_delta{delta}.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - res_fraser2: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_{q}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'     - res_fraser1: '`sm  expand(config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'   output:
#'     - comb_outliers_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK{k}_{q}_minN{n}/{implementation}/optQ/delta{delta}/combined_outliers_venn.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load packages
library(data.table)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in res tables
outliers_f2_all <- unlist(bplapply(snakemake@input$res_fraser2,
                          FUN=function(file){
                              tissue <- basename(dirname(dirname(dirname(file))))
                              res <- fread(file)
                              res[, sampleGene := paste(sampleID, hgncSymbol, tissue, sep="_#_")]
                              return(res[,sampleGene])
                          }))
outliers_f1_all <- unlist(bplapply(snakemake@input$res_fraser1,
                          FUN=function(file){
                              tissue <- basename(dirname(file))
                              tissue <- gsub("(_)+old_filter", "", tissue)
                              res <- fread(file)
                              res[, sampleGene := paste(sampleID, hgncSymbol, tissue, sep="_#_")]
                              return(res[!duplicated(sampleGene),sampleGene])
                          }))
# length(outliers_f2_all)
# length(outliers_f1_all)
all_outliers <- list(
    FRASER = outliers_f1_all,
    FRASER2 = outliers_f2_all)

#+ write combined results table
saveRDS(all_outliers, file=snakemake@output$comb_outliers_rds)
