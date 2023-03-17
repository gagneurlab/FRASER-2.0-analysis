#'---
#' title: Write FRASER2 results table
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds_fitted: '`sm expand(config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho{rho}_jaccard.h5", rho=config["rhoCutoff"], allow_missing=True)`'
#'   output:
#'     - res_table_junction: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_junction.tsv"`'
#'     - res_table_gene: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ read in fds object
fds <- loadFraserDataSet(file=snakemake@input$fds_fitted)

#+ set delta cutoff
deltaCutoff <- as.numeric(snakemake@wildcards$delta)
fdrCutoff <- as.numeric(snakemake@config$fdrCutoff)
rhoCutoff <- as.numeric(snakemake@config$rhoCutoff)
minCountCutoff <- as.numeric(snakemake@config$minCountCutoff)

#+ compute junction level results table
res_junction <- results(fds, psiType="jaccard", aggregate=FALSE, 
                        padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                        rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_junction <- as.data.table(res_junction)

#+ compute gene level results table
res_gene <- results(fds, psiType="jaccard", aggregate=TRUE, 
                    padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                    rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_gene <- as.data.table(res_gene)

#+ save results tables
fwrite(res_junction, file=snakemake@output$res_table_junction, sep="\t")
fwrite(res_gene, file=snakemake@output$res_table_gene, sep="\t")