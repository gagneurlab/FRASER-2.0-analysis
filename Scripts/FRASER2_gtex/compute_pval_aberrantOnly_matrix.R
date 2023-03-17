#'---
#' title: Get pvalue matrix for only aberrant (FDR signif) results
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fds_file:  '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho1_jaccard.h5"`'
#'   output:
#'     - pval_aberrant_matrix: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/{implementation}/optQ__newFilt/jaccard/pval_aberrant_matrix__delta{delta}.tsv.gz"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ read in fds
fds_file <- snakemake@input$fds_file
fds <- loadFraserDataSet(file=fds_file)

#+ get aberrant status
ab <- aberrant(fds, type="jaccard", 
               aggregate=TRUE, by="none",
               padjCutoff=snakemake@config$fdrCutoff,
               deltaPsiCutoff=as.numeric(snakemake@wildcards$delta),
               minCount=5,
               rhoCutoff=NA)

#+ get nominal gene pvalues and set to NA where not aberrant
pvals_gene <- pVals(fds, type="jaccard", level="gene")
stopifnot(all(dim(ab) == dim(pvals_gene)))
pvals_gene[ab == FALSE] <- NA

#+ write output
pval_dt <- data.table(geneID=rownames(pvals_gene), pvals_gene)
fwrite(pval_dt, snakemake@output$pval_aberrant_matrix)
