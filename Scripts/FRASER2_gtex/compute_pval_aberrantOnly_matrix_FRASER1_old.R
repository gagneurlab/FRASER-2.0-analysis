#'---
#' title: Get pvalue matrix for only aberrant (FDR signif) results (FRASER1 old filtering, 95-quantile)
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fraser1_fds: "/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}__old_filter--gencode34/padjBetaBinomial_theta.h5"
#'   output:
#'     - pval_aberrant_matrix: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/{psiType}/pval_aberrant_matrix__delta{delta}.tsv.gz"`'
#'   type: script
#'---

# #'   py:
# #'   - |
# #'    def get_input_fds(wildcards):
# #'     return config["datasets"][wildcards.dataset]["fds_file"]
# #'     - fraser1_fds: '`sm get_input_fds`'

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ read in fds
fds_file <- snakemake@input$fraser1_fds
fds <- loadFraserDataSet(file=fds_file)
psiTypes <- c("psi5", "psi3", "theta")
fitMetrics(fds) <- psiTypes

#+ get aberrant status
ptype <- snakemake@wildcards$psiType
ab <- aberrant(fds, type=ptype, 
               aggregate=TRUE, by="none",
               padjCutoff=snakemake@config$fdrCutoff,
               deltaPsiCutoff=as.numeric(snakemake@wildcards$delta),
               minCount=5,
               rhoCutoff=NA)

#+ get nominal gene pvalues and set to NA where not aberrant
pvals_gene <- pVals(fds, type=ptype, level="gene")
stopifnot(all(dim(ab) == dim(pvals_gene)))
pvals_gene[ab == FALSE] <- NA

#+ write output
pval_dt <- data.table(geneID=rownames(pvals_gene), pvals_gene)
fwrite(pval_dt, snakemake@output$pval_aberrant_matrix)
