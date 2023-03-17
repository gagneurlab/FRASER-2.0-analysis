#'---
#' title: Get pvalue aberant only matrices combined across psitypes
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/FRASER_pvals_aberrant__delta{dPsi}.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - pval_aberrant_matrices: '`sm expand(config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/{psiType}/pval_aberrant_matrix__delta{dPsi}.tsv.gz", psiType=config["psiTypes"], allow_missing=True)`'
#'   output:
#'     - combined_pval_matrix: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/FRASER_pval_aberrant_matrix__delta{dPsi}.tsv.gz"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load packages
library(data.table)
library(BiocParallel)

#+ input
input_files <- snakemake@input$pval_aberrant_matrices
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ read in pvalue matrices per psi type
pval_matrices <- lapply(input_files, fread)
names(pval_matrices) <- sapply(input_files, function(x) basename(dirname(x)) )
p_psi3  <- as.matrix(pval_matrices[["psi3"]][,-1])
p_psi5  <- as.matrix(pval_matrices[["psi5"]][,-1])
p_theta <- as.matrix(pval_matrices[["theta"]][,-1])

#+ combine matrices over all psi types
combined_pvals <- pmin(p_psi3, p_psi5, p_theta, na.rm=TRUE)
combined_mat <- data.table(geneID=pval_matrices[["psi3"]][,geneID], combined_pvals)

#+ write output
output_matrix <- snakemake@output$combined_pval_matrix
fwrite(combined_mat, file=output_matrix)
