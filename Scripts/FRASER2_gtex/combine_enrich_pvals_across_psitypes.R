#'---
#' title: Get pvalue matrices for all cutoffs for enrichment
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/FRASER_pvals__{dPsi}__{minCoverage}__{rho}__{act}.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - pval_matrices: '`sm expand(config["DATADIR"]  + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/{psiType}/pval_matrix__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz", psiType=config["psiTypes"], allow_missing=True)`'
#'   output:
#'     - combined_pval_matrix: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/FRASER_pval_matrix__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'
#'     - combined_pval_matrix_psiOnly: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/FRASER_pval_matrix_psiOnly__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'
#'   type: script
#'---

# #'     - pval_matrices_invertedFilters: '`sm expand(config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/{psiType}/pval_matrix_invertedFilters__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz", psiType=config["psiTypes"], allow_missing=True)`'
# #'     - combined_pval_matrix_invertedFilters: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/FRASER_pval_matrix_invertedFilters__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'
# #'     - combined_pval_matrix_invertedFilters_psiOnly: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK{k}_{quantile}_minN{n}/PCA/FRASER1_optQ/FRASER_pval_matrix_psiOnly_invertedFilters__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'

saveRDS(snakemake, snakemake@log$snakemake)
   
#+ load packages
library(data.table)
library(BiocParallel)

#+ input
input_files <- snakemake@input$pval_matrices
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

#+ combine matrices over psi3 and psi5 only
combined_pvals_psiOnly <- pmin(p_psi3, p_psi5, na.rm=TRUE)
combined_mat_psiOnly <- data.table(geneID=pval_matrices[["psi3"]][,geneID], combined_pvals_psiOnly)

#+ write output
output_matrix_psiOnly <- snakemake@output$combined_pval_matrix_psiOnly
fwrite(combined_mat_psiOnly, file=output_matrix_psiOnly)

# ### repeat for inverted filters
# 
# #+ input_invertedFilters
# input_files <- snakemake@input$pval_matrices_invertedFilters
# 
# #+ read in pvalue matrices per psi type for invertedFilters
# pval_matrices <- lapply(input_files, fread)
# names(pval_matrices) <- sapply(input_files, function(x) basename(dirname(x)) )
# p_psi3  <- as.matrix(pval_matrices[["psi3"]][,-1])
# p_psi5  <- as.matrix(pval_matrices[["psi5"]][,-1])
# p_theta <- as.matrix(pval_matrices[["theta"]][,-1])
# 
# #+ combine matrices over all psi types for invertedFilters
# combined_pvals <- pmin(p_psi3, p_psi5, p_theta, na.rm=TRUE)
# combined_mat <- data.table(geneID=pval_matrices[["psi3"]][,geneID], combined_pvals)
# 
# #+ write output for invertedFilters
# output_matrix <- snakemake@output$combined_pval_matrix_invertedFilters
# fwrite(combined_mat, file=output_matrix)
# 
# #+ combine matrices over psi3 and psi5 only
# combined_pvals_psiOnly <- pmin(p_psi3, p_psi5, na.rm=TRUE)
# combined_mat_psiOnly <- data.table(geneID=pval_matrices[["psi3"]][,geneID], combined_pvals_psiOnly)
# 
# #+ write output for invertedFilters
# output_matrix_psiOnly <- snakemake@output$combined_pval_matrix_invertedFilters_psiOnly
# fwrite(combined_mat_psiOnly, file=output_matrix_psiOnly)
