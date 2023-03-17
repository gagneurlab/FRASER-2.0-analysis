#'---
#' title: Compute pvalue matrix for enrichment
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fds_file:  '`sm config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}__size_{sampleSize}_sim{run}--gencode34/padjBetaBinomial_jaccard.h5"`'
#'   output:
#'     - pval_matrix: '`sm config["DATADIR"] + "/GTEx_v8/power_analysis/{dataset}__size{sampleSize}_sim{run}/minK20_25_minN10/PCA__pc0.1/optQ__newFilt/jaccard/pval_matrix__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
library(BBmisc)

#+ source helper functions
source("src/R/enrichment_pval.R")

#+ input
psiType  <- "jaccard"
fds_file <- snakemake@input$fds_file
deltaPsiCutoff    <- as.numeric(snakemake@wildcards$dPsi)
minCoverageCutoff <- as.numeric(snakemake@wildcards$minCoverage)
rhoCutoff         <- as.numeric(snakemake@wildcards$rho)
filterActivation  <- as.logical(snakemake@wildcards$act)
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ output
out_file <- snakemake@output$pval_matrix

#+ load fds
fds <- loadFraserDataSet(file=fds_file)

#+ annotate genes
# fds <- annotateRanges(fds, GRCh=GRCh)
txdb <- AnnotationDbi::loadDb(snakemake@config$datasets[[snakemake@wildcards$dataset]]$txdb)
orgdb <- fread(snakemake@config$datasets[[snakemake@wildcards$dataset]]$orgdb)
seqlevels_fds <- seqlevelsStyle(fds)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, 
                              feature = 'gene_name', featureName = 'hgnc_symbol', keytype = 'gene_id')

#+ compute pvalue matrix
pval_mat <- getEnrichmentPvalMatrix(fds, type=psiType, 
                                    deltaPsiCutoff=deltaPsiCutoff, 
                                    minCoverageCutoff=minCoverageCutoff, 
                                    rhoCutoff=rhoCutoff, 
                                    filterActivation=filterActivation,
                                    BPPARAM=bpparam())

#+ write output
pval_mat <- data.table(geneID=rownames(pval_mat), pval_mat)
fwrite(pval_mat, out_file)

