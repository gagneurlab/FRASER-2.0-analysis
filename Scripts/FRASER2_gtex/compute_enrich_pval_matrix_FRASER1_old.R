#'---
#' title: Compute pvalue matrix for enrichment (FRASER1 old filtering, 95-quantile)
#' author: Ines Scheller
#' wb:
#'   py:
#'   - |
#'    def get_input_fds(wildcards):
#'     return config["datasets"][wildcards.dataset]["fds_file"]
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fraser1_fds: '`sm get_input_fds`'
#'   output:
#'     - pval_matrix: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/{psiType}/pval_matrix__{dPsi}__{minCoverage}__{rho}__{act}.tsv.gz"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2_BB_loss")
library(FRASER)
library(data.table)
library(BBmisc)

#+ source helper functions
source("src/R/enrichment_pval.R")

#+ input
psiType  <- snakemake@wildcards$psiType
fds_file <- snakemake@input$fraser1_fds
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

