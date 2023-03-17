#'---
#' title: Annotate genes to introns
#' author: Ines Scheller
#' wb:
#'   threads: 5
#'   resources:
#'     - mem_mb: 30000
#'   input:
#'     - fds_fitted: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/predictedMeans_jaccard.h5"`'
#'   output:
#'     - annotation_done: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/geneAnnotation.done"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)

#+ input
fds_file <- snakemake@input$fds_fitted
implementation <- snakemake@wildcards$implementation 
txdb <- AnnotationDbi::loadDb(snakemake@config$datasets[[snakemake@wildcards$dataset]]$txdb)
orgdb <- fread(snakemake@config$datasets[[snakemake@wildcards$dataset]]$orgdb)
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ load fds
fds <- loadFraserDataSet(file=fds_file)

#+ annotate genes
seqlevels_fds <- seqlevelsStyle(fds)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, 
                              # filter=list('gene_type'='protein_coding'),
                              feature = 'gene_name', 
                              featureName = 'hgnc_symbol', 
                              keytype = 'gene_id')

#+ write output
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$annotation_done)