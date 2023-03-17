#'---
#' title: Run FRASER fit for specified q
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fds_fitted: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/predictedMeans_jaccard.h5"`'
#'     - gene_anno_done: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/geneAnnotation.done"`'
#'   output:
#'     - padj_out: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho{rho}_jaccard.h5"`'
#'   type: script
#'---

# #'     - pval_out: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)

#+ input
fds_file <- snakemake@input$fds_fitted
rhoVal <- as.numeric(snakemake@wildcards$rho)
implementation <- snakemake@wildcards$implementation
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

# set pseudocount if requested by implementation
if(grepl("pc", implementation)){
    split_impl <- strsplit(implementation, "__pc", fixed=T)
    implementation <- split_impl[[1]][1]
    pc <- as.numeric(split_impl[[1]][2])
    pseudocount(pc)
}

#+ load fds
fds <- loadFraserDataSet(file=fds_file)

#+ run FRASER fit with hyper param opt for q
for(pt in c("jaccard")){
    message(date(), ": starting BB pval calculation for type ", pt, " ...")
    fds <- calculatePvalues(fds, type=pt, implementation=implementation)
    fds <- calculatePadjValues(fds, type=pt, rhoCutoff=rhoVal)
    gc()
    fds <- calculateZscore(fds, type=pt)
    # fds <- saveFraserDataSet(fds)
}

#+ write output
fds <- saveFraserDataSet(fds)

