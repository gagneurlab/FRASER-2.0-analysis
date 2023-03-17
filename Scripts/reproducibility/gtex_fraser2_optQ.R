#'---
#' title: Optimal q obtained for the GTEx tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/reproducibility/gtex_optQ.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds: '`sm expand(config["DATADIR"] + "/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/{dataset}__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/reproducibility/gtex_optQ.html"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
register(MulticoreParam(snakemake@threads))


#+ read in fds objects
fds_ls <- lapply(snakemake@input$fds, function(f) loadFraserDataSet(file=f))

#+ obtain opt q and nr samples and junctions for each tissue
dt <- rbindlist(lapply(fds_ls, function(fds){
    return(data.table(tissue=name(fds), optQ=bestQ(fds, type="jaccard"), nr_introns=nrow(fds), nr_samples=ncol(fds)))
}
))
DT::datatable(
    dt,
    caption = 'GTEx tissue optimal q (FRASER2)',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#+ plot optQ across tissues
ggplot(dt, aes(nr_samples, optQ)) +
    geom_point() +
    theme_bw()

#+ show q fitting curve for all tissues
for(i in length(fds_ls)){
    fds <- fds_ls[[i]]
    g <- plotEncDimSearch(fds, type="jaccard") + ggtitle(name(fds))
    print(g)
}
