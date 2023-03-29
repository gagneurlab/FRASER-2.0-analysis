#'---
#' title: Compare FRASER2 filtering settings
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/filtering_infos.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 200000
#'   input:
#'     - fraser2_fds: '`sm expand(config["DATADIR"] + "/GTEx_v8/fds/minK{k}_{q}_minN{n}/PCA__pc0.1/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho1_jaccard.h5", dataset=config["tissues_for_reproducibility"], k=[20], q=[5, 25, 75, 95], n=[1, 5, 10])`'
#'   output:
#'     - junction_infos: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/filtering_info.tsv"`'
#'   type: script
#'---

# q=[5, 25, 75, 95], n=[1, 5, 10]
saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ read in fds files
info_dt <- rbindlist(bplapply(snakemake@input$fraser2_fds, 
                 FUN=function(fdsFile){
    # load fds and extract info
    fds <- loadFraserDataSet(file=fdsFile)
    nsamples <- ncol(fds)
    njunc <- nrow(fds)
    
    # load gene level pvals and extract infos
    pvals <- pVals(fds, type="jaccard", level="gene")
    ngenes <- nrow(pvals)
    
    # extract filter settings from path
    filter_setting <- basename(dirname(dirname(dirname(dirname(fdsFile)))))
    k <- as.numeric(sapply(regmatches(filter_setting, 
                                      regexec(pattern="minK([0-9]+)", text=filter_setting)), 
                           "[", 2))
    q <- as.numeric(sapply(regmatches(filter_setting, 
                                      regexec(pattern="_([0-9]+)_", text=filter_setting)), 
                           "[", 2))
    n <- as.numeric(sapply(regmatches(filter_setting, 
                                      regexec(pattern="minN([0-9]+)", text=filter_setting)), 
                           "[", 2))
    t <- basename(dirname(fdsFile))
    t <- gsub("__optQ__newFilt", "", t)
    
    # return data.table with all infos
    data.table(tissue=t, k=k, q=q, n=n, njuncs=njunc, ngenes=ngenes, nsamples=nsamples)
}))

#+ write output file
fwrite(info_dt, file=snakemake@output$junction_infos)
