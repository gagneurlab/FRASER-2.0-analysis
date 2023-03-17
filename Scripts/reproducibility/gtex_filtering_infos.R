#'---
#' title: Compare FRASER2 filtering settings
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/filtering_infos.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 125000
#'   input:
#'     - fraser2_pvals: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK{k}_{q}_minN{n}/PCA__pc0.1/optQ__newFilt/jaccard/pval_matrix__0.3__5__1__FALSE.tsv.gz", dataset=config["tissues_for_detailed_analysis"], k=[20], q=[5, 25, 75, 95], n=[1,5,10])`'
#'     - fraser2_fds: '`sm expand(config["DATADIR"] + "/GTEx_v8/fds/minK{k}_{q}_minN{n}/PCA__pc0.1/savedObjects/{dataset}__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5", dataset=config["tissues_for_detailed_analysis"], k=[20], q=[5, 25, 75, 95], n=[1,5,10])`'
#'   output:
#'     - junction_infos: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/filtering_info.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
.libPaths("~/R/4.1/FRASER2_BB_loss")
library(FRASER)
library(data.table)

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ read in fds files
info_dt <- rbindlist(bpmapply(snakemake@input$fraser2_fds, snakemake@input$fraser2_pvals, 
                 FUN=function(fdsFile, pvalMat){
    # load fds and extract info
    fds <- loadFraserDataSet(file=fdsFile)
    nsamples <- ncol(fds)
    njunc <- nrow(fds)
    
    # load gene level pvals and extract infos
    pvals <- fread(pvalMat)
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
    t <- basename(dirname(dirname(dirname(dirname(dirname(pvalMat))))))
    
    # return data.table with all infos
    data.table(tissue=t, k=k, q=q, n=n, njuncs=njunc, ngenes=ngenes, nsamples=nsamples)
}, SIMPLIFY=FALSE))

#+ write output file
fwrite(info_dt, file=snakemake@output$junction_infos)
