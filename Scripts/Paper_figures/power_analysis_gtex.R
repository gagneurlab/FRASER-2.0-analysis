#'---
#' title: Power analysis
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/power_analysis_gtex.Rds"`'
#'   threads: 8
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'   - res_simulations: '`sm expand(config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/results/" + 
#'                "gencode34/fraser/{tissue}__size_{sampleSize}_sim{run}/results.tsv", tissue=["Muscle_-_Skeletal"], 
#'                sampleSize=[50, 100, 200, 300, 400, 500, 600, 700], run=[2,3,4,5])`'
#'   output:
#'    - comb_results: '`sm config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/combined_results.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
# .libPaths("~/R/4.1/FRASER2")
# library(FRASER)
library(data.table)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ get results
res_all <- rbindlist(bplapply(snakemake@input$res_simulations, FUN=function(res_file){
    res <- fread(res_file)
    res[, tmp := paste(sampleID, hgncSymbol, sep="_")]
    res <- res[, .(sampleID, hgncSymbol, pValue, pValueGene, padjustGene, psiValue, deltaPsi, tmp)]
    group <- basename(dirname(res_file))
    fsplit <- strsplit(group, "__", fixed=TRUE)[[1]]
    res[, group := group]
    res[, tissue := fsplit[1]]
    fsplit2 <- strsplit(fsplit[2], "_", fixed=TRUE)[[1]]
    res[, size := as.numeric(fsplit2[2])]
    res[, sim  := gsub("sim", "run", fsplit2[3])]
    return(res)
}))
setorder(res_all, size, sim)

#+ write combined results table of solved cases
fwrite(res_all, file=snakemake@output$comb_results, sep="\t")
