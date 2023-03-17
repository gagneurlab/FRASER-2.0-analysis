#'---
#' title: Get all jaccard variant enrichment results
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/gtex_nrOutliers_all_delta_cutoffs.Rds"`'
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - nr_outliers_fraser2:  '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_25_minN10/PCA__pc0.1/optQ/delta{delta}/FRASER2_nrOutliers.tsv", dataset=config["tissues_for_reproducibility"], delta=[0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7])`'
#'   output:
#'     - combined_tsv: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/nr_outliers_FRASER_FRASER2_all_deltas.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in all tsvs and combine
dt_combined <- rbindlist(bplapply(snakemake@input$nr_outliers_fraser2, FUN=function(file_f2){
    dt <- fread(file_f2)
    dt <- dt[, .(sampleID, FRASER2_nrOutlierGenes, FRASER2_nrOutlierJunctions_jaccard, tissue)]
    deltaCutoff <- as.numeric(gsub("delta", "", basename(dirname(file_f2))))
    dt[, delta := deltaCutoff]
    return(dt)
}))

#+ write combined file
fwrite(dt_combined, file=snakemake@output$combined_tsv)

# plot
ggplot(dt, aes(factor(delta), FRASER2_nrOutlierGenes+1)) + 
        geom_boxplot(outlier.size=0.5) + 
        scale_y_log10() + 
        annotation_logticks(side="l") +
        labs(x="delta cutoff", y="Outliers per sample + 1") +
        theme_pubr() + 
        cowplot::background_grid(major="y", minor="y")

