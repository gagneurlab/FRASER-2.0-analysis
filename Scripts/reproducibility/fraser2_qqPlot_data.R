#'---
#' title: Create data for FRASER2 global QQ plot
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/fraser2_improvements/minK{k}_{q}_minN{n}/{dataset}/fraser2_qqPlot_data_{implementation}.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 100000
#'   input:
#'     - fraser2_fds: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{q}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'
#'   output:
#'     - qq_plot_table: '`sm config["DATADIR"] + "/{dataset_group}/fraser2_improvements/minK{k}_{q}_minN{n}/optQ/{implementation}/{dataset}/fraser2_qqPlot_data.tsv.gz"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

.libPaths("~/R/4.1/FRASER2")
library(FRASER)
register(MulticoreParam(snakemake@threads))

#+ read in FRASER2 and FRASER1 fds
fds <- loadFraserDataSet(file=snakemake@input$fraser2_fds)

#+ get qq-plot for fraser2
(f2_qq <- plotQQ(fds, global=TRUE, aggregate=FALSE, type="jaccard") )
# f2_qq_data <- as.data.table(f2_qq$data)
f2_qq_data <- as.data.table(f2_qq$layers[[2]]$data)

#+ save qq plot data as table
fwrite(f2_qq_data, file=snakemake@output$qq_plot_table)
