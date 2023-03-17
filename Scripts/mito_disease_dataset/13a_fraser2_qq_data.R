#'---
#' title: Create data for FRASER2 global QQ plot
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/13a_globalQQ_F2_minK{k}_{q}_minN{n}_{implementation}/.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 75000
#'   input:
#'     - fraser2_fds: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{k}-quantile0.{q}-quantCoverage{n}" + 
#'                "--" + config["mito_annotation"] + 
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   output:
#'     - qq_plot_table: '`sm config["mito_processed_data"] + "/qqPlot_data/minK{k}_{q}_minN{n}/{implementation}/fraser2_qqPlot_data.tsv.gz"`'
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
