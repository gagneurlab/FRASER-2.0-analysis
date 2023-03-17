#'---
#' title: Create data for FRASER1 global QQ plot
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/fraser2_improvements/minK20_5_minN10/{dataset}/fraser1_qqPlot_data_PCA.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 150000
#'   input:
#'     - fraser1_fds: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_splicing/datasets/savedObjects/{dataset}/pvaluesBetaBinomial_junction_theta.h5"`'
#'   output:
#'     - qq_plot_table: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_5_minN10/optQ/PCA/{dataset}/fraser1_qqPlot_data.tsv.gz"`'
#'   type: script
#'---

# #'     - fraser1_fds: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_splicing/datasets/savedObjects/{dataset}_old_filter/pvaluesBetaBinomial_junction_theta.h5"`'

saveRDS(snakemake, snakemake@log$snakemake)

library(FRASER)
register(MulticoreParam(snakemake@threads))

#+ read in FRASER2 and FRASER1 fds
fds <- loadFraserDataSet(file=snakemake@input$fraser1_fds)
    
#+ get qq-plot for fraser2
(f1_qq <- plotQQ(fds, global=TRUE, aggregate=FALSE, type=c("psi3", "psi5", "theta")) )
f1_qq_data <- as.data.table(f1_qq$data)
conv_int_data <- as.data.table(f1_qq$layers[[2]]$data)
f1_qq_data <- f1_qq_data[type != unique(conv_int_data$type), ]
f1_qq_data <- rbind(f1_qq_data, conv_int_data, fill=TRUE)

#+ save qq plot data as table
fwrite(f1_qq_data, file=snakemake@output$qq_plot_table)
