#'---
#' title: Create data for FRASER1 global QQ plot
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/fraser2_improvements/minK20_95_minN1/{dataset}/fraser1_qqPlot_data_PCA.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 150000
#'   input:
#'     - fraser1_fds: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_splicing/datasets/savedObjects/{dataset}_old_filter/pvaluesBetaBinomial_junction_theta.h5"`'
#'   output:
#'     - qq_plot_table: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_95_minN1/optQ/PCA/{dataset}/fraser1_qqPlot_data.tsv.gz"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

library(FRASER)
register(MulticoreParam(snakemake@threads))

#+ read in FRASER2 and FRASER1 fds
fds <- loadFraserDataSet(file=snakemake@input$fraser1_fds)

#+ lambda inflation
inflation <- function(p) {
    chisq <- qchisq(p, df=1, lower.tail=FALSE)
    lambda <- median(chisq) / qchisq(0.5, 1)
    lambda
}
lambda_psi5 <- inflation(c(pVals(fds, type="psi5", level="junction")))
lambda_psi3 <- inflation(c(pVals(fds, type="psi3", level="junction")))
lambda_theta <- inflation(c(pVals(fds, type="theta", level="junction")))
lambda_dt <- data.table(lambda=c(lambda_psi5, lambda_psi3, lambda_theta),
                        type=c("psi5", "psi3", "theta"))

#+ get qq-plot for fraser2
(f1_qq <- plotQQ(fds, global=TRUE, aggregate=FALSE, type=c("psi3", "psi5", "theta")) )
f1_qq_data <- as.data.table(f1_qq$data)
conv_int_data <- as.data.table(f1_qq$layers[[2]]$data)
f1_qq_data <- f1_qq_data[type != unique(conv_int_data$type), ]
f1_qq_data <- rbind(f1_qq_data, conv_int_data, fill=TRUE)

#+ add lambdas
f1_qq_data <- merge(f1_qq_data, lambda_dt, by="type", all.x=TRUE)

#+ save qq plot data as table
fwrite(f1_qq_data, file=snakemake@output$qq_plot_table)
