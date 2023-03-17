#'---
#' title: Plot global QQ plot for FRASER2 vs FRASER1
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/mito/fraser2_improvements/minK{k}_{q}_minN{n}/mito_-_fib/qqPlot_comparison_{implementation}.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 30000
#'   input:
#'     - fraser1_qq_data: '`sm config["DATADIR"] + "/mito/processed_data/qqPlot_data/fraser1_qqPlot_data.tsv.gz"`'
#'     - fraser2_qq_data: '`sm config["DATADIR"] + "/mito/processed_data/qqPlot_data/minK{k}_{q}_minN{n}/{implementation}/fraser2_qqPlot_data.tsv.gz"`'
#'   output:
#'     - out_rds: '`sm config["DATADIR"] + "/mito/processed_data/qqPlot_data/minK{k}_{q}_minN{n}/optQ/{implementation}/mito_-_fib/global_qqPlots.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

.libPaths("~/R/4.1/FRASER2_BB_loss")
library(FRASER)
library(ggplot2)
library(RColorBrewer)
register(MulticoreParam(snakemake@threads))

#+ load qq plot data
f1_qq_data <- fread(snakemake@input$fraser1_qq_data)
f2_qq_data <- fread(snakemake@input$fraser2_qq_data)
dt2p <- rbind(f1_qq_data, f2_qq_data, fill=TRUE)

#+ create combined global QQ plot
main <- snakemake@wildcards$dataset
type <- c("psi3", "psi5", "theta", "jaccard")
g <- ggplot(dt2p[plotPoint == TRUE, ], 
            aes(x = exp, y = obs, col = aberrant, label = sampleID, 
                text = paste("<br>SampleID: ", sampleID, "<br>K: ", k, "<br>N: ", n))) + 
    geom_point() + 
    theme_bw() + ggtitle(main) + 
    xlab(expression(-log[10] ~ "(expected P)")) + 
    ylab(expression(-log[10] ~ "(observed P)"))

g$mapping$colour <- quote(type)
# g <- g + scale_color_brewer(palette = "Dark2", name = "Splice metric", 
#                             labels = FRASER:::ggplotLabelPsi(unique(dt2p$type)))

ggplotLabelPsi <- function (type){
    type <- as.character(type)
    vapply(type, FUN = function(x) switch(x, 
                                          jaccard = c(bquote(atop(Intron ~ Jaccard, Index))), 
                                          psi5 = c(bquote(psi[5])), 
                                          psi3 = c(bquote(psi[3])), 
                                          theta = c(bquote(theta))), 
           FUN.VALUE = c(bquote(psi[3])))
}

g <- g + scale_color_manual(values=c("darkorchid4", rev(RColorBrewer::brewer.pal(4, "Blues")[-1])), 
                            name = "Splice metric", 
                            labels = ggplotLabelPsi(unique(dt2p$type)))

# conf.alpha <- 0.05
# dt2p[, `:=`(rank, seq_len(.N)), by = type]
# dt2p[plotPoint == TRUE, `:=`(upper, -log10(qbeta(conf.alpha/2, 
#                                             rank, max(rank) - rank))), by = type]
# dt2p[plotPoint == TRUE, `:=`(lower, -log10(qbeta(1 - conf.alpha/2, 
#                                             rank, max(rank) - rank))), by = type]
typeOrder <- c("theta", "psi5", "psi3", "jaccard")
type2take <- min(which(typeOrder %in% unique(dt2p$type)))
dt2p[type != typeOrder[type2take], `:=`(c("upper", "lower"), list(NA, NA))]

g <- g + 
    geom_ribbon(data = dt2p[plotPoint == TRUE & !is.na(upper), ], 
                aes(x = exp, ymin = lower, ymax = upper, text = NULL), 
                alpha = 0.2, color = "gray") + 
    geom_abline(intercept = 0, slope = 1, col = "firebrick")

#+ combined_qq_plot, fig.width=6, fig.height=4
# g

#+ save plot
saveRDS(object=g, file=snakemake@output$out_rds)
