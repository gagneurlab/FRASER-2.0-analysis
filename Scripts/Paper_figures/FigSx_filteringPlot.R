#'---
#' title: Paper figure Sx (intron filtering plot)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_filteringPlot_{dataset}_minK{k}_{quantile}_minN{n}.Rds"`'
#'   threads: 8
#'   resources:
#'     - mem_mb: 64000
#'   input:
#'     - fds_file: '`sm config["DATADIR"] + "/GTEx_v8/fds/savedObjects/raw-{dataset}__jaccard/jaccard.h5"`'
#'     - ss_update_done: '`sm config["DATADIR"] + "/GTEx_v8/fds/savedObjects/raw-{dataset}__jaccard/ss_update.done"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_filteringPlot_{dataset}_minK{k}_{quantile}_minN{n}.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_filteringPlot_{dataset}_minK{k}_{quantile}_minN{n}.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ input
fds_file <- snakemake@input$fds_file
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ load fds
fds <- loadFraserDataSet(file=fds_file)
message("Loaded fds.")

#+ filter junctions based on jaccard index
minExpressionInOneSample <- as.integer(snakemake@wildcards$k)
quantile <- 1 - (as.integer(snakemake@wildcards$quantile)/100)
quantileMinExpression <- as.integer(snakemake@wildcards$n)
minFilterDelta <- snakemake@config$minFilterDelta

#+ filter introns with low read support and corresponding splice sites
fds <- filterExpressionAndVariability(fds,
                                      minExpressionInOneSample=minExpressionInOneSample,
                                      quantile=quantile,
                                      quantileMinExpression=quantileMinExpression,
                                      minDeltaPsi=minFilterDelta, 
                                      filterOnJaccard=TRUE,
                                      filter=FALSE)
# fds <- FRASER:::filterExpression_jaccard(fds,
#                                          minExpressionInOneSample=minExpressionInOneSample,
#                                          quantile=quantile,
#                                          quantileMinExpression=quantileMinExpression,
#                                          filter=FALSE)
# fds <- FRASER:::filterVariability_jaccard(fds, minDelta=minFilterDelta, filter=FALSE)
message(date(), ": Filtering done!")

#+ visualize filtering effect
p <- plotFilterExpression(fds) +
    theme_pubr() + 
    theme(legend.position="bottom",
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) +
    cowplot::background_grid(major="xy", minor="xy")


#+ save figure as png and pdf
ggsave(plot=p, filename=snakemake@output$outPng, width=page_width, height=0.5*page_width, unit=width_unit, dpi=300)
ggsave(plot=p, filename=snakemake@output$outPdf, width=page_width, height=0.5*page_width, unit=width_unit, dpi=300)
