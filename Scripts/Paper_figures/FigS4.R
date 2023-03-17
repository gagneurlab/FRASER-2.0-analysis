#'---
#' title: Paper figure S4 (FRASER2 param opt)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS4_recallAt{x}.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 20000
#'   input:
#'     - optimization_plots: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/parameter_optimization_rv_recallAt{x}_ggplots.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS4_recallAt{x}.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS4_recallAt{x}.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake) 

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
# library(gridExtra)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
recall_at <- snakemake@wildcards$x

#+ read in plots and data for the different panels
param_optimizations <- readRDS(snakemake@input$optimization_plots)

#+ extract joined delta and filtering opt results
joined_opt_x_delta <- param_optimizations[["joined_delta_filter_opt_spliceAIvars"]] +
    xlab(expression(Delta~J~" cutoff")) +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
    cowplot::background_grid(major="y", minor="y")

gg_sup <- ggarrange(joined_opt_x_delta)


#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=0.75*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=0.75*page_width, unit=width_unit, dpi=300) 