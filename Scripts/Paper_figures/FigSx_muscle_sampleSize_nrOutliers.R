#'---
#' title: Paper figure Sx (subsampling on muscle skeletal, outliers per sample)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_muscle_sample_size_nrOultiers.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - power_analysis_gtex: '`sm config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/combined_results.tsv"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_muscle_sample_size_nrOultiers.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_muscle_sample_size_nrOultiers.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm) 
library(data.table)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ power analysis results (transcriptome_wide)
res_gtex <- fread(snakemake@input$power_analysis_gtex)

res_gtex[, aberrant_0_0 := padjustGene <= 0.1 & abs(deltaPsi) >= 0]
res_gtex[, aberrant_0_1 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.1]
res_gtex[, aberrant_0_2 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.2]
res_gtex[, aberrant_0_3 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.3]
res_gtex_melt <- melt(res_gtex, 
                      id.vars=c("sampleID","hgncSymbol","size","sim"), 
                      measure.vars=c("aberrant_0_0", "aberrant_0_1", "aberrant_0_2", "aberrant_0_3"),
                      variable.name="dJ", 
                      value.name="aberrantStatus")
res_gtex_plot <- res_gtex_melt[, .SD[aberrantStatus == TRUE, .N], by=.(size=factor(size), sim, sampleID, dJ)]
res_gtex_plot <- res_gtex_plot[,.(med=median(V1)),by=.(size, sim, dJ)]
res_gtex_plot[, dJ := gsub("_", ".", gsub("aberrant_", "", dJ))]
gtotal_gtex <- ggplot(res_gtex_plot, aes(size, med, col=dJ)) + 
    geom_beeswarm(size=0.5) +
    geom_line(data=res_gtex_plot[, .(mean=mean(med)), by="size,dJ"], aes(size, mean, color=dJ, group=dJ)) +
    ylim(0, 15) + 
    labs(
        x = 'Sample size', 
        y = 'Median of splicing\noutliers per sample') +
    scale_color_brewer(palette="Dark2") +
    guides(color=guide_legend(title=bquote(Delta~J))) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    cowplot::background_grid(major="y", minor="y")
gtotal_gtex

#+ combine panels into figure, width=15, height=12
gg <- ggarrange(
    gtotal_gtex,
    common.legend=TRUE, legend = "bottom"
)
# gg

#+ save figure as png and pdf
ggsave(plot=gg, filename=snakemake@output$outPng, width=0.66*page_width, height=0.4*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg, filename=snakemake@output$outPdf, width=0.66*page_width, height=0.4*page_width, unit=width_unit)