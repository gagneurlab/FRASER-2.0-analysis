#'---
#' title: Paper figure S7 (GTEx Venn)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_gtex_venn_recall.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - comb_outliers_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/optQ/delta0.1/combined_outliers_venn.Rds"`'
#'     - pr_curves_venn: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_vs_FRASER2_venn_rv_recall_plots_rareSpliceAI.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_gtex_venn_recall.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_gtex_venn_recall.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(ggplot2)
library(ggpubr)
library(cowplot)
library(eulerr)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in outliers for venn plotting
all_outliers <- readRDS(snakemake@input$comb_outliers_rds)

#+ enrichment of splice variants for venn diagram categories
pr_curves <- readRDS(snakemake@input$pr_curves_venn)
maxRank <- 50000
pr_curve_venn <- pr_curves[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare\nsplice affecting variants") +
    # guides(col=guide_legend(order=2, nrow=1)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="bottom", 
          legend.box="vertical",
          legend.title=element_blank())  + 
    cowplot::background_grid(major="xy", minor="xy") +
    guides(col=guide_legend(nrow=3, order=1))
pr_curve_venn

#+ compile figure
gg_figure <- ggarrange(
            pr_curve_venn
          )
# gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=0.5*page_width, height=0.6*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=0.5*page_width, height=0.6*page_width, unit=width_unit) 
