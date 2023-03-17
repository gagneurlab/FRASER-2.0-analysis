#'---
#' title: Paper figure S10 (rare disease analysis, additional plots)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_mito_venn.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - prokisch_f2_vs_f1_plots: '`sm config["DATADIR"] + "/mito/processed_results/ggplot_rds/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/FRASER2_vs_FRASER1_old_filter_ggplots.Rds"`'
#'     - prokisch_fdr_comparison: '`sm config["DATADIR"] + "/mito/processed_results/ggplot_rds/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/FDR_subset_MAF0.001_no_utr_ggplots_F1_old_filter.Rds"`'
#'     - udn_nr_outliers_comparison: '`sm config["DATADIR"] + "/udn/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/nrOutliers_comparison_ggplot.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_venn.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_venn.pdf"`'
#'   type: script
#'---

# #'     - qq_plot_mito: '`sm config["DATADIR"] + "/mito/processed_data/qqPlot_data/minK20_25_minN10/optQ/PCA__pc0.1/mito_-_fib/global_qqPlots.Rds"`'

# #'     - venn_diagram: '`sm "/s/project/fraser/fraser2/figures/FRASER_vs_FRASER2/PCA__pc0.1/fib-minExpr20-quantile0.25-quantCoverage10/gencode34/venn_qopt_deltaJaccard0.1_new_filter.png"`'
# #'     - venn_diagram_fdr_subset: '`sm "/s/project/fraser/fraser2/figures/FRASER_vs_FRASER2/PCA__pc0.1/fib-minExpr20-quantile0.25-quantCoverage10/gencode34/venn_qopt_deltaJaccard0.1_FDR_subset_MAF0.001_no_utr.png"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggvenn)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in plots and data for the different panels
f2_vs_f1_plots <- readRDS(snakemake@input$prokisch_f2_vs_f1_plots)
fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
udn_plots <-readRDS(snakemake@input$udn_nr_outliers_comparison)
# qq_plot_mito <- readRDS(snakemake@input$qq_plot_mito)
# 
# qq_plot_mito <- qq_plot_mito + ggtitle("") + 
#     theme_pubr() +
#     guides(color=guide_legend(nrow=2)) +
#     theme(axis.title=element_text(face="bold"),
#     text=element_text(size=font_size))

nr_outliers_comparison <- f2_vs_f1_plots[["g_numOut_gene"]] +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))
udn_nr_outliers_comparison <- udn_plots[["g_numOut_gene"]] +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))
fdr_comparison <- fdr_comparison_plots[["FDR_comparison_gene"]] +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))
g_venn_w_F1 <- fdr_comparison_plots[["FDR_venn_w_F1"]] + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))
g_venn_wo_F1 <- fdr_comparison_plots[["FDR_venn_wo_F1"]] + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))
g_venn_all <- fdr_comparison_plots[["FDR_venn_all"]] + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size))

#+ combine panels into figure, width=15, height=12
gg_12 <- ggarrange(nr_outliers_comparison,
                   g_venn_w_F1,
                   labels=letters[1:2],
                   nrow=1, ncol=2,
                   widths=c(1, 1.5))
gg_34 <- ggarrange(fdr_comparison,
                   g_venn_all,
                   labels=letters[3:4],
                   nrow=1, ncol=2,
                   widths=c(1, 1.5))
gg_5 <- ggarrange(udn_nr_outliers_comparison,
                  labels=letters[5],
                  nrow=1, ncol=1)
gg_figure_sup <- ggarrange(gg_12,
                           gg_34,
                           gg_5,
                           nrow=3, ncol=1,
                           heights=c(1,1,1))
# gg_figure_sup

#+ save figure as png and pdf
ggsave(plot=gg_figure_sup, filename=snakemake@output$outPng, width=page_width, height=1.0*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure_sup, filename=snakemake@output$outPdf, width=page_width, height=1.0*page_width, unit=width_unit, dpi=300) 
