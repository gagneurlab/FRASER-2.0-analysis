#'---
#' title: Paper figure S3 (rho)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_rhoOpt_recallAt{x}.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - fraser2_fds_pc01: '`sm config["DATADIR"] + "/GTEx_v8/fds/minK20_5_minN10/PCA__pc0.1/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'
#'     - pr_curve: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER2_goodnessOfFit_rv_recall_plots_rareAbSplice.Rds"`'
#'     - rho_optimization: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/parameter_optimization_rv_recallAt{x}_ggplots.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_rhoOpt_recallAt{x}.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_rhoOpt_recallAt{x}.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake) 

#+ echo=FALSE
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
library(ggpubr)
library(cowplot)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in plots and data for the different panels
inputfile_pr_curve_rho <- snakemake@input$pr_curve
pr_curves <- readRDS(inputfile_pr_curve_rho)
param_optimizations <- readRDS(snakemake@input$rho_optimization)

#+ extract needed plots from input
maxRank <- 50000
snptype_rho <- strsplit(gsub(".Rds", "", basename(inputfile_pr_curve_rho)), "_", fixed=TRUE)[[1]][6]
snptype_rho <- gsub("rare", "", snptype_rho)
colors_all <- RColorBrewer::brewer.pal(8, "Greens")
colors_sub <- colors_all[c(2, 5, 8)]
pr_curve <- pr_curves[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare splice-disrupting\ncandidate variants") +
    scale_color_manual(values=colors_sub) +
    # scale_x_log10() +
    # scale_x_continuous(labels=scales::scientific, limits=c(0, maxRank)) +
    guides(color=guide_legend(title=expression(rho))) +
    ggtitle(paste0("Variant annotation tool: ", snptype_rho)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="top", 
          legend.box="vertical")  + 
    cowplot::background_grid(major="xy", minor="xy")
rho_optimization <- param_optimizations[["rho_opt_fixed_pc"]] +
    xlab(expression(rho~" cutoff")) +
    # guides(fill=guide_legend(title=expression(rho))) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(
          # axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
          axis.text.x=element_text(angle=90)) + 
    cowplot::background_grid(major="y", minor="y")

#+ read in fds objects
fds <- loadFraserDataSet(file=snakemake@input$fraser2_fds_pc01)
jidx <- 850 # 836 (quantile=25%) # 835 (quantile=95%)
rho_example <- plotExpression(fds, idx=jidx, type="jaccard", padjCutoff=0.1)
rho_example <- rho_example +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          legend.position="none")
rho_example_predObs <- plotExpectedVsObservedPsi(fds, idx=jidx, type="jaccard", padjCutoff=0.1)
rho_example_predObs <- rho_example_predObs +
    ggtitle("") +
    labs(
        x = "Predicted Intron Jaccard Index",
        y = "Observed Intron Jaccard Index"
    ) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="none")

#+ get ecdf for rho
rho_vals <- rho(fds, type="jaccard")
rho_ecdf <- ggplot(data.table(rho=rho_vals), aes(rho)) + 
    stat_ecdf(geom = "step") + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(title="", #"Empirical Cumulative  Density",
         y = expression("F("~rho~")"), 
         x = expression(rho)) +
    annotation_logticks(sides="b") +
    theme_manuscript(fig_font_size=font_size, fig_font=font)
# rho_ecdf

#+ get legend
legend_1 <- get_legend(pr_curve + guides(col="none"))
legend_2 <- get_legend(pr_curve + guides(shape="none"))
legends <- ggarrange(legend_1, legend_2, nrow=1, ncol=2)

#+ combine panels into figure, width=20, height=15
gg_12 <- ggarrange(pr_curve + guides(shape="none", col="none"),
                       rho_optimization,
                       labels=LETTERS[1:2],
                       nrow=1, ncol=2, widths=c(2, 3))
gg_345 <- ggarrange(rho_example + ggtitle(""), 
                       rho_example_predObs,
                       rho_ecdf,
                       labels=LETTERS[3:5],
                       nrow=1, ncol=3)
gg_figure <- ggarrange(
    legends,
    gg_12, 
    gg_345,
    nrow=3, ncol=1,
    heights=c(1,4,4))
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=1*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=1*page_width, unit=width_unit, dpi=300) 
