#'---
#' title: Paper figure 2 (FRASER2 param opt)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_paramOpt_recallAt{x}.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - pr_curve_pseudocount: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_rv_recall_plots_rho1_rareSpliceAI.Rds"`'
#'     - pr_curve_pseudocount_absplice: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_rv_recall_plots_rho1_rareAbSplice.Rds"`'
#'     - optimization_plots: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/parameter_optimization_rv_recallAt{x}_ggplots.Rds"`'
#'     - pr_curve_deltaJaccard: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/minK20_95_minN1/FRASER2_deltaJaccard_rv_recall_plots_pc0.1_rho1__rareSpliceAI.Rds"`'
#'     - pr_curve_deltaJaccard_absplice: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/minK20_95_minN1/FRASER2_deltaJaccard_rv_recall_plots_pc0.1_rho1__rareAbSplice.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_paramOpt_recallAt{x}.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_paramOpt_recallAt{x}.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
recall_at <- snakemake@wildcards$x

#+ read in plots and data for the different panels
inputfile_pr_curve_pc <- snakemake@input$pr_curve_pseudocount_absplice
inputfile_pr_curve_delta <- snakemake@input$pr_curve_deltaJaccard_absplice
pr_curves_pseudocount <- readRDS(inputfile_pr_curve_pc)
param_optimizations <- readRDS(snakemake@input$optimization_plots)
pr_curves_deltaJaccard <- readRDS(inputfile_pr_curve_delta)

#+ extract needed plots from input for pseudocount opt (panels a and b)
maxRank <- 20000
snptype_pc <- strsplit(gsub(".Rds", "", basename(inputfile_pr_curve_pc)), "_", fixed=TRUE)[[1]][7]
snptype_pc <- gsub("rare", "", snptype_pc)
pr_curve_pc <- pr_curves_pseudocount[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare splice-disrupting\ncandidate variants") +
    guides(color=guide_legend(order=1, title="Pseudocount", nrow=2), 
           shape=guide_legend(order=2, title="Nominal\np-value cutoff", nrow=2)) + 
    ggtitle(paste0("Variant annotation tool: ", snptype_pc)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(
        legend.position=c(0.6, 0.3), 
        legend.direction="horizontal",
        legend.background=element_rect(fill='transparent'))  + 
    cowplot::background_grid(major="xy", minor="xy")
# pr_curve_pc
pc_optimization <- param_optimizations[["pc_opt_fixed_rho"]] +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(axis.text.x=element_text(angle=90)) +
    # theme(axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
    cowplot::background_grid(major="y", minor="y")
# pc_optimization

#+ extract deltaJaccard opt results
# delta_opt <- param_optimizations[["joined_delta_filter_opt_spliceAIvars"]] +
delta_opt <- param_optimizations[["delta_filter_fixed_n1q95"]] +    
    xlab(expression(Delta~J~"cutoff")) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(axis.text.x=element_text(angle=90)) +
    # theme(axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
    cowplot::background_grid(major="y", minor="y")
# delta_opt

#+ extract PR curve for diff deltaJaccard values
# colors_all <- RColorBrewer::brewer.pal(8, "Oranges")
# colors_sub <- colors_all[c(3, 5, 8)]
snptype_delta <- strsplit(gsub(".Rds", "", basename(inputfile_pr_curve_delta)), "_", fixed=TRUE)[[1]][9]
snptype_delta <- gsub("rare", "", snptype_delta)
pr_curve_dJ <- pr_curves_deltaJaccard[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare splice-disrupting\ncandidate variants") +
    guides(color=guide_legend(title=expression(Delta~J), nrow=2), 
           shape="none") +
    ggtitle(paste0("Variant annotation tool: ", snptype_delta)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(
        legend.position=c(0.7, 0.2), 
        legend.direction="horizontal",  
        legend.background=element_rect(fill='transparent')) + 
    cowplot::background_grid(major="xy", minor="xy")
# pr_curve_dJ

#+ extract filtering opt results
filter_opt_plot_name <- "filter_opt_fixed_delta0.1+AbSplicevars"
snptype_filter <- strsplit(filter_opt_plot_name, "+", fixed=TRUE)[[1]][2]
snptype_filter <- gsub("vars", "", snptype_filter)
filter_opt <- param_optimizations[[filter_opt_plot_name]] +
    ggtitle(paste0("Variant annotation tool: ", snptype_filter)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    # theme(axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
    cowplot::background_grid(major="y", minor="y")

#+ extract filtering stats (nr junctions / genes)
filter_stats <- param_optimizations[["filter_stats"]] + 
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    cowplot::background_grid(major="y", minor="y")


# Merging legends
legend_1 <- get_legend(pr_curve_pc + guides(shape="none"))
legend_2 <- get_legend(pr_curve_pc + guides(col="none"))
legend_3 <- get_legend(pr_curve_dJ + guides(shape="none")) 
legends <- ggarrange(legend_1, legend_2, legend_3, nrow=3, ncol=1)

#+ combine panels into figure, width=12, height=15
gg_figure <- ggarrange(
    pr_curve_pc, #+ guides(shape="none", col="none"), # + scale_x_log10(),
    pc_optimization,
    pr_curve_dJ, # + guides(shape="none", col="none"), # + scale_x_log10(),
    delta_opt,
    filter_opt, 
    filter_stats,
    labels=LETTERS[1:6],
    font.label=list(size=12, color = "black", face = "bold", family = font),
    nrow=3, ncol=2, widths=c(2.2, 3)
)
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=1.5*page_width, unit=width_unit) 
