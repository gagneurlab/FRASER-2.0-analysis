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
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig2_recallAt{x}.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/Fig2_recallAt{x}.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
# .libPaths("~/R/4.1/FRASER2")
# library(FRASER)
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
pr_curves_pseudocount <- readRDS(snakemake@input$pr_curve_pseudocount_absplice)
param_optimizations <- readRDS(snakemake@input$optimization_plots)
pr_curves_deltaJaccard <- readRDS(snakemake@input$pr_curve_deltaJaccard_absplice)

#+ extract needed plots from input for pseudocount opt (panels a and b)
maxRank <- 20000
# colors_all <- RColorBrewer::brewer.pal(8, "Blues")
# colors_sub <- colors_all[c(3, 5, 8)]
pr_curve_pc <- pr_curves_pseudocount[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare splice-disrupting\ncandidate variants") +
    # scale_color_manual(values=colors_sub) + 
    # geom_line(size=1.5) +
    # geom_point(size=4) +
    # scale_x_continuous(labels=scales::scientific, limits=c(0, maxRank)) +
    guides(color=guide_legend(order=1, title="Pseudocount", nrow=2), 
           shape=guide_legend(order=2, title="Nominal\np-value cutoff", nrow=2)) + 
    theme_pubr() + 
    theme(
        # legend.position="top", 
        legend.position=c(0.6, 0.3), 
        # legend.box="vertical",
        legend.direction="horizontal",
        legend.background=element_rect(fill='transparent'),
        axis.title=element_text(face="bold"),
        text=element_text(size=font_size))  + 
    cowplot::background_grid(major="xy", minor="xy")
pc_optimization <- param_optimizations[["pc_opt_fixed_rho"]] +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          # axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
          axis.text.x=element_text(angle=90)) + 
    cowplot::background_grid(major="y", minor="y")

#+ extract filtering opt results
# filter_opt <- param_optimizations[["filter_opt_fixed_delta0.3+spliceAIvars"]] +
filter_opt <- param_optimizations[["filter_opt_fixed_delta0.1+spliceAIvars"]] +
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
    cowplot::background_grid(major="y", minor="y")

#+ extract deltaJaccard opt results
# delta_opt <- param_optimizations[["joined_delta_filter_opt_spliceAIvars"]] +
delta_opt <- param_optimizations[["delta_filter_fixed_n1q95"]] +    
    xlab(expression(Delta~J~"cutoff")) +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          # axis.text.x=element_text(angle=45, vjust = 1, hjust=1)) + 
          axis.text.x=element_text(angle=90)) + 
    cowplot::background_grid(major="y", minor="y")

#+ extract PR curve for diff deltaJaccard values
# colors_all <- RColorBrewer::brewer.pal(8, "Oranges")
# colors_sub <- colors_all[c(3, 5, 8)]
pr_curve_dJ <- pr_curves_deltaJaccard[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare splice-disrupting\ncandidate variants") +
    # scale_color_manual(values=colors_sub) + 
    # geom_line(size=1.5) +
    # geom_point(size=4) +
    guides(color=guide_legend(title=expression(Delta~J), nrow=2), 
           shape="none") +
    theme_pubr() + 
    theme(
        # legend.position="top", 
        legend.position=c(0.7, 0.2), 
        legend.direction="horizontal",  
        # legend.box="vertical",
        legend.background=element_rect(fill='transparent'),
        axis.title=element_text(face="bold"),
        text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")

#+ extract filtering stats (nr junctions / genes)
filter_stats <- param_optimizations[["filter_stats"]] + 
    theme_pubr() + 
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")


# Merging legends
legend_1 <- get_legend(pr_curve_pc + guides(shape="none"))
legend_2 <- get_legend(pr_curve_pc + guides(col="none"))
legend_3 <- get_legend(pr_curve_dJ + guides(shape="none")) 
legends <- ggarrange(legend_1, legend_2, legend_3, nrow=3, ncol=1)
# legends <- ggarrange(legends, ggplot() + theme_nothing(), ncol=2) 

#+ combine panels into figure, width=12, height=15
gg_main <- ggarrange(
    pr_curve_pc, #+ guides(shape="none", col="none"), # + scale_x_log10(),
    pc_optimization,
    pr_curve_dJ, # + guides(shape="none", col="none"), # + scale_x_log10(),
    delta_opt,
    filter_opt, 
    filter_stats,
    labels=LETTERS[1:6],
    nrow=3, ncol=2, widths=c(2.2, 3)
)
gg_main
gg_figure <- gg_main

# gg_figure <- ggarrange(
#     pr_curve_pc + guides(shape="none", col="none"), # + scale_x_log10(),
#     pc_optimization,
#     legends, 
#     ggplot() + theme_nothing(),
#     pr_curve_dJ + guides(shape="none", col="none"), # + scale_x_log10(),
#     delta_opt,
#     filter_opt, 
#     filter_stats,
#     labels=c(LETTERS[1:2], c("", ""), LETTERS[3:6]),
#     nrow=4, ncol=2, 
#     widths=c(2.2, 3),
#     heights=c(1,0.3, 1, 1)
# )

gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300) 
