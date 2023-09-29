#'---
#' title: Paper figure Sx (RV recall analysis for splice vicinity set)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_rv_recall_spliceVicinity.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - variant_recall_all_tissues_spliceVicinity: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissues_rv_recall_plots_{varType}.Rds", varType=["rareSplicing"])`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_rv_recall_spliceVicinity.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_rv_recall_spliceVicinity.pdf"`'
#'    - outSvg: '`sm config["PAPER_FIGDIR"] + "/FigSx_rv_recall_spliceVicinity.svg"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/FigSx_rv_recall_spliceVicinity.tiff"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5
scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", 
                                              scales::scientific_format()(x))) }

#+ precision rank plot jointly on all tissues
var_recall_all_VEP <- readRDS(snakemake@input$variant_recall_all_tissues_spliceVicinity[[1]])
maxRank <- 5e6
maxRankForPlot <- 3e6

recall_rank_dt <- (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare splice site vicinity\n(VEP)"]
dt4cutoffs <- (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare splice site vicinity\n(VEP)"]
var_sets_to_show <- c("rare splice site vicinity\n(VEP)")
recall_rank_dt[Method == "FRASER2", Method := "FRASER 2.0"]
dt4cutoffs[Method == "FRASER2", Method := "FRASER 2.0"]
recall_rank_dt[Method == "LeafcutterMD", Method := "LeafCutterMD"]
dt4cutoffs[Method == "LeafcutterMD", Method := "LeafCutterMD"]
methods_to_show <- c("LeafCutterMD", "SPOT", "FRASER", "FRASER 2.0")
recall_rank_dt[, Method:=factor(Method, levels=methods_to_show)]
dt4cutoffs[, Method:=factor(Method, levels=methods_to_show)]
recall_rank_dt <- recall_rank_dt[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
dt4cutoffs <- dt4cutoffs[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
g_var_rank_rec  <- ggplot(recall_rank_dt, aes(rank, recall, col=Method)) +
    # facet_wrap(~snptype) +
    geom_line() +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[3]]$data$slope, 
                col="firebrick", linetype="dashed") + 
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    scale_color_manual(values=c("LeafCutterMD"="orange", "SPOT"="darkolivegreen", "FRASER"="dodgerblue3", "FRASER 2.0"="purple4")) +
    scale_x_continuous(breaks=seq(0, maxRankForPlot, by=1e6),
                       labels=function(x) ifelse(x == 3e6, "", scientific_10(x)),
                       # labels=scientific_10,
                       limits=c(0, maxRankForPlot)) +
    labs(title="", 
         x="Top N outliers", y="Recall of rare splice-disrupting\ncandidate variants") + 
    guides(linetype = "none",
           color=guide_legend(order=1, nrow=2, title=""),
           shape=guide_legend(order=2,nrow=2, title="Nominal\np-value\ncutoff")) + 
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    cowplot::background_grid(major="xy", minor="xy") 
g_var_rank_rec

#+ combine panels into figure
gg_sup <- ggarrange(g_var_rank_rec )
# gg_sup

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=0.5*page_width, height=0.5*page_width, unit=width_unit, dpi=350)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=0.5*page_width, height=0.5*page_width, unit=width_unit, dpi=350)
ggsave(plot=gg_sup, filename=snakemake@output$outSvg, width=0.5*page_width, height=0.5*page_width, unit=width_unit, dpi=350)
ggsave(plot=gg_sup, filename=snakemake@output$outTiff, width=0.5*page_width, height=0.5*page_width, unit=width_unit, dpi=350)
