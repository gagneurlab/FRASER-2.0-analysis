#'---
#' title: Paper figure Sx (precision recall plot of FDR signficant results)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_precision_recall_FDRsignif.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 32000
#'   input:
#'     - variant_recall_all_tissues: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_fdrSignif_allTissues_rv_recall_plots_{snptype}.Rds", snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_precision_recall_FDRsignif.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_precision_recall_FDRsignif.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
library(eulerr)
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
var_recall_all_VEP <- readRDS(snakemake@input$variant_recall_all_tissues[[1]])
var_recall_all_SpliceAI <- readRDS(snakemake@input$variant_recall_all_tissues[[2]])
var_recall_all_MMSplice <- readRDS(snakemake@input$variant_recall_all_tissues[[3]])
var_recall_all_AbSplice <- readRDS(snakemake@input$variant_recall_all_tissues[[4]])
# maxRank <- 5e6
# maxRankForPlot <- 3e6

recall_rank_dt <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                         (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare SpliceAI"],
                         (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare MMSplice"],
                         (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare AbSplice"])
dt4cutoffs <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                     (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare SpliceAI"],
                     (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare MMSplice"],
                     (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare AbSplice"])
recall_rank_dt[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
dt4cutoffs[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
# var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare AbSplice")
var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice")
recall_rank_dt[Method == "FRASER2", Method := "FRASER 2.0"]
dt4cutoffs[Method == "FRASER2", Method := "FRASER 2.0"]
methods_to_show <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER 2.0")
recall_rank_dt <- recall_rank_dt[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
dt4cutoffs <- dt4cutoffs[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
g_var_prec_rec  <- ggplot(recall_rank_dt[Method != "totalPossibleRank"], 
                          aes(x=recall, y=precision, col=Method)) +
    facet_wrap(~snptype) +
    geom_line() +
    geom_point(data=dt4cutoffs, aes(x=recall, y=precision, color=Method, shape=Cutoff), size=3) +
    labs(title=paste(dataset), 
         x="Recall of rare splice-disrupting candidate variants", y="Precision") +
    grids(color="white") +
    scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                       labels=metric_labels) +
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
           color=guide_legend(title="Method", order = 1)) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    cowplot::background_grid(major="xy", minor="xy") 
g_var_prec_rec

#+ save figure as png and pdf
ggsave(plot=g_var_prec_rec, filename=snakemake@output$outPng, width=0.75*page_width, height=0.5*page_width, unit=width_unit)
ggsave(plot=g_var_prec_rec, filename=snakemake@output$outPdf, width=0.75*page_width, height=0.5*page_width, unit=width_unit)
