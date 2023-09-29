#'---
#' title: Paper figure 3 (GTEx analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/fig2.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - global_qq: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/optQ/PCA__pc0.1/Skin_-_Not_Sun_Exposed_Suprapubic/global_qqPlots.Rds"`'
#'     - variant_recall_all_tissues: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissues_rv_recall_plots_{varType}.Rds", varType=["rareSplicing", "rareSpliceSite", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'     - nr_outliers_per_sample:  '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/nrOutliers_comparison_ggplot.Rds"`'
#'     - comb_outliers_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/optQ/delta0.1/combined_outliers_venn.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig2.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/Fig2.pdf"`'
#'    - outSvg: '`sm config["PAPER_FIGDIR"] + "/Fig2.svg"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/Fig2.tiff"`'
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

#+ read in plots and data for the different panels
global_qq_plot <- readRDS(snakemake@input$global_qq)
global_qq_plot <- global_qq_plot + ggtitle("")  + 
    labs(
        x = expression(bold(-log[10]("expected"~"P"))),
        y = expression(bold(-log[10]("observed"~"P")))
    ) + 
    # guides(col = guide_legend(nrow = 4, title="")) +
    guides(col = guide_legend(nrow = 2, title="")) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    cowplot::background_grid(major="xy", minor="xy")
global_qq_plot$layers[[1]]$aes_params$size <- point_size
# global_qq_plot
global_qq_plot$data[, .(lambda=unique(lambda)), by="type"]

#+ read in outliers for venn plotting
all_outliers <- readRDS(snakemake@input$comb_outliers_rds)
names(all_outliers) <- c("FRASER", "FRASER 2.0")
#+ create proportional venn diagram
venn_plot <- plot(euler(all_outliers), 
                  quantities = list(TRUE, type=c("counts", "percent"), 
                                    fontsize=font_size, 
                                    fontfamily=font),
                  labels=list(fontsize=font_size, 
                              fontfamily=font),
                  fills=c("white", "white"),
                  edges=c("dodgerblue3", "purple4"),
                
)
venn_plot

#+ precision rank plot jointly on all tissues
var_recall_all_VEP <- readRDS(snakemake@input$variant_recall_all_tissues[[1]])
var_recall_all_VEP_spliceSite <- readRDS(snakemake@input$variant_recall_all_tissues[[2]])
var_recall_all_SpliceAI <- readRDS(snakemake@input$variant_recall_all_tissues[[3]])
var_recall_all_MMSplice <- readRDS(snakemake@input$variant_recall_all_tissues[[4]])
var_recall_all_AbSplice <- readRDS(snakemake@input$variant_recall_all_tissues[[5]])
maxRank <- 5e6
maxRankForPlot <- 3e6

recall_rank_dt <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                        (var_recall_all_VEP_spliceSite[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare direct splice site\nvariant (VEP)"],
                        (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare SpliceAI"],
                        (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare MMSplice"],
                        (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare AbSplice"])
dt4cutoffs <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                     (var_recall_all_VEP_spliceSite[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare direct splice site\nvariant (VEP)"],
                     (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare SpliceAI"],
                     (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare MMSplice"],
                     (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare AbSplice"])
recall_rank_dt[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare direct splice site\nvariant (VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
dt4cutoffs[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare direct splice site\nvariant (VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
# var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare AbSplice")
var_sets_to_show <- c("rare direct splice site\nvariant (VEP)", #"rare splice site vicinity\n(VEP)",
                        "rare MMSplice", "rare SpliceAI", "rare AbSplice")
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
    facet_wrap(~snptype) +
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

### nr of outliers per sample comparison
nr_outliers_per_sample_dt_gtex <- readRDS(snakemake@input$nr_outliers_per_sample)$data

# nr of outliers per sample (GTEx)
nr_outliers_per_sample_dt_gtex[, dataset_group := "GTEx"]
nr_outliers_per_sample_dt_gtex[, dataset_label := paste0(tissue, " (N=", uniqueN(sampleID), ")"), by=tissue]
nr_outliers_per_sample_dt_gtex[, tissue := NULL]
nr_outliers_per_sample_dt_gtex[grepl("Brain", dataset_label), dataset_label:=paste0("Brain (N=", uniqueN(sampleID), ")")]
tissue_order <- nr_outliers_per_sample_dt_gtex[,median(nr_outliers_per_sample),by="dataset_label,method"][method == "FRASER"][order(-V1), dataset_label]
nr_outliers_per_sample_dt_gtex[, dataset_label := factor(dataset_label, levels=tissue_order)]
nr_outliers_per_sample_dt_gtex[method == "FRASER2", method := "FRASER 2.0"]

nr_outliers_per_sample_plot <- ggplot(nr_outliers_per_sample_dt_gtex, 
                                      aes(dataset_label, nr_outliers_per_sample+1, col=method)) +
    # facet_grid(~dataset_group, scales="free_x", space="free_x") +
    geom_boxplot(outlier.size=point_size) +
    labs(x="", y="Splicing outliers per sample + 1") +
    scale_y_log10() + 
    scale_color_manual(values=c("dodgerblue3", "purple4")) +
    # coord_flip() + 
    annotation_logticks(sides="l") +
    # annotation_logticks(sides="b") + 
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme(
        legend.title = element_blank(),
        # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        # plot.margin=unit(c(0.5,0.5,0.5,2.75), "cm"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=1),
        plot.margin=unit(c(0.5,0.5,1,0.5), "cm")
    ) +
    cowplot::background_grid(major="y", minor="y")
# nr_outliers_per_sample_plot

#+ combine panels into figure, width=12, height=20
col1 <- ggarrange(global_qq_plot, 
                  venn_plot,
                  labels=LETTERS[c(1,3)], 
                  font.label = list(size = 12, color = "black", face = "bold"),#, family = "Arial"),
                  align="hv",
                  nrow=2, ncol=1,
                  heights=c(1.8,1))
row1 <- ggarrange(col1, 
                ggarrange(
                    g_var_rank_rec,
                    labels=LETTERS[2]),
                    font.label = list(size = 12, color = "black", face = "bold"),#, family = "Arial"),
                align="hv",
                nrow=1, ncol=2,
                widths=c(1,2))
row2 <- ggarrange(nr_outliers_per_sample_plot,
                  labels=LETTERS[4], 
                  font.label = list(size = 12, color = "black", face = "bold"),#, family = "Arial"),
                  nrow=1, ncol=1)
gg_figure <- ggarrange(row1, 
                       row2,
                       align="hv",
                       nrow=2, ncol=1,
                       heights=c(1,1))
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_figure, filename=snakemake@output$outSvg, width=page_width, height=1.5*page_width, unit=width_unit, dpi=350)
ggsave(plot=gg_figure, filename=snakemake@output$outTiff, width=page_width, height=1.5*page_width, unit=width_unit, dpi=350)
