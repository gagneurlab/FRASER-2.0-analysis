#'---
#' title: Paper figure S6 (GTEx expression outliers analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_expOutliers.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'    - ae_proportions_table: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/optQ/delta0.1/AE_proportions.tsv"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_expOutliers.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_expOutliers.pdf"`'
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

#+ read in data for plotting
plot_dt <- fread(snakemake@input$ae_proportions_table)
plot_dt[, method:=factor(method, levels=c("LeafCutterMD", "SPOT", "FRASER", "FRASER 2.0", "OUTRIDER"))]
plot_dt <- plot_dt[ae_type == "Underexpression outliers"]
# plot_dt[, ae_type:=factor(ae_type, levels=c("Underexpression outliers", "Overexpression outliers"))]

#+ plot proportions of expression outliers among splice outliers
g_prop <- ggplot(plot_dt[method != "OUTRIDER"], aes(method, prop, fill=method)) +
    # facet_wrap(~ae_type) +
    geom_boxplot() +
    scale_fill_manual(values=c("LeafCutterMD"="orange", "SPOT"="darkolivegreen", "FRASER"="dodgerblue3", "FRASER 2.0"="purple4")) +
    labs(x="", y="Proportion of underexpression\noutliers among splicing outliers") +
    # stat_compare_means(paired=TRUE, size=3) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="none",
          axis.text.x=element_text(size=font_size-2)) +
    cowplot::background_grid(major="y", minor="y")
# g_prop

g_nr <- ggplot(plot_dt, aes(method, nr+1, fill=method)) +
    # facet_wrap(~ae_type) +
    # geom_bar(stat="identity") +
    geom_boxplot() +
    scale_fill_manual(values=c("LeafCutterMD"="orange", "SPOT"="darkolivegreen", "FRASER"="dodgerblue3", "FRASER 2.0"="purple4", "OUTRIDER"="violetred")) +
    labs(x="", y="Underexpression outliers + 1\namong outlier calls") +
    # stat_compare_means(data=plot_dt[method != "OUTRIDER",], 
    #                    paired=TRUE, 
    #                    comparisons=list(c("FRASER", "FRASER2")),
    #                    size=3) +
    scale_y_log10() +
    annotation_logticks(sides="l") +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="none",
          axis.text.x=element_text(size=font_size-2)) +
    cowplot::background_grid(major="y", minor="y")
# g_nr

#+ combine into one figure
gg_sup <- ggarrange(g_prop,
                    g_nr,
                    nrow=1, ncol=2,
                    labels=LETTERS[1:2], 
                    font.label = list(size = 12, color = "black", face = "bold"),
                    widths=c(1,1.25)
)
# gg_sup <- ggarrange(g_prop)
gg_sup

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=1.0*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=1.0*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
# ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=0.5*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
# ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=0.5*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
