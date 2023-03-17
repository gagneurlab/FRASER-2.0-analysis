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

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in data for plotting
plot_dt <- fread(snakemake@input$ae_proportions_table)
plot_dt[, method:=factor(method, levels=c("FRASER", "FRASER2", "OUTRIDER"))]
plot_dt[, ae_type:=factor(ae_type, levels=c("Underexpression outliers", "Overexpression outliers"))]

#+ plot proportions of expression outliers among splice outliers
g_prop <- ggplot(plot_dt[method != "OUTRIDER"], aes(method, prop, fill=method)) +
    facet_wrap(~ae_type) +
    geom_boxplot() +
    scale_fill_manual(values=c("lightblue", "purple4")) +
    labs(x="", y="Proportion of expression outliers\namong splicing outliers") +
    stat_compare_means(paired=TRUE, size=2) +
    theme_pubr() +
    theme(legend.position="none", 
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          axis.text.x=element_text(size=font_size-2)) +
    cowplot::background_grid(major="y", minor="y")
# g_prop

g_nr <- ggplot(plot_dt, aes(method, nr, fill=method)) +
    facet_wrap(~ae_type) +
    # geom_bar(stat="identity") +
    geom_boxplot() +
    scale_fill_manual(values=c("lightblue", "purple4", "darkgreen")) +
    labs(x="", y="Number of expression outliers\namong splicing outliers") +
    stat_compare_means(data=plot_dt[method != "OUTRIDER",], 
                       paired=TRUE, 
                       comparisons=list(c("FRASER", "FRASER2")),
                       size=2) +
    scale_y_log10() +
    annotation_logticks(sides="l") +
    theme_pubr() +
    theme(legend.position="none", 
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          axis.text.x=element_text(size=font_size-2)) +
    cowplot::background_grid(major="y", minor="y")
# g_nr

#+ combine into one figure
gg_sup <- ggarrange(g_prop,
                    g_nr, 
                    # nrow=2, ncol=1,
                    nrow=1, ncol=2,
                    labels=letters[1:2],
                    widths=c(1,1.25)
)
gg_sup

# TODO 
# - enrichment of stop/frameshift variants (higher in FRASER2?)
# - for venn diagram categories (F1 only, F2 only, intersection)
# - separately for over/under expression


#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=1.0*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=0.75*page_width, height=0.4*page_width, unit=width_unit, dpi=300)
