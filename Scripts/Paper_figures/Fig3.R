#'---
#' title: Paper figure 3 (GTEx analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/fig3.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - seq_depth_cor_tissue: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/seq_depth_cor_single_ggplot.Rds"`'
#'     - seq_depth_cor_all: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/seq_depth_cor_all_ggplot.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig3.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/Fig3.pdf"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/Fig3.tiff"`'
#'   type: script
#'---


saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
source("src/R/ggplot_theme_for_manuscript.R") 

# library(extrafont)
# font_import(prompt=FALSE)
# loadfonts(device="all")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5
scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", 
                                              scales::scientific_format()(x))) }

# sequencing depth correlation plots
seq_depth_cor_tissue_plot <- readRDS(snakemake@input$seq_depth_cor_tissue)
seq_depth_dt <- seq_depth_cor_tissue_plot$data
seq_depth_dt[method == "FRASER2", method := "FRASER 2.0"]
seq_depth_dt[method == "LeafcutterMD", method := "LeafCutterMD"]
seq_depth_dt[, method := factor(method, levels=c("LeafCutterMD", "SPOT", "FRASER", "FRASER 2.0"))]
seq_depth_dt[, cor_coef := round(cor(record_count, nr_outliers, method="spearman"), digits=2), by="method"]
seq_depth_cor_tissue_plot <- ggscatter(seq_depth_dt, 
                                       x="record_count", y="nr_outliers_plus_one", 
                                       size=point_size,
                                       add="reg.line", 
                                       conf.int=TRUE, 
                                       add.params=list(color="blue", fill="lightgray"), 
                                       cor.coef=FALSE, # TRUE
                                       cor.method="spearman",
                                       cor.coef.size=3) + 
    facet_wrap(~ method, nrow=2) +
    geom_text(data=seq_depth_dt[, .(record_count=2.25*1e8, 
                                    nr_outliers_plus_one = 6500, 
                                    cor_coef = paste0("rho == ", unique(cor_coef)) ), by="method"], 
              mapping=aes(label=cor_coef), parse=TRUE, size = font_size/.pt) +
    scale_y_log10(limits=c(1, 10000)) + 
    labs(x="Mapped reads", y="Splicing outliers + 1") + 
    scale_x_continuous(breaks=c(7.5e7, 1.5e8, 2.25e8), labels=scientific_10(c(7.5e7, 1.5e8, 2.25e8))) +
    # annotation_logticks(sides="l") +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    cowplot::background_grid(major="xy", minor="xy")
seq_depth_cor_tissue_plot

seq_depth_cor_all_plot <- readRDS(snakemake@input$seq_depth_cor_all)
cor_all_dt <- seq_depth_cor_all_plot$data
cor_all_melt <- melt(cor_all_dt, id.vars="tissue", value.name="cor_coef", variable.name="method")
cor_all_melt[method == "FRASER2", method := "FRASER 2.0"]
cor_all_melt[method == "LeafcutterMD", method := "LeafCutterMD"]
cor_all_melt[, method := factor(method, levels=c("LeafCutterMD", "SPOT", "FRASER", "FRASER 2.0"))]
seq_depth_cor_all_plot <- ggplot(cor_all_melt[tissue != "Minor_Salivary_Gland"], aes(method, cor_coef, fill=method)) +
    geom_boxplot(outlier.size=point_size) +
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4", "violetred")) +
    stat_compare_means( 
                       # aes(label = paste0("p = ", after_stat(p.format))),
                       paired=TRUE,
                       comparisons=list( c("FRASER 2.0", "LeafCutterMD"), c("FRASER 2.0", "SPOT"), c("FRASER 2.0", "FRASER") ),
                       step.increase=0.2, 
                       vjust=-0.1,
                       size = font_size/.pt
                       ) +
    labs(x="", y="Spearman correlation\ncoefficent") +
    ylim(-0.5,1.5) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme(legend.position="none") + 
    cowplot::background_grid(major="y", minor="y")
seq_depth_cor_all_plot

#+ combine panels into figure, width=12, height=20
gg_figure <- ggarrange(seq_depth_cor_tissue_plot,
                  seq_depth_cor_all_plot,
                  labels=LETTERS[1:2],
                  font.label = list(size = 12, color = "black", face = "bold"),
                  align="hv",
                  nrow=2, ncol=1,
                  heights=c(1.75,1)
                  )
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=114, height=1*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=114, height=1*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_figure, filename=snakemake@output$outTiff, width=114, height=1*page_width, unit=width_unit, dpi=300)
# width=page_width
