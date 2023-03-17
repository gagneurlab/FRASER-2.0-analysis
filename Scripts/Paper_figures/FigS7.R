#'---
#' title: Paper figure S7 (GTEx Venn)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS7.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - comb_outliers_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/optQ/delta0.1/combined_outliers_venn.Rds"`'
#'     - pr_curves_venn: '`sm config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_vs_FRASER2_venn_rv_recall_plots_rareSpliceAI.Rds"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS7.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS7.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(ggplot2)
library(ggpubr)
library(cowplot)
library(eulerr)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size - 1
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in outliers for venn plotting
all_outliers <- readRDS(snakemake@input$comb_outliers_rds)

# #+ create ggvenn
# library(ggvenn)
# venn_dt <- data.table(outlier_id = unique(unlist(all_outliers)))
# venn_dt[, FRASER := outlier_id %in% all_outliers$FRASER]
# venn_dt[, FRASER2 := outlier_id %in% all_outliers$FRASER2]
# ggplot(venn_dt) +
#     geom_venn(aes(A=FRASER, B=FRASER2), 
#               fill_color = c("lightblue", "purple4"), 
#               stroke_size = 0.5,
#               text_size=4) +
#     coord_fixed() +
#     theme_void()

#+ create proportional venn diagram
venn_plot <- plot(euler(all_outliers), 
     quantities = list(TRUE, type=c("counts", "percent"), fontsize=font_size),
     labels=list(fontsize=font_size),
     fills=c("white", "white"),
     edges=c("black", "black")
     # fills=c("lightblue", "purple4")
     )
venn_plot


#+ enrichment of splice variants for venn diagram categories
pr_curves <- readRDS(snakemake@input$pr_curves_venn)
maxRank <- 50000
pr_curve_venn <- pr_curves[[paste0('recall_n=', maxRank)]] + 
    labs(title="", y="Recall of rare\nsplice affecting variants") +
    # guides(col=guide_legend(order=2, nrow=1)) +
    theme_pubr() + 
    theme(legend.position="bottom", 
          legend.box="vertical",
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size))  + 
    cowplot::background_grid(major="xy", minor="xy") +
    guides(col=guide_legend(nrow=3, order=1))
pr_curve_venn

#+ compile figure
gg_figure <- ggarrange(
            # venn_plot,
            pr_curve_venn
            # nrow=1, ncol=2,
            # widths=c(1,1),
            # labels=LETTERS[1]
          )
# gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=0.5*page_width, height=0.6*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=0.5*page_width, height=0.6*page_width, unit=width_unit, dpi=300) 
