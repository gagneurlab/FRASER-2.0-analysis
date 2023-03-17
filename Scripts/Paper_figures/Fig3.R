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
#'   type: script
#'---

# #'     - reproducibility_all: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/tissue_reproducibility_FDR/all_tissue_reproducibility_ggplot.Rds"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
# library(gridExtra)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5
scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", 
                                              scales::scientific_format()(x))) }

# # reproducibility across tissues
# reprod_plots <- readRDS(snakemake@input$reproducibility_all)
# reprod_perc <- reprod_plots[["reprod_perc"]]
# dodge <- position_dodge(width=0.9)
# reprod_data <- as.data.table(reprod_perc$data)
# reprod_data <- reprod_data[snptype == "rare SpliceAI",]
# # reprod_data[Method == "FRASER\n+ Intron Jaccard Index", Method := "FRASER\n+ Intron\nJaccard Index"]
# # reprod_data[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER\n+ Intron\nJaccard Index", "FRASER2"))]
# reprod_data <- reprod_data[Method != "FRASER\n+ Intron Jaccard Index",]
# reprod_data[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER2"))]
# reprod_data <- reprod_data[, .(nr_outliers=sum(nr_outliers)), by="Method,nr_tissues"]
# reprod_data[, percentage := nr_outliers / sum(nr_outliers) , by="Method"]
# # conf interval
# reprod_data[, row_idx := 1:.N]
# reprod_data[, c("lower", "upper") := as.list(
#     prop.test(nr_outliers, 
#               reprod_data[Method == .SD[, unique(Method)], sum(nr_outliers)]
#     )$conf.int[1:2]), 
#     by="row_idx"]
# reprod_main <-  ggplot(reprod_data, #[pval_cutoff == "italic(P) < 10^-07"],
#                        aes(nr_tissues, percentage, fill=Method)) +
#     geom_col(position=dodge) +
#     scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
#     scale_y_continuous(labels = scales::percent) +
#     labs(x="Tissues in which outlier is present", y="Percentage of outliers") +
#     geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = 0.5) +
#     # guides(fill=guide_legend(nrow=2)) +
#     theme_pubr() + 
#     theme(
#         legend.position="right",
#         # legend.position=c(0.7,0.8),
#         # legend.spacing.x=unit(0.25, "cm"),
#         # legend.title=element_text(size=font_size),
#         legend.title=element_blank(),
#         # legend.text=element_text(size=font_size-1),
#         axis.title=element_text(face="bold"),
#         text=element_text(size=font_size)) + 
#     cowplot::background_grid(major="y", minor="y") 
# # reprod_main
# reprod_sup <- reprod_plots[["reprod_both"]] +
#     labs(x="Tissues in which outlier is present", y="") +
#     scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) +
#     # guides(fill=guide_legend(nrow=2)) +
#     theme_pubr() + 
#     theme(
#         legend.position="right",
#         # legend.position=c(0.8,0.9),
#         # legend.spacing.x=unit(0.25, "cm"),
#         # legend.title=element_text(size=font_size),
#         legend.title=element_blank(),
#         legend.text=element_text(size=font_size),
#         axis.title=element_text(face="bold"),
#         text=element_text(size=font_size)) + 
#     cowplot::background_grid(major="y", minor="y") 
# # reprod_sup

# sequencing depth correlation plots
seq_depth_cor_tissue_plot <- readRDS(snakemake@input$seq_depth_cor_tissue)
seq_depth_dt <- seq_depth_cor_tissue_plot$data
seq_depth_dt[method == "FRASER2", method := "FRASER 2.0"]
seq_depth_dt[, method := factor(method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER 2.0"))]
seq_depth_cor_tissue_plot <- ggscatter(seq_depth_dt, 
                                       x="record_count", y="nr_outliers_plus_one", 
                                       size=point_size,
                                       add="reg.line", 
                                       conf.int=TRUE, 
                                       add.params=list(color="blue", fill="lightgray"), 
                                       cor.coef=TRUE, 
                                       cor.method="spearman",
                                       cor.coef.size=3) + 
    # facet_wrap(~ method, scales="free_y") +
    facet_wrap(~ method, nrow=2) +
    scale_y_log10(limits=c(1, 10000)) + 
    labs(x="Mapped reads", y="Splicing outliers + 1") + 
    # seq_depth_cor_tissue_plot <- seq_depth_cor_tissue_plot + 
    scale_x_continuous(breaks=c(7.5e7, 1.5e8, 2.25e8), labels=scientific_10(c(7.5e7, 1.5e8, 2.25e8))) +
    # annotation_logticks(sides="l") +
    # labs(x="Mapped reads") +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# seq_depth_cor_tissue_plot$layers[[1]]$aes_params$size <- point_size
seq_depth_cor_tissue_plot

seq_depth_cor_all_plot <- readRDS(snakemake@input$seq_depth_cor_all)
# seq_depth_cor_all_plot <- seq_depth_cor_all_plot +
#     labs(x="Correlation coefficent\nFRASER", y="Correlation coefficient\nFRASER 2.0") +
#     theme_pubr() +
#     theme(axis.title=element_text(face="bold"),
#           text=element_text(size=font_size)) +
#     cowplot::background_grid(major="xy", minor="xy")
# seq_depth_cor_all_plot$layers[[1]]$aes_params$size <- point_size
cor_all_dt <- seq_depth_cor_all_plot$data
cor_all_melt <- melt(cor_all_dt, id.vars="tissue", value.name="cor_coef", variable.name="method")
cor_all_melt[method == "FRASER2", method := "FRASER 2.0"]
cor_all_melt[, method := factor(method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER 2.0"))]
seq_depth_cor_all_plot <- ggplot(cor_all_melt[tissue != "Minor_Salivary_Gland"], aes(method, cor_coef, fill=method)) +
    geom_boxplot(outlier.size=point_size) +
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4", "violetred")) +
    stat_compare_means( 
                       # aes(label = paste0("p = ", after_stat(p.format))),
                       paired=TRUE,
                       comparisons=list( c("FRASER 2.0", "LeafcutterMD"), c("FRASER 2.0", "SPOT"), c("FRASER 2.0", "FRASER") ),
                       step.increase=0.2, 
                       vjust=-0.1
                       ) +
    labs(x="", y="Correlation coefficent (R)") +
    ylim(-0.5,1.5) +
    theme_pubr() +
    theme(legend.position="none",
        axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")
seq_depth_cor_all_plot

#+ combine panels into figure, width=12, height=20
# row1 <- ggarrange(
#                 seq_depth_cor_tissue_plot,
#                   labels=LETTERS[1],
#                   nrow=1, ncol=1,
#                   # widths=c(3,1),
#                   align="hv")
# row2 <- ggarrange(nr_outliers_per_sample_plot, 
#                   labels=LETTERS[2], 
#                   nrow=1, ncol=1)
row2 <- ggarrange(seq_depth_cor_tissue_plot,
                  seq_depth_cor_all_plot,
                  labels=LETTERS[1:2],
                  align="hv",
                  nrow=2, ncol=1,
                  heights=c(1.75,1)
                  # widths=c(2,1)
                  )
gg_figure <- row2
# gg_figure <- ggarrange(row1, row2, nrow=2)
# gg_figure <- ggarrange(
#     row1,
#     row2, 
#     nrow=2, ncol=1,
#     align="hv", 
#     heights=c(1,1))
gg_figure

#+ save figure as png and pdf
# ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=15, height=13)
# ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=15, height=13)
# ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=0.75*page_width, unit=width_unit)
# ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=0.75*page_width, unit=width_unit)
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=0.7*page_width, height=1*page_width, unit=width_unit)
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=0.7*page_width, height=1*page_width, unit=width_unit)
