#'---
#' title: Paper figure 4 (rare disease analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_power_analysis_withUDN.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - power_analysis_full_res_all: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/combined_results_full.tsv"`'
#'     - power_analysis_full_res_patho: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_full.tsv"`'
#'     - power_analysis_full_res_all_withUDN: '`sm config["DATADIR"] + "/power_analysis/mito_with_extUDN/processed_results/aberrant_splicing/combined_results_full.tsv"`'
#'     - power_analysis_full_res_patho_withUDN: '`sm config["DATADIR"] + "/power_analysis/mito_with_extUDN/processed_results/aberrant_splicing/patho_results_full.tsv"`'
#'     - patho_sample_anno: '`sm config["mito_sample_anno"]`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_power_analysis_withUDN.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_power_analysis_withUDN.pdf"`'
#'   type: script
#'---

# #'     - power_analysis_gtex: '`sm config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/combined_results.tsv"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm) 
library(data.table)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5

#+ pathogenic cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj),]
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "C19orf70"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]

#+ power analysis results (transcriptome_wide)
res_all_prokisch <- fread(snakemake@input$power_analysis_full_res_all)
res_tps_prokisch <- fread(snakemake@input$power_analysis_full_res_patho)
res_all_withUDN <- fread(snakemake@input$power_analysis_full_res_all_withUDN)
res_tps_withUDN <- fread(snakemake@input$power_analysis_full_res_patho_withUDN)

res_all_prokisch[, dataset := "Yépez et al. samples only"]
res_tps_prokisch[, dataset := "Yépez et al. samples only"]
res_all_withUDN[, dataset := "Yépez et al. samples with external UDN samples"]
res_tps_withUDN[, dataset := "Yépez et al. samples with external UDN samples"]
res_all <- rbind(res_all_prokisch, res_all_withUDN)
res_tps <- rbind(res_tps_prokisch, res_tps_withUDN)
res_all[, dataset := factor(dataset, levels=c(unique(res_all_prokisch$dataset), unique(res_all_withUDN$dataset)))]
res_tps[, dataset := factor(dataset, levels=c(unique(res_tps_prokisch$dataset), unique(res_tps_withUDN$dataset)))]

# dont show sizes 70 and 90 in final plot
res_all <- res_all[!size %in% c(70, 90)]
res_tps <- res_tps[!size %in% c(70, 90)]

# plot proportion of recovered pathogenic cases
gprop <- ggplot(res_tps[padjustGene <= 0.1, .(Noutliers=.N, prop=.N/length(patho_tmp)), by=.(size, sim, FDR_set, dataset)], 
                aes(as.factor(size), prop, color=dataset)) +
    geom_beeswarm(groupOnX=TRUE, size=point_size, dodge.width=0.75) +  
    labs(
        y = paste0('Fraction of recovered\n splicing outliers (N=', length(patho_tmp), ')'),
        x = 'Sample size') +
    scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
    scale_color_brewer(palette="Set1") +
    guides(color=guide_legend(title="", nrow=2)) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.9),
          legend.position=c(0.51, 0.18)) + 
          # legend.background=element_blank()) + 
    cowplot::background_grid(major="xy", minor="xy")
gprop

# plot p values
gpval <- ggplot(res_tps_withUDN, aes(as.factor(size), -log10(pValue))) +
    geom_violin() + 
    geom_boxplot(width=0.1, alpha=0.2) +
    labs(
        y = expression(bold(-log[10](italic(P)~"value"))), 
        x = 'Sample size', 
        col = 'Gene') + 
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position = "bottom") + 
    cowplot::background_grid(major="xy", minor="xy") +
    guides(col=guide_legend(nrow=3)) 
gpval

# plot # of outliers
gtotal <- ggplot(res_all_withUDN[padjustGene <= 0.1, .N, by=.(size=factor(size), sim, sampleID, dataset)][,.(med=median(N)),by=.(size, sim, dataset)], 
                 aes(size, med)) + 
    geom_beeswarm(size=0.5) +
    # ylim(0, 25) +
    labs(
        x = 'Sample size', 
        y = 'Median of splicing\noutliers per sample') +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.9)) + 
    cowplot::background_grid(major="y", minor="y")
gtotal

nr_out_both <- merge(res_all_prokisch[size == 300, uniqueN(hgncSymbol), by="sampleID,sim"][,median(V1), by="sampleID"][sampleID %in% patho_sa$RNA_ID],
                     res_all_withUDN[size == 300, uniqueN(hgncSymbol), by="sampleID,sim"][,median(V1), by="sampleID"][sampleID %in% patho_sa$RNA_ID],
                     by="sampleID")
nr_out_both

g_nr_out_comp <- ggplot(nr_out_both, aes(V1.x, V1.y)) +
    geom_point(size=point_size) +
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    xlim(0, nr_out_both[, max(V1.x)]) + 
    ylim(0, nr_out_both[, max(V1.y)]) + 
    labs(x="Outliers per sample\n(combined with Yépez et al. samples only)",
         y="Outliers per sample\n(combined with UDN samples)") +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    cowplot::background_grid(major="xy", minor="xy")
g_nr_out_comp

#+ combine panels into figure, width=15, height=12
gg <- ggarrange(
    gprop,
    g_nr_out_comp,
    gtotal,
    gpval,
    # nrow=1, ncol=2,
    nrow=2, ncol=2,
    widths=c(4,5),
    heights=c(1, 0.9),
    labels=LETTERS[1:4],
    font.label=list(size=12, color = "black", face = "bold", family = font),
    align="hv"
)
gg

#+ save figure as png and pdf
ggsave(plot=gg, filename=snakemake@output$outPng, width=page_width, height=0.9*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg, filename=snakemake@output$outPdf, width=page_width, height=0.9*page_width, unit=width_unit)