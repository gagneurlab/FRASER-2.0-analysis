#'---
#' title: Paper figure 4 (rare disease analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS9.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - power_analysis_full_res_all: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/combined_results_full.tsv"`'
#'     - power_analysis_full_res_patho: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_full.tsv"`'
#'     - power_analysis_gtex: '`sm config["DATADIR"] + "/power_analysis/GTEx_v8/processed_results/aberrant_splicing/combined_results.tsv"`'
#'     - patho_sample_anno: '`sm config["mito_sample_anno"]`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS9.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS9.pdf"`'
#'   type: script
#'---

# #'     - power_analysis_omimRv_res_all: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/combined_results_OMIM_RV.tsv"`'
# #'     - power_analysis_omimRv_res_patho: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_OMIM_RV.tsv"`'
saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm) 
library(data.table)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ pathogenic cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj) & !grepl("deletion|cnv", VARIANT_EFFECT),]
# patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] 
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "C19orf70"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]

#+ power analysis results (transcriptome_wide)
res_all <- fread(snakemake@input$power_analysis_full_res_all)
res_tps <- fread(snakemake@input$power_analysis_full_res_patho)
res_gtex <- fread(snakemake@input$power_analysis_gtex)

# dont show sizes 70 and 90 in final plot
res_all <- res_all[!size %in% c(70, 90)]
res_tps <- res_tps[!size %in% c(70, 90)]

# plot p values
gpval <- ggplot(res_tps, aes(as.factor(size), -log10(pValue))) +
    geom_violin() + 
    geom_boxplot(width=0.1, alpha=0.2) +
    # geom_beeswarm(aes(col = gsub(".*_", "", tmp))) + 
    labs(
        y = expression(bold(-log[10](italic(P)~"value"))), 
        x = 'Sample size', 
        col = 'Gene') + 
    theme_pubr() +
    theme(legend.position = "bottom",
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy") +
    guides(col=guide_legend(nrow=3)) 
# gpval
# plot delta jaccard
gdelta <- ggplot(res_tps, aes(size, abs(deltaPsi))) +
    facet_wrap(~tmp) +
    geom_point(size=0.5) +
    geom_line(data=res_tps[, median(abs(deltaPsi)), by="tmp,size"], aes(size, V1)) +
    labs(
        y = expression(bquote(Delta(Intron~Jaccard~Index))),
        x = 'Sample size',
        col = 'Gene') +
    theme_pubr() +
    theme(legend.position = "bottom",
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) +
    cowplot::background_grid(major="xy", minor="xy")
gdelta <- ggplot(res_tps[, median(deltaPsi), by="tmp,size"], aes(size, V1, col=tmp)) + 
    geom_line() +
    geom_point(data=res_tps, aes(size, deltaPsi, col=tmp), size=0.5) +
    labs(
        # y = expression(paste("|", Delta, "(Intron Jaccard Index)", "|")),
        y = expression(paste(Delta, "(Intron Jaccard Index)")),
        x = 'Sample size',
        col = '') +
    theme_pubr() +
    theme(legend.position = "bottom",
        axis.title=element_text(face="bold"),
        text=element_text(size=font_size)) +
    cowplot::background_grid(major="xy", minor="xy")
# gdelta
# plot # of outliers
gtotal <- ggplot(res_all[padjustGene <= 0.1, .N, by=.(size=factor(size), sim, sampleID)][,.(med=median(N)),by=.(size, sim)], 
                 aes(size, med)) + 
    geom_beeswarm(size=0.5) +
    ylim(0, 15) + 
    labs(
        x = 'Sample size', 
        y = 'Median of splicing\noutliers per sample') +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")
gtotal

res_gtex[, aberrant_0_0 := padjustGene <= 0.1 & abs(deltaPsi) >= 0]
res_gtex[, aberrant_0_1 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.1]
res_gtex[, aberrant_0_2 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.2]
res_gtex[, aberrant_0_3 := padjustGene <= 0.1 & abs(deltaPsi) >= 0.3]
res_gtex_melt <- melt(res_gtex, 
                      id.vars=c("sampleID","hgncSymbol","size","sim"), 
                      measure.vars=c("aberrant_0_0", "aberrant_0_1", "aberrant_0_2", "aberrant_0_3"),
                      variable.name="dJ", 
                      value.name="aberrantStatus")
res_gtex_plot <- res_gtex_melt[, .SD[aberrantStatus == TRUE, .N], by=.(size=factor(size), sim, sampleID, dJ)]
res_gtex_plot <- res_gtex_plot[,.(med=median(V1)),by=.(size, sim, dJ)]
res_gtex_plot[, dJ := gsub("_", ".", gsub("aberrant_", "", dJ))]
gtotal_gtex <- ggplot(res_gtex_plot, aes(size, med, col=dJ)) + 
    geom_beeswarm(size=0.5) +
    ylim(0, 15) + 
    labs(
        x = 'Sample size', 
        y = 'Median of splicing\noutliers per sample') +
    guides(color=guide_legend(title=bquote(Delta~J))) +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")
# gtotal_gtex

#+ combine panels into figure, width=15, height=12
gg <- ggarrange(
                gtotal,
                gpval,
                nrow=1, ncol=2,
                widths=c(3,4),
                labels=LETTERS[1:2],
                common.legend=TRUE, legend = "bottom"
            )
# gg

#+ save figure as png and pdf
ggsave(plot=gg, filename=snakemake@output$outPng, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg, filename=snakemake@output$outPdf, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300) 
