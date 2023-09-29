#'---
#' title: Paper figure 4 (rare disease analysis)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/fig4.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - prokisch_f2_vs_f1_plots: '`sm config["mito_processed_results"] + "/ggplot_rds/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/FRASER2_vs_FRASER1_old_filter_ggplots.Rds"`'
#'     - prokisch_fdr_comparison: '`sm config["mito_processed_results"] + "/ggplot_rds/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/FDR_subset_MAF0.001_no_utr_ggplots_F1_old_filter.Rds"`'
#'     - prokisch_omim_res_table: '`sm config["mito_processed_results"] + "/results/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/deltaJaccard0.1/results_gene_FDRomim.tsv"`'
#'     - udn_nr_outliers_comparison: '`sm config["DATADIR"] + "/udn/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/nrOutliers_comparison_ggplot.Rds"`'
#'     - power_analysis_full_res_all: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/combined_results_full.tsv"`'
#'     - power_analysis_full_res_patho: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_full.tsv"`'
#'     - patho_sample_anno: '`sm config["mito_sample_anno"]`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig4.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/Fig4.pdf"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/Fig4.tiff"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggupset)
library(ggbeeswarm)
library(data.table)
library(tidyverse)
library(magrittr)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5

#+ read in plots and data for the different panels
f2_vs_f1_plots <- readRDS(snakemake@input$prokisch_f2_vs_f1_plots)
fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
udn_plots <- readRDS(snakemake@input$udn_nr_outliers_comparison)

# boxplot of number of outliers in different cohorts
mito_dt <- fdr_comparison_plots[["boxplots_numOut_gene"]]$data
mito_dt[, dataset_group:="Yépez et al. dataset"] # "Mitochondrial disease cohort"]
mito_dt[, dataset_label:="Fibroblasts"]
mito_dt
mito_dt_omim <- fread(snakemake@input$prokisch_omim_res_table)
mito_dt_omim <- mito_dt_omim[!duplicated(mito_dt_omim, by=c("sampleID","seqnames","start","end","strand")), ]
# get nr of outliers per sample (gene-level)
mito_dt_omim[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
mito_dt_omim[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
# add name of dataset
mito_dt_omim[, dataset_group:="Yépez et al. dataset"] # "Mitochondrial disease cohort"]
mito_dt_omim[, dataset_label:="Fibroblasts"]
mito_dt_omim[, method:="FRASER2 (OMIM)"]
mito_dt_omim <- mito_dt_omim[, .(num_out = unique(numOutliersPerSample)), by="sampleID,method,dataset_label,dataset_group"]
mito_dt_omim

mito_dt <- rbind(mito_dt, mito_dt_omim)
dcast_dt <- dcast(mito_dt, sampleID + dataset_label + dataset_group ~ method, value.var="num_out")
dcast_dt[is.na(`FRASER2 (OMIM)`), `FRASER2 (OMIM)` := 0]
mito_dt <- melt(dcast_dt, id.vars=c("sampleID", "dataset_label","dataset_group"), value.name="num_out", variable.name="method")
mito_dt

udn_dt <- udn_plots[["boxplots_numOut_gene"]]$data
udn_dt[, dataset_group := "Undiagnosed Disease Network (UDN)"]
udn_dt[, dataset_label := NULL]
setnames(udn_dt, "dataset_name", "dataset_label")

num_out_dt <- rbind(mito_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)], 
                    udn_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)])
tissue_order <- num_out_dt[dataset_group == "Undiagnosed Disease Network (UDN)", median(num_out),by="dataset_label,method"][method == "FRASER"][order(-V1), dataset_label]
num_out_dt[, dataset_label := factor(dataset_label, levels=tissue_order)]
num_out_dt[method == "FRASER2", method := "FRASER 2.0"]
num_out_dt[method == "FRASER2 (OMIM)", method := "FRASER 2.0\n(OMIM)"]
num_out_dt[method == "FRASER2\n(OMIM + RV)", method := "FRASER 2.0\n(OMIM + RV)"]
num_out_dt[, method := factor(method, levels=c("FRASER", "FRASER 2.0", "FRASER 2.0\n(OMIM)", "FRASER 2.0\n(OMIM + RV)"))]
num_out_dt[, dataset_group := factor(dataset_group, levels=c("Undiagnosed Disease Network (UDN)", "Yépez et al. dataset"))]

g_boxplot_all <- ggplot(num_out_dt, aes(dataset_label, num_out+1, col=method)) +
    facet_grid(~dataset_group, scales="free_x", space="free_x") +
    geom_boxplot(outlier.size=point_size) +
    labs(x="", y="Splicing outliers\nper sample + 1") +
    scale_y_log10() +
    scale_color_manual(values=c("dodgerblue3", "purple4", "purple1", "violetred")) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme(legend.title = element_blank()) +
    cowplot::background_grid(major="y", minor="y") #+
    # annotate("segment", x=0, xend=0, y=0, yend=Inf, linewidth = 1) # add y axis line so that logticks dont look weird in middle facets
a <- annotation_logticks(sides='l')
a$data <- data.frame(x=NA, dataset_group=levels(num_out_dt[, dataset_group])[1]) # have logticks only on left facet (changes facet order :( ))
g_boxplot_all <- g_boxplot_all + a
g_boxplot_all

# fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
upset_dt <- fdr_comparison_plots[["gg_upset"]]$data
upset_dt[, tmp := sapply(comb, FUN=function(x) paste(unlist(x), collapse="-"))]
set_order <- c(
    "FRASER-FRASER2-FRASER2 (OMIM + RV)-pathogenic splice defect",
    "FRASER-FRASER2-pathogenic splice defect",
    "FRASER-pathogenic splice defect",
    "FRASER-FRASER2 (OMIM + RV)-pathogenic splice defect",
    "FRASER2-FRASER2 (OMIM + RV)-pathogenic splice defect",
    "FRASER-FRASER2 (OMIM + RV)",
    "FRASER2-FRASER2 (OMIM + RV)",
    "FRASER-FRASER2-FRASER2 (OMIM + RV)",
    "FRASER2 (OMIM + RV)",
    "FRASER-FRASER2",
    "FRASER2",
    "FRASER"
)
upset_plot <- ggplot(upset_dt, aes(x=comb)) +
    geom_bar() +
    theme_pubr() +
    scale_x_mergelist(sep = "-", limits=set_order) +
    axis_combmatrix(sep = "-", override_plotting_function = function(df){
        # print(df)
        df %>%
            mutate(patho_set = case_when(
                ! observed ~ "not observed",
                map_lgl(labels_split, ~ "pathogenic splice defect" %in% .x) ~ "Patho",
                observed ~ "Non-patho"
            )) %>%
            mutate(single_label = fct_recode(single_label, "Pathogenic events" = "pathogenic splice defect"))  %>%
            mutate(single_label = fct_recode(single_label, "FRASER 2.0" = "FRASER2"))  %>%
            mutate(single_label = fct_recode(single_label, "FRASER 2.0 (OMIM + RV)" = "FRASER2 (OMIM + RV)"))  %>%
            mutate(single_label = factor(single_label, levels=rev(c("Pathogenic events", "FRASER 2.0", "FRASER 2.0 (OMIM + RV)", "FRASER")), ordered=TRUE)) %>%
            ggplot(aes(x = at, y = single_label)) +
            geom_rect(aes(fill = index %% 2 == 0), ymin=df$index-0.5, ymax=df$index+0.5, xmin=0, xmax=1) +
            geom_point(aes(color = patho_set), size = 3) +
            geom_line(data= function(dat) dat[dat$observed, ,drop=FALSE], aes(group = labels, color = patho_set), size= 1.2) +
            ylab("") + xlab("") +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_fill_manual(values= c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
            scale_color_manual(values= c("Patho" = "firebrick", "Non-patho" = "black", "not observed" = "lightgrey")) +
            guides(fill="none") +
            theme(
                legend.position = "none",
                panel.background = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.length = unit(0, "pt"),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.line = element_blank(),
                panel.border = element_blank()
            )
    }) + 
    labs(x="", y="Event count") +
    geom_text(aes(label = ..count..), stat = 'count', nudge_y = 500, 
              size=font_size/.pt) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme_combmatrix() +
    # theme_combmatrix(combmatrix.label.width=unit(3, "cm")) +
    theme(axis.title.y=element_text(vjust=-45)) +
    cowplot::background_grid(major="y", minor="y")
upset_plot

#+ power analysis results
# pathogenic cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj),]
# patho_sa <- patho_sa[!grepl("deletion|cnv", VARIANT_EFFECT),]
# patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] 
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "C19orf70"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]

# transcriptome-wide results
res_tps <- fread(snakemake@input$power_analysis_full_res_patho)
# dont show sizes 70 and 90 in final plot
res_tps <- res_tps[!size %in% c(70, 90)]
res_tps[, dataset := "Yepez et al. samples only"]

# plot proportion of recovered pathogenic cases
gprop <- ggplot(res_tps[padjustGene <= 0.1, .(Noutliers=.N, prop=.N/length(patho_tmp)), by=.(size, sim, FDR_set, dataset)], 
                aes(as.factor(size), prop)) + 
    geom_beeswarm(groupOnX=TRUE, size=point_size) +  
    labs(
        y = paste0('Fraction of recovered\n splicing outliers (n=', length(patho_tmp), ')'),
        x = 'Sample size') +
    scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.9)) + 
    cowplot::background_grid(major="xy", minor="xy")
# gprop

#+ combine panels into figure, width=15, height=12
gg_row2 <- ggarrange(
    upset_plot,
    gprop,
    nrow=1, ncol=2,
    labels=LETTERS[2:3],
    font.label = list(size = 12, color = "black", face = "bold"),#, family = "Arial"),
    widths=c(1.9,1)
)
gg_figure <- ggarrange(g_boxplot_all,
                       gg_row2,
                       nrow=2, ncol=1,
                       heights=c(1,1), 
                       labels=c(LETTERS[1], ""),
                       font.label = list(size = 12, color = "black", face = "bold"))#, family = "Arial"))
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outTiff, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300) 
