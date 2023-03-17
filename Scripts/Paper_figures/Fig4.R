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
#'     - patho_sample_anno: '`sm config["mito_sample_anno"]'`
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig4.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/Fig4.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
# library(ggvenn)
library(ggupset)
library(ggbeeswarm)
library(data.table)
library(tidyverse)
library(magrittr)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5

#+ read in plots and data for the different panels
f2_vs_f1_plots <- readRDS(snakemake@input$prokisch_f2_vs_f1_plots)
fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
udn_plots <- readRDS(snakemake@input$udn_nr_outliers_comparison)


# g_venn_wo_F1 <- fdr_comparison_plots[["FDR_venn_wo_F1"]] +
#     theme(axis.title=element_text(face="bold"),
#           text=element_text(size=font_size))

#+ create proportional venn diagram
# venn_dt <- fdr_comparison_plots[["FDR_venn_wo_F1"]]$data
# venn_ls <- list(FRASER2 = venn_dt[FRASER2 == TRUE, values],
#                 "FRASER2 (OMIM + RV)" = venn_dt[`FRASER2 (OMIM + RV)` == TRUE , values],
#                 "Pathogenic events" = venn_dt[`pathogenic
# splice defect` == TRUE, values])
# library(eulerr)
# venn_plot <- plot(euler(venn_ls), 
#                   quantities = list(TRUE, type=c("counts", "percent"), fontsize=font_size),
#                   labels=list(fontsize=font_size)#,
#                   # fills=c("lightblue", "purple4")
# )
# venn_plot

# # draw venn manual for main figure plot
# # get set sizes
# venn_dt <- g_venn_wo_F1$data
# venn_dt[, FRASER := NULL]
# area_dt <- melt(venn_dt, id.vars="values", variable.name="Set", value.name="present")
# area_dt <- area_dt[present == TRUE,]
# area_dt[, comb:=paste(Set, collapse=","), by="values"]
# area_dt <- area_dt[!duplicated(values)]
# area_sizes <- area_dt[, table(comb)]
# area_labels <- as.data.table(area_sizes)
# area_labels[, x:=c(0, 2.5, 1.5, 1.5, 2.25, 0.75)]
# area_labels[, y:=c(0.2, 0.2, 0.2, -0.7, -0.7, -0.7)]
# 
# gen_circle <- function(group, x_offset = 0, y_offset = 0, radius = 1,
#                        radius_b = radius, theta_offset = 0, length.out = 100) {
#     tibble(group = group,
#            theta = seq(0, 2 * pi, length.out = length.out)) %>%
#         mutate(x_raw = radius * cos(theta),
#                y_raw = radius_b * sin(theta),
#                x = x_offset + x_raw * cos(theta_offset) - y_raw * sin(theta_offset),
#                y = y_offset + x_raw * sin(theta_offset) + y_raw * cos(theta_offset))
# }
# 
# g_venn_wo_F1_manual <- ggplot() +
#     geom_polygon(aes(x = x, y = y),
#                  data = gen_circle(group = 0L, x_offset = 0, y_offset = 0, radius = 2),
#                  color = "purple4", fill = "purple4", alpha = .5) +
#     geom_polygon(aes(x = x, y = y),
#                  data = gen_circle(group = 0L, x_offset = 2, y_offset = 0, radius = 1, radius_b = 1),
#                  color = "purple1" , fill = "purple1", alpha = .5) +
#     geom_polygon(aes(x = x, y = y),
#                  data = gen_circle(group = 0L, x_offset = 1.5, y_offset = -0.7, radius = 1.2, radius_b=0.3),
#                  color = "darkgreen" , fill = "darkgreen", alpha = .5) +
#     geom_text(data=area_labels, aes(label=N, x=x, y=y), size=3) +
#     geom_text(data=area_labels[c(1, 2,6)], aes(label=comb, x=x-1, y=ifelse(y>0, y+1, y-1)), size=3) +
#     theme_nothing() 
# # g_venn_wo_F1_manual

boxplot_prokisch <- fdr_comparison_plots[["boxplots_numOut_gene"]] + 
    labs(y="Splicing outliers\nper sample + 1",
         title="Mitochondrial disease cohort") +
    scale_color_manual(values=c("dodgerblue3", "purple4", "purple1")) +
    theme(axis.title=element_text(face="bold"),
          plot.title=element_text(face="bold", size=font_size),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")

udn_dt <- dcast(udn_plots[["boxplots_numOut_gene"]]$data[, median(num_out), by="method,dataset_name"], dataset_name ~ method)[, fc := FRASER/FRASER2]
udn_dt
g_boxplot_udn <- udn_plots[["boxplots_numOut_gene"]] +
    labs(y="Splicing outliers\nper sample + 1",
         title="Undiagnosed Disease Network (UDN)") +
    scale_color_manual(values=c("dodgerblue3", "purple4", "purple1")) +
    theme(axis.title=element_text(face="bold"),
          plot.title=element_text(face="bold", size=font_size),
          text=element_text(size=font_size)) +
    cowplot::background_grid(major="y", minor="y")

mito_dt <- fdr_comparison_plots[["boxplots_numOut_gene"]]$data
mito_dt[, dataset_group:="Mitochondrial disease cohort"]
mito_dt[, dataset_label:=paste0("Fibroblasts (N=", uniqueN(sampleID),")")]
mito_dt
mito_dt_omim <- fread(snakemake@input$prokisch_omim_res_table)
mito_dt_omim <- mito_dt_omim[!duplicated(mito_dt_omim, by=c("sampleID","seqnames","start","end","strand")), ]
# get nr of outliers per sample (gene-level)
mito_dt_omim[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
mito_dt_omim[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
# add name of dataset
mito_dt_omim[, dataset_group:="Mitochondrial disease cohort"]
mito_dt_omim[, dataset_label:=paste0("Fibroblasts (N=", uniqueN(mito_dt$sampleID),")")]
mito_dt_omim[, method:="FRASER2 (OMIM)"]
mito_dt_omim <- mito_dt_omim[, .(num_out = unique(numOutliersPerSample)), by="sampleID,method,dataset_label,dataset_group"]
mito_dt_omim

mito_dt <- rbind(mito_dt, mito_dt_omim)
dcast_dt <- dcast(mito_dt, sampleID + dataset_label + dataset_group ~ method, value.var="num_out")
dcast_dt[is.na(`FRASER2 (OMIM)`), `FRASER2 (OMIM)` := 0]
mito_dt <- melt(dcast_dt, id.vars=c("sampleID", "dataset_label","dataset_group"), value.name="num_out", variable.name="method")

udn_dt <- udn_plots[["boxplots_numOut_gene"]]$data
udn_dt[, dataset_group := "Undiagnosed Disease Network (UDN)"]

num_out_dt <- rbind(mito_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)], 
                    udn_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)])
tissue_order <- num_out_dt[,median(num_out),by="dataset_label,method"][method == "FRASER"][order(-V1), dataset_label]
num_out_dt[, dataset_label := factor(dataset_label, levels=tissue_order)]
num_out_dt[method == "FRASER2", method := "FRASER 2.0"]
num_out_dt[method == "FRASER2 (OMIM)", method := "FRASER 2.0\n(OMIM)"]
num_out_dt[method == "FRASER2\n(OMIM + RV)", method := "FRASER 2.0\n(OMIM + RV)"]
num_out_dt[, method := factor(method, levels=c("FRASER", "FRASER 2.0", "FRASER 2.0\n(OMIM)", "FRASER 2.0\n(OMIM + RV)"))]

# fc_dt <- dcast(num_out_dt, dataset_label + dataset_group ~ method, value.var="num_out", fun.aggregate=function(x) as.numeric(median(x)))[,fc := FRASER / `FRASER 2.0`]
# fc_dt

g_boxplot_all <- ggplot(num_out_dt, aes(dataset_label, num_out+1, col=method)) +
    facet_grid(~dataset_group, scales="free_x", space="free_x") +
    # geom_violin() +
    # geom_boxplot(width=0.25, alpha=0.2, outlier.size=point_size) +
    geom_boxplot(outlier.size=point_size) +
    labs(x="", y="Splicing outliers\nper sample + 1") +
    scale_y_log10() +
    scale_color_manual(values=c("dodgerblue3", "purple4", "purple1", "violetred")) +
    theme_pubr() +
    theme(
        legend.title = element_blank(),
        legend.position = "top",
        # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        # axis.text.x=element_text(angle=90, hjust=1, vjust=1),
        # plot.margin=unit(c(0.5,0.5,1,0.5), "cm"),
        axis.title=element_text(face="bold"),
        text=element_text(size=font_size)) +
    cowplot::background_grid(major="y", minor="y")
a <- annotation_logticks(sides='l')
a$data <- data.frame(x=NA, dataset_group="Mitochondrial disease cohort") # have logticks only on left facet
g_boxplot_all <- g_boxplot_all + a
g_boxplot_all

# fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
upset_dt <- fdr_comparison_plots[["gg_upset"]]$data
upset_dt[, tmp := sapply(comb, FUN=function(x) paste(unlist(x), collapse="-"))]
set_order <- c(
    "FRASER-FRASER2-FRASER2 (OMIM + RV)-pathogenic splice defect",
    "FRASER-FRASER2-pathogenic splice defect",                    
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
    labs(x="", y="") +
    geom_text(aes(label = ..count..), stat = 'count', nudge_y = 500, 
              size=font_size/.pt) +
    theme(text=element_text(size=font_size)) +
    theme_combmatrix(combmatrix.label.text = element_text(size=font_size)) + 
    cowplot::background_grid(major="y", minor="y")
upset_plot

#+ power analysis results
# pathogenic cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj) & !grepl("deletion|cnv", VARIANT_EFFECT),]
# patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] 
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "C19orf70"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]

# transcriptome-wide results
res_all <- fread(snakemake@input$power_analysis_full_res_all)
res_tps <- fread(snakemake@input$power_analysis_full_res_patho)
# dont show sizes 70 and 90 in final plot
res_all <- res_all[!size %in% c(70, 90)]
res_tps <- res_tps[!size %in% c(70, 90)]

# plot proportion of recovered pathogenic cases
gprop <- ggplot(res_tps[padjustGene <= 0.1, .(Noutliers=.N, prop=.N/length(patho_tmp)), by=.(size, sim, FDR_set)], 
                aes(as.factor(size), 100*prop)) +
    geom_beeswarm(groupOnX=TRUE, size=0.75) +  
    labs(
        y = paste0('Percentage of recovered\n splicing outliers (n=', length(patho_tmp), ')'),
        x = 'Sample size') +
    scale_y_continuous(limits=c(0, 100)) +
    theme_pubr() +
    theme(axis.title=element_text(face="bold"),
          text=element_text(size=font_size),
          axis.text.x=element_text(size=font_size,
                                   angle=45, hjust=1, vjust=0.9)
        ) + 
    cowplot::background_grid(major="xy", minor="xy")
# gprop

#+ combine panels into figure, width=15, height=12
# gg_row1 <- ggarrange(
#     boxplot_prokisch,
#     upset_plot,
#     nrow=1, ncol=2,
#     labels=LETTERS[1:2],
#     widths=c(1.1,2)
#     )
# gg_row2 <- ggarrange(
#     gprop,
#     # g_boxplot_udn,
#     nrow=1, ncol=2,
#     labels=LETTERS[3],
#     widths=c(2,1)
#     )
# gg_figure <- ggarrange(gg_row1,
#                        gg_row2,
#                        # gg_row3,
#                        nrow=2, ncol=1,
#                        heights=c(1,1))

gg_row2 <- ggarrange(
    upset_plot,
    gprop,
    nrow=1, ncol=2,
    labels=LETTERS[2:3],
    widths=c(1.5,1)
)
gg_figure <- ggarrange(g_boxplot_all,
                       gg_row2,
                       nrow=2, ncol=1,
                       heights=c(1,1), 
                       labels=c(LETTERS[1], ""))
gg_figure

#+ save figure as png and pdf
ggsave(plot=gg_figure, filename=snakemake@output$outPng, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_figure, filename=snakemake@output$outPdf, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300) 
