#'---
#' title: Paper figure 1
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/fig1.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - variant_enrich_comparison: '`sm expand(config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_vs_jaccard_rv_recall_plots_rare{snptype}.Rds", snptype=["Splicing", "MMSplice", "SpliceAI", "AbSplice"])`'
#'     - variant_enrich_comparison_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_types_vs_jaccard_rv_recall_data_rare{snptype}.Rds", snptype=["Splicing", "MMSplice", "SpliceAI", "AbSplice"])`'
#'   output:
#'    - outSvg: '`sm config["PAPER_FIGDIR"] + "/Fig1d.svg"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/Fig1d.tiff"`'
#'   type: script
#'---

# #'     - f1_res_table: '/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/Skin_-_Not_Sun_Exposed_Suprapubic/results_per_junction.tsv'
# #'     - f1_fds: '/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/datasets/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic--gencode34/fds-object.RDS'
# #'     - f2_fds : '`sm config["DATADIR"] + "/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic__optQ__newFilt/fds-object.RDS"`'

# # bash command for creating pysashimi plot
# cd Projects/external_tools/pysashimi/
# conda activate pysashimi_is
# python main.py -g ../../Test/gencode.v34.annotation.gtf -b /s/project/fraser/fraser2/figures/paper_figures/pysashimi/bam_sa_sashimi.tsv -e 12:52677250-52677800:- \
# --config ./settings.ini -o /s/project/fraser/fraser2/figures/paper_figures/pysashimi/ch12_52677250_5267800_psi5_v3.png --remove-empty-gene --threshold 20 --color-factor 3 --indicator-lines 52677645
# add anno in slides
# then: convert ch12_52677250_5267800_psi5_v2_annotated.png -trim ch12_52677250_5267800_psi5_v2_annotated_trim.png
# or:
# python main.py -g ../../Test/gencode.v34.annotation.gtf -b /s/project/fraser/fraser2/figures/paper_figures/pysashimi/bam_sa_sashimi.tsv -e 12:52677250-52677800:- \
# --config ./settings.ini -o /s/project/fraser/fraser2/figures/paper_figures/pysashimi/ch12_52677250_5267800_psi5.png --transcripts-to-show KRT1|KRT1-201 \
# --indicator-lines 52677645 --threshold 20 --color-factor 3

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(data.table)
# .libPaths("~/R/4.1/FRASER2")
# library(FRASER)
library(ggplot2)
library(ggpubr)
# library(gridExtra)
library(cowplot)

#+ read in figure font size and width params from config
font_size <- 12 # snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5


# #+ variant enrichment comparison to FRASER1
# var_enrich_comp_plot <- readRDS(snakemake@input$variant_enrich_comparison[1])
# maxRank <- 100000
# g_var_enrich <- var_enrich_comp_plot[[paste0('recall_n=', maxRank)]] +
#     labs(title="", y="Recall of rare splice-disrupting\ncandidate variants", x="Top N outliers")+ 
#     # geom_line(size=1.5) +
#     theme_pubr() + 
#     theme(legend.position="right", 
#           legend.title=element_text(size=font_size),
#           legend.text=element_text(size=font_size-2),
#           axis.title=element_text(face="bold")) #, size=14))

#+ var recall for all variant sets
var_recall_all_VEP <- readRDS(snakemake@input$variant_enrich_comparison[[1]])
var_recall_all_SpliceAI <- readRDS(snakemake@input$variant_enrich_comparison[[2]])
var_recall_all_MMSplice <- readRDS(snakemake@input$variant_enrich_comparison[[3]])
var_recall_all_AbSplice <- readRDS(snakemake@input$variant_enrich_comparison[[4]])
maxRank <- 1e5

getMetricLabels <- function(methods){
    vapply(methods, FUN = function(x) switch(x, 
                                             FRASER = c(bquote(psi[5]~","~psi[3]~","~theta)), 
                                             IntronJaccardIndex = c(bquote(Intron~Jaccard~Index))),
           FUN.VALUE = c(bquote(psi[3])))
}

recall_rank_dt <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                         (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare SpliceAI"],
                         (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare MMSplice"],
                         (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare AbSplice"])
dt4cutoffs <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare splice site vicinity\n(VEP)"],
                     (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare SpliceAI"],
                     (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare MMSplice"],
                     (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare AbSplice"])
recall_rank_dt[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
dt4cutoffs[, snptype := factor(snptype, levels=c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))]
# var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare AbSplice")
var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice")
methods_to_show <- c("FRASER", "IntronJaccardIndex")
recall_rank_dt <- recall_rank_dt[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
dt4cutoffs <- dt4cutoffs[snptype %in% var_sets_to_show & Method %in% methods_to_show,]
g_var_rank_rec  <- ggplot(recall_rank_dt, aes(rank, recall, col=Method)) +
    facet_wrap(~snptype) +
    geom_line() +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[3]]$data$slope, 
                col="firebrick", linetype="dashed") + 
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    scale_color_brewer(palette="Paired", 
                       labels=getMetricLabels(as.character(unique(recall_rank_dt[, Method])))) +
    # scale_x_continuous(breaks=seq(0, maxRank, by=1250000),
    #                    # labels=function(x) ifelse(x == 2.5e6 | x == 7.5e6, "", scientific_10(x)),
    #                    labels=c(0, scientific_10(1250000), "", scientific_10(3750000), ""),
    #                    limits=c(0, maxRank)) +
    labs(title="", 
         x="Top N outliers", y="Recall of rare splice-disrupting\ncandidate variants") + 
    guides(linetype = "none",
           color=guide_legend(order=1, nrow=2, title="Splice metric"),
           shape=guide_legend(order=2,nrow=2, title="Nominal\np-value\ncutoff")) + 
    theme_pubr() + 
    cowplot::background_grid(major="xy", minor="xy") +
    theme(
        legend.position="top",
        # legend.position=c(0.17, 0.95), # for rank-recall plot of all nominal pvals
        # legend.background=element_rect(fill='transparent'),
        # legend.position=c(0.85, 0.3), # for rank-recall plot of FDR signif
        legend.title=element_text(size=font_size),
        legend.text=element_text(size=font_size-2),
        text=element_text(size=font_size),
        axis.title=element_text(face="bold"),
        plot.margin=unit(c(0,1,0.25,0.25), units="cm")
    ) 

# #+ motivation example
# res_f1 <- fread(snakemake@input$f1_res_table)
# example <- res_f1[sampleID == "GTEX-1117F-2926-SM-5GZYI" & seqnames == "chr12" & hgncSymbol == "KRT1" & start == 52677481 & end == 52677645]
# fds_f1 <- loadFraserDataSet(file=snakemake@input$f1_fds)
# idx <- FRASER:::getIndexFromResultTable(fds_f1, example)
# psi5_vals <- assay(fds_f1, "psi5")[idx,example$sampleID]
# delta_psi5_vals <- deltaPsiValue(fds_f1, type="psi5")[idx,example$sampleID]
# psi3_vals <- assay(fds_f1, "psi3")[idx,example$sampleID]
# delta_psi3_vals <- deltaPsiValue(fds_f1, type="psi3")[idx,example$sampleID]
# # pred_psi_vals <- predictedMeans(fds_f1, type=example$type)[idx,example$sampleID]
# 
# fds_f2 <- loadFraserDataSet(file="/s/project/fraser/fraser2/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic__optQ__newFilt/fds-object.RDS")
# idx_f2 <- FRASER:::getIndexFromResultTable(fds_f2, example)
# jacc_vals <- assay(fds_f2, 'jaccard')[idx_f2,example$sampleID]
# delta_jacc_vals <- deltaPsiValue(fds_f2, type='jaccard')[idx_f2,example$sampleID]
# # pred_jacc_vals <- predictedMeans(fds_f2, type='jaccard')[idx_f2,example$sampleID]
# 
# dt <- data.table(metric=c("psi5", "psi3"), observed=c(psi5_vals, psi3_vals), delta=c(delta_psi5_vals, delta_psi3_vals))
# dt <- rbind(dt, data.table(metric="Intron~Jaccard~Index", observed=jacc_vals, delta=delta_jacc_vals))
# dt <- melt(dt, id.vars="metric")
# psi5_label <- FRASER:::ggplotLabelPsi("psi5", asCharacter = TRUE)[[1]]
# psi3_label <- FRASER:::ggplotLabelPsi("psi3", asCharacter = TRUE)[[1]]
# dt[metric == "psi5", metric := psi5_label]
# dt[metric == "psi3", metric := psi3_label]
# dt[, metric := factor(metric, levels=c(psi5_label, "Intron~Jaccard~Index", psi3_label))]
# setnames(dt, "metric", "Splice metric")
# 
# # g_bar_example <- ggbarplot(dt, x="variable", y="value", fill="Splice metric", 
# #                            position=position_dodge(0.7)) +
# #     scale_y_continuous(limits=c(0, 1)) + 
# #     scale_fill_brewer(palette="Paired", labels=c(expression(psi[5]), "Intron Jaccard Index")) +
# #     labs(y="", x="") 
# g_bar_example <- ggbarplot(dt, x="variable", y="value", fill="Splice metric", 
#                            position=position_dodge(0.7)) +
#     facet_wrap(~`Splice metric`, labeller=label_parsed) +
#     geom_hline(yintercept=0.3, color="black", linetype="dotted") +
#     scale_y_continuous(limits=c(0, 1)) + 
#     scale_fill_manual(values=c("firebrick3", "darkorchid4", "dodgerblue3")) +
#     labs(y="", x="") +
#     theme(legend.position="none")
# g_bar_example

#+ save figure as png and pdf
ggsave(plot=g_var_rank_rec, filename=snakemake@output$outSvg, width=1.33*page_width, height=page_width, unit=width_unit, dpi=350)
ggsave(plot=g_var_rank_rec, filename=snakemake@output$outTiff, width=1.33*page_width, height=page_width, unit=width_unit, dpi=350)
