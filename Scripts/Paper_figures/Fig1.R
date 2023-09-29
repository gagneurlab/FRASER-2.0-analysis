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
#'     - variant_enrich_comparison: '`sm expand(config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_vs_jaccard_rv_recall_plots_rare{snptype}.Rds", snptype=["SpliceSite", "MMSplice", "SpliceAI", "AbSplice"])`'
#'     - variant_enrich_comparison_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/Skin_-_Not_Sun_Exposed_Suprapubic/plot_rds/FRASER2_enrichment/FRASER_types_vs_jaccard_rv_recall_data_rare{snptype}.Rds", snptype=["SpliceSite", "MMSplice", "SpliceAI", "AbSplice"])`'
#'   output:
#'    - outSvg: '`sm config["PAPER_FIGDIR"] + "/Fig1d.svg"`'
#'    - outTiff: '`sm config["PAPER_FIGDIR"] + "/Fig1d.tiff"`'
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/Fig1d.png"`'
#'   type: script
#'---

# # bash command for creating pysashimi plot
# cd Projects/external_tools/pysashimi/
# conda activate pysashimi_is
# python main.py -g ../../Test/gencode.v34.annotation.gtf -b ~/paper_figures/pysashimi/bam_sa_sashimi.tsv -e 12:52677250-52677800:- \
# --config ./settings.ini -o ~/paper_figures/pysashimi/ch12_52677250_5267800_psi5_v3.png --remove-empty-gene --threshold 20 --color-factor 3 --indicator-lines 52677645
# add anno in slides
# then: convert ch12_52677250_5267800_psi5_v2_annotated.png -trim ch12_52677250_5267800_psi5_v2_annotated_trim.png
# or:
# python main.py -g ../../Test/gencode.v34.annotation.gtf -b /~/paper_figures/pysashimi/bam_sa_sashimi.tsv -e 12:52677250-52677800:- \
# --config ./settings.ini -o ~/paper_figures/pysashimi/ch12_52677250_5267800_psi5.png --transcripts-to-show KRT1|KRT1-201 \
# --indicator-lines 52677645 --threshold 20 --color-factor 3

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- 14 # snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5

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

recall_rank_dt <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare direct splice site\nvariant (VEP)"], # "rare splice site vicinity\n(VEP)"],
                         (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare SpliceAI"],
                         (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare MMSplice"],
                         (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$data)[,snptype := "rare AbSplice"])
dt4cutoffs <- rbind( (var_recall_all_VEP[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare direct splice site\nvariant (VEP)"], # "rare splice site vicinity\n(VEP)"],
                     (var_recall_all_SpliceAI[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare SpliceAI"],
                     (var_recall_all_MMSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare MMSplice"],
                     (var_recall_all_AbSplice[[paste0('recall_n=', maxRank)]]$layers[[2]]$data)[,snptype := "rare AbSplice"])
recall_rank_dt[, snptype := factor(snptype, levels=c("rare direct splice site\nvariant (VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))] # "rare splice site vicinity\n(VEP)", 
dt4cutoffs[, snptype := factor(snptype, levels=c("rare direct splice site\nvariant (VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice"))] # "rare splice site vicinity\n(VEP)", 
# var_sets_to_show <- c("rare splice site vicinity\n(VEP)", "rare AbSplice")
var_sets_to_show <- c("rare direct splice site\nvariant (VEP)", "rare MMSplice", "rare SpliceAI", "rare AbSplice")
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
    theme_manuscript(fig_font=font, fig_font_size=font_size) +
    cowplot::background_grid(major="xy", minor="xy") +
    theme(
        plot.margin=unit(c(0,1,0.25,0.25), units="cm") 
    ) 

#+ save figure as png and pdf
ggsave(plot=g_var_rank_rec, filename=snakemake@output$outSvg, width=1.33*page_width, height=page_width, unit=width_unit, dpi=350)
ggsave(plot=g_var_rank_rec, filename=snakemake@output$outTiff, width=1.33*page_width, height=page_width, unit=width_unit, dpi=350)
ggsave(plot=g_var_rank_rec, filename=snakemake@output$outPng, width=1.33*page_width, height=page_width, unit=width_unit, dpi=350)
