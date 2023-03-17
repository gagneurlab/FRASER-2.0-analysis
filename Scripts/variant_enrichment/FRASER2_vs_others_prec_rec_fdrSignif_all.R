#'---
#' title: Compare final FRASER2 AUPRC values across tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_var_enrich_auprc_all_tissues.Rds"`'
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fraser2_vs_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_rv_recall_data_rareSpliceAI.Rds", dataset=["Skin_-_Not_Sun_Exposed_Suprapubic"])`'
#'     - auprc: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_auprc_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/all_final_prec_rec.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/all_final_prec_rec_ggplots.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))
source("src/R/enrichment_helper.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size

#+ read in recall at 10/20 outliers data
eval_dt <- rbindlist(bplapply(snakemake@input$auprc, function(file){
    eval_dt <- readRDS(file)
    # eval_dt <- melt(eval_dt, id.vars="Method")
    t <- basename(dirname(dirname(dirname(file))))
    eval_dt[, tissue:=t]
    snptype <- grep("rare(.)+",
                    strsplit(basename(file), "_", fixed=T)[[1]],
                    value=T)
    snptype <- gsub(".Rds", "", snptype)
    eval_dt[, snptype:=snptype]
    return(eval_dt)
}))

eval_dt[snptype == "rareSpliceAI", snptype:="rare SpliceAI"]
eval_dt[snptype == "rareMMSplice", snptype:="rare MMSplice"]
eval_dt[snptype == "rareAbSplice", snptype:="rare AbSplice"]
eval_dt[snptype == "rareSplicing", snptype:="rare splice site vicinity (VEP)"]
eval_dt[, snptype:=factor(snptype, levels=c("rare AbSplice", "rare MMSplice",
                                              "rare SpliceAI",
                                              "rare splice site vicinity (VEP)"))]
eval_dt <- eval_dt[Method != "IntronJaccardIndex",]
method_labels <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER2")
eval_dt[, Method := factor(Method, levels=method_labels)]
# eval_dt[Method == "IntronJaccardIndex", Method := "FRASER\n+ Intron Jaccard Index"]
# eval_dt[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER\n+ Intron Jaccard Index", "FRASER2"))]
# method_labels <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER\n+ Intron Jaccard Index", "FRASER2")

#+ create list of ggplots to save
ggplots <- list()

#+ aucpr_across_tissues, fig.width=21, fig.height=14
g_auprc_across_tissues <- ggplot(eval_dt, aes(Method, AUPRC, fill=Method)) +
    # facet_grid(mean_outliers_per_sample ~ snptype, labeller=label_both, scales="free_y") +
    facet_grid(. ~ snptype) +
    geom_boxplot() +
    labs(x="", y="Area under the precision-recall curve") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    # scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_auprc_across_tissues

#+ average_precision_across_tissues, fig.width=21, fig.height=14
g_ap_across_tissues <- ggplot(eval_dt, aes(Method, average_precision, fill=Method)) +
    # facet_grid(mean_outliers_per_sample ~ snptype, labeller=label_both, scales="free_y") +
    facet_grid(. ~ snptype) +
    geom_boxplot() +
    labs(x="", y="Average precision") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    # scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_ap_across_tissues

#+ auprc_across_tissues_spliceAI, fig.width=21, fig.height=14
g_auprc_across_tissues_spliceAI <- ggplot(eval_dt[snptype == "rare SpliceAI",], 
                                                   aes(Method, AUPRC, fill=Method)) +
    geom_boxplot() +
    labs(x="", y="Area under the precision-recall curve") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    # scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_auprc_across_tissues_spliceAI

#+ auprc_across_tissues_VEP, fig.width=21, fig.height=14
g_auprc_across_tissues_VEP <- ggplot(eval_dt[snptype == "rare splice site vicinity (VEP)",], 
                                              aes(Method, AUPRC, fill=Method)) +
    geom_boxplot() +
    labs(x="", y="Area under the precision-recall curve") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    # scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_auprc_across_tissues_VEP

#+ average_precision_across_tissues_spliceAI, fig.width=21, fig.height=14
g_ap_across_tissues_spliceAI <- ggplot(eval_dt[snptype == "rare SpliceAI",], 
                                          aes(Method, average_precision, fill=Method)) +
    geom_boxplot() +
    labs(x="", y="Average precision") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    # scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_ap_across_tissues_spliceAI

ggplots[["auprc_across_tissues"]] <- g_auprc_across_tissues
ggplots[["ap_across_tissues"]] <- g_ap_across_tissues
ggplots[["auprc_across_tissues_spliceAI"]] <- g_auprc_across_tissues_spliceAI
ggplots[["ap_across_tissues_spliceAI"]] <- g_ap_across_tissues_spliceAI
ggplots[["auprc_across_tissues_VEP"]] <- g_auprc_across_tissues_VEP


# recall vs rank plots
rrdt <- readRDS(snakemake@input$fraser2_vs_all)$recallData
dt4cutoffs <- rrdt[!is.na(Cutoff) & Method != "IntronJaccardIndex"]
rrdt <- rrdt[is.na(Cutoff) & Method != "IntronJaccardIndex"]
method_labels <- levels(eval_dt[, Method])

#'
#' The recall plots
#'
#+ recall plots 1, fig.height=6, fig.width=10
snptype <- "rareSpliceAI"
dataset <- "Skin not sun-exposed"
for(maxRank in c(10000, 20000, 50000, 100000, 200000)){
    
    ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
        rrdt, maxRank=maxRank, maxPoints=1e4) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall of rare \nsplice affecting variants") +
        grids(color="white") +
        xlim(0, maxRank) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        geom_abline(intercept=0, slope=rrdt[Method == "totalPossibleRank",1/rank], col="firebrick", linetype="dashed") + 
        # scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
        #                    labels=method_labels) +
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                           labels=method_labels) +
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Method", order = 1)) +
        theme_bw() + theme(text=element_text(size=16)) +
        xlim(0, maxRank) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) 
    ggplots[[paste0('recall_n=', maxRank)]]
    
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]] <- plotRecallRankForEnrichment(
        rrdt, maxRank=maxRank, maxPoints=1e4) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall of rare \nsplice affecting variants") +
        grids(color="white") +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        # scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
        #                    labels=method_labels) +
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                           labels=method_labels) +
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Splice metric", order = 1)) +
        theme_bw() + theme(text=element_text(size=16)) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) +
        scale_x_log10(limits=c(1, maxRank)) +
        annotation_logticks(sides="b")
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]]
}

ggplots[["precision_recall"]] <- ggplot(rrdt[Method != "totalPossibleRank"], aes(x=recall, y=precision, col=Method)) +
    geom_line() +
    geom_point(data=dt4cutoffs, aes(x=recall, y=precision, color=Method, shape=Cutoff), size=3) +
    labs(title=paste(dataset), 
         x="Recall of rare splice affecting variants", y="Precision") +
    grids(color="white") +
    # scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
    #                    labels=method_labels) +
    scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                       labels=method_labels) +
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
           color=guide_legend(title="Method", order = 1)) +
    theme_bw() + 
    theme(text=element_text(size=12)) 
ggplots[["precision_recall"]]

ggplots[["precision_recall"]] + 
    geom_line(size=0.7) +
    labs(title="") + 
    theme_pubr() + 
    theme(legend.position=c(0.8, 0.7), 
          legend.title=element_text(size=font_size),
          legend.text=element_text(size=font_size-2),
          axis.title=element_text(face="bold"))

ggplots[[paste0('recall_n=', maxRank)]] + 
    geom_line(size=0.7) +
    labs(title="") + 
    theme_pubr() + 
    theme(legend.position=c(0.7, 0.3), 
          legend.title=element_text(size=font_size),
          legend.text=element_text(size=font_size-2),
          axis.title=element_text(face="bold"))

ggplots[[paste0('recall_n=', maxRank, "_xlog")]] + 
    geom_line(size=0.7) +
    labs(title="") + 
    theme_pubr() + 
    theme(legend.position=c(0.2, 0.7), 
          legend.title=element_text(size=font_size),
          legend.text=element_text(size=font_size-2),
          axis.title=element_text(face="bold"))

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
