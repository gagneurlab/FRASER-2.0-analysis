#'---
#' title: Compare final FRASER2 variant enrichment across tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_var_enrich_all_tissues.Rds"`'
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fraser2_vs_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER_vs_jaccard_rv_recall_data_{varType}.Rds", dataset=["Skin_-_Not_Sun_Exposed_Suprapubic"], varType=["rareSplicing", "rareMMSplice", "rareSpliceAI"])`'
#'     - recall_at_x_outliers: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER_vs_jaccard_recall_at_x_outliers_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice"])`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/all_final_rv_recall.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/all_final_rv_recall_ggplots.Rds"`'
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
recall_dt <- rbindlist(bplapply(snakemake@input$recall_at_x_outliers, function(file){
    recall_dt <- readRDS(file)[["recallAt20Outliers"]]
    recall_dt[, mean_outliers_per_sample:=20]
    recall_dt10 <- readRDS(file)[["recallAt10Outliers"]]
    recall_dt10[, mean_outliers_per_sample:=10]
    recall_dt <- rbind(recall_dt, recall_dt10)
    t <- basename(dirname(dirname(dirname(file))))
    recall_dt[, tissue:=t]
    snptype <- grep("rare(.)+",
                    strsplit(basename(file), "_", fixed=T)[[1]],
                    value=T)
    snptype <- gsub(".Rds", "", snptype)
    recall_dt[, snptype:=snptype]
    return(recall_dt)
}))

recall_dt[snptype == "rareSpliceAI", snptype:="rare SpliceAI"]
recall_dt[snptype == "rareMMSplice", snptype:="rare MMSplice"]
recall_dt[snptype == "rareSplicing", snptype:="rare splice site vicinity (VEP)"]
recall_dt[, snptype:=factor(snptype, levels=c("rare MMSplice",
                                              "rare SpliceAI",
                                              "rare splice site vicinity (VEP)"))]
recall_dt <- recall_dt[Method != "IntronJaccardIndex",]
recall_dt[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER2"))]
metric_labels <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER2")
# recall_dt[Method == "IntronJaccardIndex", Method := "FRASER\n+ Intron Jaccard Index"]
# recall_dt[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER\n+ Intron Jaccard Index", "FRASER2"))]
# metric_labels <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER\n+ Intron Jaccard Index", "FRASER2")

#+ create list of ggplots to save
ggplots <- list()

#+ recall_across_tissues, fig.width=21, fig.height=14
g_F2_across_tissues <- ggplot(recall_dt, aes(Method, recall, fill=Method)) +
    # facet_grid(mean_outliers_per_sample ~ snptype, labeller=label_both, scales="free_y") +
    facet_grid(snptype ~ mean_outliers_per_sample, labeller=label_both) +
    geom_boxplot() +
    labs(x="", y="Recall at x outliers/sample") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_F2_across_tissues

#+ recall_across_tissues_10_outliers, fig.width=21, fig.height=14
g_F2_across_tissues_10_outliers <- ggplot(recall_dt[mean_outliers_per_sample == 10,], 
                                          aes(Method, recall, fill=Method)) +
    facet_grid(~snptype) +
    geom_boxplot() +
    labs(x="", y="Recall at 10 outliers/sample") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_F2_across_tissues_10_outliers

#+ recall_across_tissues_20_outliers, fig.width=21, fig.height=14
g_F2_across_tissues_20_outliers <- ggplot(recall_dt[mean_outliers_per_sample == 20,], 
                                          aes(Method, recall, fill=Method)) +
    facet_grid(~snptype) +
    geom_boxplot() +
    labs(x="", y="Recall at 20 outliers/sample") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_F2_across_tissues_20_outliers

#+ recall_across_tissues_20_outliers_spliceAI, fig.width=21, fig.height=14
g_F2_across_tissues_20_outliers_spliceAI <- ggplot(recall_dt[mean_outliers_per_sample == 20 &
                                                        snptype == "rare SpliceAI",], 
                                          aes(Method, recall, fill=Method)) +
    geom_boxplot() +
    labs(x="", y="Recall of rare \nsplice affecting variants \nat 20 outliers/sample") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_F2_across_tissues_20_outliers_spliceAI

#+ recall_across_tissues_20_outliers_VEP, fig.width=21, fig.height=14
g_F2_across_tissues_20_outliers_VEP <- ggplot(recall_dt[mean_outliers_per_sample == 20 &
                                                                 snptype == "rare splice site vicinity (VEP)",], 
                                                   aes(Method, recall, fill=Method)) +
    # geom_violin() +
    # geom_boxplot(width=0.1, color="black", alpha=0.2) +
    geom_boxplot() +
    labs(x="", y="Recall of rare variants in \nthe splice site vicinity (VEP)\nat 20 outliers/sample") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) + 
    theme_pubr() +
    cowplot::background_grid(major="y", minor="y") +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text = element_text(size = font_size)) 
g_F2_across_tissues_20_outliers_VEP

ggplots[["recall_across_tissues"]] <- g_F2_across_tissues
ggplots[["recall_across_tissues_10_outliers"]] <- g_F2_across_tissues_10_outliers
ggplots[["recall_across_tissues_20_outliers"]] <- g_F2_across_tissues_20_outliers
ggplots[["recall_across_tissues_20_outliers_spliceAI"]] <- g_F2_across_tissues_20_outliers_spliceAI
ggplots[["recall_across_tissues_20_outliers_VEP"]] <- g_F2_across_tissues_20_outliers_VEP


# recall vs rank plots
for(recallData_file in snakemake@input$fraser2_vs_all){
    
    rrdt <- readRDS(recallData_file)$recallData
    dt4cutoffs <- rrdt[!is.na(Cutoff) & Method != "IntronJaccardIndex"]
    rrdt <- rrdt[is.na(Cutoff) & Method != "IntronJaccardIndex"]
    metric_labels <- levels(recall_dt[, Method])
    
    # recall plots 1, fig.height=6, fig.width=10
    snptype <- gsub("FRASER_vs_jaccard_rv_recall_data_", "", basename(recallData_file))
    snptype <- gsub(".Rds", "", snptype)
    message(date(), "snptype: ", snptype)
    dataset <- "Skin not sun-exposed"
    for(maxRank in c(10000, 20000, 50000, 100000, 200000)){
        
        ggplots[[paste0('recall_n=', maxRank, "__", snptype)]] <- plotRecallRankForEnrichment(
            rrdt, maxRank=maxRank, maxPoints=1e4) +
            labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
                 x="Top N outliers", y="Recall of rare \nsplice affecting variants") +
            grids(color="white") +
            xlim(0, maxRank) +
            ylim(0, rrdt[rank < maxRank, max(recall)]) +
            geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
            geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
            # scale_color_brewer(palette="Paired", labels=metric_labels) + 
            scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                               labels=metric_labels) +
            scale_shape_discrete(labels=function(x)parse(text=x)) +
            guides(linetype = "none") + 
            guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
                   color=guide_legend(title="Method", order = 1)) +
            theme_bw() + theme(text=element_text(size=16)) +
            xlim(0, maxRank) +
            ylim(0, rrdt[rank < maxRank, max(recall)]) 
        ggplots[[paste0('recall_n=', maxRank, "__", snptype)]]
        
        ggplots[[paste0('recall_n=', maxRank, "_xlog", "__", snptype)]] <- plotRecallRankForEnrichment(
            rrdt, maxRank=maxRank, maxPoints=1e4) +
            labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
                 x="Top N outliers", y="Recall of rare \nsplice affecting variants") +
            grids(color="white") +
            geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
            geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
            # scale_color_brewer(palette="Paired", labels=metric_labels) + 
            scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                               labels=metric_labels) +
            scale_shape_discrete(labels=function(x)parse(text=x)) +
            guides(linetype = "none") + 
            guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
                   color=guide_legend(title="Method", order = 1)) +
            theme_bw() + theme(text=element_text(size=16)) +
            ylim(0, rrdt[rank < maxRank, max(recall)]) +
            scale_x_log10(limits=c(1, maxRank)) 
        ggplots[[paste0('recall_n=', maxRank, "_xlog", "__", snptype)]]
    }

    # maxRank <- 1e+05
    # ggplots[[paste0('recall_n=', maxRank)]] + 
    #     # geom_line(size=1.5) + 
    #     labs(title="") + 
    #     theme_pubr() + 
    #     theme(legend.position="right", 
    #           legend.title=element_text(size=font_size),
    #           legend.text=element_text(size=font_size-2),
    #           axis.title=element_text(face="bold"))
}

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
