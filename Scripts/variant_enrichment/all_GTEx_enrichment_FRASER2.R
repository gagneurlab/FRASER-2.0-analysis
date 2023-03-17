#'---
#' title: Compare jaccard to old version, on opt q, for all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_var_enrich_optQ_newFilt.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 25000
#'   input:
#'     - rareSplice: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_final_enrichment_data_rareSplicing.Rds", dataset=config["tissues_for_reproducibility"])`'
#'     - rareMMSplice: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_final_enrichment_data_rareMMSplice.Rds", dataset=config["tissues_for_reproducibility"])`'
#'     - rareSpliceAI: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_final_enrichment_data_rareSpliceAI.Rds", dataset=config["tissues_for_reproducibility"])`'
#'     - detailed_enrichment_done: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/selected_GTEx_tissues/FRASER2_enrichment_investigation.done"`'
#'     - reproducibility_done_k20_q95_n1: '`sm config["DATADIR"] + "/GTEx_v8/reproducibility/minK20_95_minN1/optQ/PCA__pc0.1/tissue_reproducibility.Rds"`' 
#'     - reproducibility_done_k20_q25_n10: '`sm config["DATADIR"] + "/GTEx_v8/reproducibility/minK20_25_minN10/optQ/PCA__pc0.1/tissue_reproducibility.Rds"`' 
#'     - reproducibility_done_k10_q25_n10: '`sm config["DATADIR"] + "/GTEx_v8/reproducibility/minK10_25_minN10/optQ/PCA__pc0.1/tissue_reproducibility.Rds"`' 
#'   output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/GTEx_tissues/FRASER2_enrichment_full.html"`'
#'    - outPng: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/GTEx_tissues/FRASER2_enrichment_full.png"`'
#'   type: noindex
#' output:
#'   html_document
#'---

# #'     - bestQ_PCA: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/minK20_95_minN1/PCA__pc01/jaccard/0.3__5__0.1__FALSE/all_q_enrichment_plots_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareMMSplice", "rareSpliceAI"])`'
# #'     - bestQ_BB: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/minK20_95_minN1/PCA-BB-Decoder/jaccard/0.3__5__0.1/all_q_enrichment_plots_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareMMSplice", "rareSpliceAI"])`'
# , "rareIntronic"
# #'     - rareIntronic: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_enrichment_data_rareIntronic.Rds", dataset=config["tissues_for_reproducibility"])`'
saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
.libPaths("~/R/4.1/FRASER2_BB_loss")
# message("libPath is ", .libPaths())
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(FRASER)
library(BiocParallel)
library(cowplot)

#+ input
rareSpliceLs   <- snakemake@input$rareSplice
rareMMSpliceLs <- snakemake@input$rareMMSplice
rareSpliceAILs <- snakemake@input$rareSpliceAI
outPng         <- snakemake@output$outPng
threads        <- snakemake@threads

FRASER1_IMPLEMENTATION <- "FRASER1_PCA_k20_q95_n1_p"
FRASER2_IMPLEMENTATION <- "FRASER2_PCA_k20_q25_n10_pc01_blacklist_p"
# METHODS_2_PLOT <- c("FRASER_p", "FRASER_psi5_psi3_p", "FRASER_theta_p")
# METHODS_2_PLOT <- c("FRASER_p", "FRASER2_delta>=0.1_p", "FRASER2_delta>=0.2_p", "FRASER2_delta>=0.3_p", "FRASER2_delta>=0.3+blacklist_p", "FRASER2_delta>=0.4_p", "FRASER2_delta>=0.5_p")
METHODS_2_PLOT <- c(FRASER1_IMPLEMENTATION, "FRASER2_PCA_k20_q95_n1_pc01_blacklist_p", "FRASER2_PCA_k10_q25_n10_pc01_blacklist_p")
# METHODS_2_PLOT <- c(FRASER1_IMPLEMENTATION)
PVALUES_2_PLOT <- c(-5, -7, -9)
FDR_LIMIT      <- 0.1
BPPARAM        <- MulticoreParam(threads, 100, progre=TRUE)

length(rareSpliceLs)
rareSpliceLs[1:5]
rareMMSpliceLs[1:5]
rareSpliceAILs[1:5]
outPng


#'
#' FUNCTIONS
#'
readEnrichmentFiles <- function(file){
    x <- readRDS(file)[['enrich_final_ls']]
    rbindlist(lapply(names(x)[1:4], function(name){
        tissue <- strsplit(name, ": ")[[1]][1]
        score <- strsplit(name, ": ")[[1]][2]
        
        tmp_ans <- x[[name]]
        tmp_ans[,tissue:=tissue]
        tmp_ans[,score:=gsub("Pval", "P <", score)]
        tmp_ans[,score:=gsub("Zscore", "Z score", score)]
        tmp_ans
    }))
}


#' 
#' read enrichment data
#' and only extract the final enrichment scores
dt2plotSplice <- rbindlist(
    bplapply(rareSpliceLs, readEnrichmentFiles, BPPARAM=BPPARAM))
dt2plotMMSplice <- rbindlist(
    bplapply(rareMMSpliceLs, readEnrichmentFiles, BPPARAM=BPPARAM))
dt2plotSpliceAI <- rbindlist(
    bplapply(rareSpliceAILs, readEnrichmentFiles, BPPARAM=BPPARAM))

dt2plot <- rbind(
    dt2plotSplice[,snptype:="Splice region"],
    dt2plotMMSplice[,snptype:="MMSplice"],
    dt2plotSpliceAI[,snptype:="SpliceAI"])
dt2plot[,snptype:=factor(snptype, levels=c("Splice region", "MMSplice", "SpliceAI"))]


#'
#' Create P value enrichment plot
#'
#' FRASER2 vs old FRASER (and others)
#'
# frdt   <- dt2plot[cutoff == TRUE & Method == "FRASER2_delta>=0.3+blacklist_p",      .(
frdt   <- dt2plot[cutoff == TRUE & Method == FRASER2_IMPLEMENTATION,      .(
    FRASER2=enrichment, FRASER2_min=min.ci, FRASER2_max=max.ci, tissue, score, snptype)]
otherdt <- dt2plot[cutoff == TRUE & Method %in% METHODS_2_PLOT, .(
    Method, enrich=enrichment, enrich_min=min.ci, 
    enrich_max=max.ci, tissue, score, snptype)]
# otherdt[, Method:=mName4Plot(Method, removeTest=TRUE, AE_Name=AE_METHOD)]
# otherdt[, Method:=factor(Method, levels=mName4Plot(METHODS_2_PLOT, removeTest=TRUE, AE_Name=AE_METHOD))]
otherdt[, Method:=gsub('_p$', '', Method)]
otherdt[, Method:=factor(Method)]

dt <- merge(frdt, otherdt)
dt <- dt[grepl(paste0("(", paste(PVALUES_2_PLOT, collapse="|"), ")$"), score)]
dt[,score:=gsub(" 1e-", " 10^-", score)]
dt[,score:=gsub("P ", "italic(P)", score)]
dt[,snptype:=factor(gsub(" ", "~", snptype), levels=c("Splice~region", "MMSplice", "SpliceAI"))]
# dt[,snptype:=factor(gsub(" ", "~", snptype), levels=c("Splice~region"))]
levels(dt$Method) <- gsub(" ", "~", levels(dt$Method))

#+ enrichment_per_snptype , fig.height=6, fig.width=10
g <- ggplot(dt[Method == gsub('_p$', '', FRASER1_IMPLEMENTATION),], 
                aes(enrich, FRASER2)) +
    facet_grid(facets=score ~ snptype, scales="free", labeller=label_parsed) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER2)") +
    xlab("Enrichment (FRASER)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    scale_x_log10() +
    scale_y_log10()
g

#+ enrichment_per_snptype_all_methods , fig.height=8, fig.width=15
g1 <- ggplot(dt[snptype == "Splice~region"], aes(enrich, FRASER2)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab(paste0("Enrichment (", FRASER2_IMPLEMENTATION, ")")) +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g1

g2 <- ggplot(dt[snptype == "MMSplice"], aes(enrich, FRASER2)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab(paste0("Enrichment (", FRASER2_IMPLEMENTATION, ")")) +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g2

g3 <- ggplot(dt[snptype == "SpliceAI"], aes(enrich, FRASER2)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab(paste0("Enrichment (", FRASER2_IMPLEMENTATION, ")")) +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g3

# # only FRASER2 against old FRASER version
# g1_vs_fraser_only <- ggplot(dt[snptype == "Splice~region" & Method == "FRASER"], aes(enrich, FRASER2)) +
#     geom_abline(slope=1, intercept=0) +
#     geom_point(color="gray40", alpha=0.6) +
#     ylab("Enrichment (FRASER2)") +
#     xlab("Enrichment (FRASER)") +
#     cowplot::theme_cowplot() +
#     geom_point(aes(x=1,y=1), col="white", alpha=0) +
#     theme_cowplot() +
#     grids() +
#     facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
#     scale_x_log10() +
#     scale_y_log10()
# g1_vs_fraser_only

#'
#' Arrange the plots
#'
# g <- ggarrange(labels=letters[1:3], align="hv", ncol=1,
#                g1,
#                g2,
#                g3)
# g <- g1
# g


#+ save figure
outPng
ggsave(outPng, g, width=12, height=6)
# ggsave(outPdf, g, width = 16*factor, height = 10*factor)

