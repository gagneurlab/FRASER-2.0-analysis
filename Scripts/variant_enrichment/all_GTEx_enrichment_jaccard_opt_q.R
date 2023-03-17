#'---
#' title: Compare jaccard to old version, on opt q, for all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_var_enrich_optQ.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 25000
#'   input:
#'     - rareSplice: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/optQ/jaccard_enrichment_data.Rds", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/variant_enrichment/GTEx_tissues/optQ_enrichment_jaccard.html"`'
#'    - outPng: '`sm config["htmlOutputPath"] + "/variant_enrichment/GTEx_tissues/optQ_enrichment_jaccard.png"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
.libPaths("~/R/4.1/FRASER2")
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
outPng         <- snakemake@output$outPng
threads        <- snakemake@threads

# METHODS_2_PLOT <- c("FRASER_p", "FRASER_psi5_psi3_p", "FRASER_theta_p")
METHODS_2_PLOT <- c("FRASER_p", "FRASER2_delta>=0.1_p", "FRASER2_delta>=0.2_p", "FRASER2_delta>=0.3_p", "FRASER2_delta>=0.3+blacklist_p", "FRASER2_delta>=0.4_p", "FRASER2_delta>=0.5_p")
PVALUES_2_PLOT <- c(-5, -7, -9)
FDR_LIMIT      <- 0.1
BPPARAM        <- MulticoreParam(threads, 100, progre=TRUE)

length(rareSpliceLs)
rareSpliceLs[1:5]
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
dt2plotSplice[,snptype:="Splice region"]

dt2plot <- dt2plotSplice

# dt2plot[,snptype:=factor(snptype, levels=c("Splice region", "MMSplice"))]
dt2plot[,snptype:=factor(snptype, levels=c("Splice region"))]

#'
#' Create P value enrichment plot
#'
#' FRASER2 vs old FRASER (and others)
#'
frdt   <- dt2plot[cutoff == TRUE & Method == "FRASER2_delta>=0.3+blacklist_p",      .(
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
# dt[,snptype:=factor(gsub(" ", "~", snptype), levels=c("Splice~region", "MMSplice"))]
dt[,snptype:=factor(gsub(" ", "~", snptype), levels=c("Splice~region"))]
levels(dt$Method) <- gsub(" ", "~", levels(dt$Method))

g1 <- ggplot(dt[snptype == "Splice~region"], aes(enrich, FRASER2)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER2)") +
    xlab("Enrichment (other method)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g1

# g2 <- ggplot(dt[snptype == "MMSplice"], aes(enrich, FRASER2)) +
#     geom_abline(slope=1, intercept=0) +
#     geom_point(color="gray40", alpha=0.6) +
#     ylab("Enrichment (FRASER2)") +
#     xlab("Enrichment (other method)") +
#     cowplot::theme_cowplot() +
#     geom_point(aes(x=1,y=1), col="white", alpha=0) +
#     theme_cowplot() +
#     grids() +
#     facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
#     scale_x_log10() +
#     scale_y_log10()
# g2

# only FRASER2 against old FRASER version
g1_vs_fraser_only <- ggplot(dt[snptype == "Splice~region" & Method == "FRASER"], aes(enrich, FRASER2)) +
    geom_abline(slope=1, intercept=0) +
    geom_point(color="gray40", alpha=0.6) +
    ylab("Enrichment (FRASER2)") +
    xlab("Enrichment (FRASER)") +
    cowplot::theme_cowplot() +
    geom_point(aes(x=1,y=1), col="white", alpha=0) +
    theme_cowplot() +
    grids() +
    facet_grid(facets=score + snptype ~ Method, scales="free", labeller=label_parsed) +
    scale_x_log10() +
    scale_y_log10()
g1_vs_fraser_only

#'
#' Arrange the plots
#'
# g <- ggarrange(labels=letters[1:2], align="hv", ncol=1,
#                g1,
#                g2)
g <- g1
# g


#+ save figure
factor <- 0.55
outPng
ggsave(outPng, g, width = 16*factor, height = 14*factor)
# ggsave(outPdf, g, width = 16*factor, height = 10*factor)

