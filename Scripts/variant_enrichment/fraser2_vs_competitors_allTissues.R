#'---
#' title: Compare FRASER2 to FRASER1, LeafcutterMD and SPOT across all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/FRASER2_vs_others_allTissues_rv_recall_{snptype}.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/GTEx_v8/variant_extraction/" + wildcards.snptype + "_filtered_VariantsTable.tsv.gz"
#'   threads: 5
#'   resources:
#'     - mem_mb: 600000
#'   input:
#'     - SPOT: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/SPOT_allTissues_pvals_for_rv.Rds"`'
#'     - LeafcutterMD: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/LeafcutterMD_allTissues_pvals_for_rv.Rds"`'
#'     - FRASER: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER_allTissues_pvals_for_rv.Rds"`'
#'     - FRASER2: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER2_allTissues_pvals_for_rv.Rds"`'
#'     - FRASER2_noBlacklist: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER2_noBlacklist_allTissues_pvals_for_rv.Rds"`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/FRASER2_vs_others_allTissues_rv_recall_{snptype}.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissues_rv_recall_plots_{snptype}.Rds"`'
#'     - rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissues_rv_recall_data_{snptype}.Rds"`'
#'     - auprc_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissuesd_auprc_{snptype}.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

# #'     - recall_at_x_outliers_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissuesd_recall_at_x_outliers_{snptype}.Rds"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
.libPaths("~/R/4.1/FRASER2")
# message("libPath is ", .libPaths())
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(FRASER)
library(BiocParallel)
library(GenomicRanges)
source("src/R/enrichment_helper.R")

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ input
variants_file <- snakemake@input$variant_table
dataset <- "all tissues"
input_files <- snakemake@input[!names(snakemake@input) %in% c("", "variant_table", "RScript")]

#+ read in variants
variants <- fread(variants_file)

#+ read input files
res_for_comparision <- lapply(input_files, readRDS)

#+ merge to one res table
res_all <- Reduce(function(...) merge(..., all = TRUE), res_for_comparision)
setnames(res_all, "sampleID", "RNA_ID")
res_all[, subjectID:= gsub("^", "GTEX-", gsub("-.*", "", gsub("GTEX-", "", res_all$RNA_ID)))]
# res_all[, tissue:=dataset]
res_all

#' Dataset: `r dataset`  
#' Number of samples: `r res_all[,uniqueN(subjectID)]`  
#' Number of genes: `r res_all[,uniqueN(geneID)]`  

#+ comparison methods
methodP <- paste0(names(input_files), "_p")
total_test <- sapply(methodP, function(x) sum(!is.na(res_all[,get(x)])))
sort(total_test)

#+ replace NAs with 1s
for(i in methodP){
    res_all[is.na(get(i)), c(i):=list(1)]
}

#+ set unified factors across tables
usamples <- unique(c(res_all$subjectID))
ugenes   <- unique(c(res_all$geneID))
res_all[ ,subjectID :=factor(subjectID, levels=usamples)]
res_all[ ,geneID    :=factor(geneID,    levels=ugenes)]

variants[,subjectID :=factor(as.character(subjectID, levels=usamples))]
variants[,geneID    :=factor(as.character(SYMBOL,    levels=ugenes))]
variants[,IMPACT    :=factor(IMPACT, levels=c("HIGH", "MODERATE", "LOW"))]

###########################################
#'
#' # Merge variants with outlier calls
#'
###########################################

#'
#' ## Get tested genes for all methods (FRASER, SPOT, LEAFCUTTER)
#' 
overlapGenesTested <- !rowAnyNAs(as.matrix(res_all[,c(methodP), with=FALSE]))
res_all <- res_all[overlapGenesTested]
length(unique(res_all[,subjectID]))
length(unique(res_all[,geneID]))

#'
#' # Merge data sets into one table
#+ merge data
var2merge <- unique(variants[, .(subjectID, geneID, simple_conseq, MAF, IMPACT)])
setkey(var2merge, subjectID, geneID, IMPACT)
var2merge <- var2merge[!duplicated(var2merge, by=c("subjectID", "geneID"))]

featuresSubset <- merge(
    res_all, var2merge, all.x=TRUE, by=c("subjectID", "geneID"))
featuresSubset[, back_simple_conseq:=simple_conseq]
featuresSubset

ggplots <- list()

#'
#' # Enrichments (all)
#'
enrich_table <- list(
    # list(name="Pval 1e-3", methods=methodP, isZscore=FALSE, cutOff=1e-3),
    list(name="Pval 1e-5", methods=methodP, isZscore=FALSE, cutOff=1e-5),
    list(name="Pval 1e-7", methods=methodP, isZscore=FALSE, cutOff=1e-7),
    list(name="Pval 1e-9", methods=methodP, isZscore=FALSE, cutOff=1e-9)
)

enrich_final_ls <- list()
for(et in enrich_table){
    enrichdt <- rbindlist(lapply(et$methods, function(x){
        dt1 <- calculateEnrichment(featuresSubset, cols=x, cutOff=et$cutOff,
                                   isZscore=et$isZscore)
        ans <- dt1$dt[, .(cutoff, nRareEvent, total, fraction, nNA, Method=x,
                          enrichment=dt1$enrichment, 
                          min.ci=dt1$min.ci, max.ci=dt1$max.ci,
                          multiCall=1)]
        ans
    }))
    enrich_final_ls[[paste0(dataset, ": ", et$name)]] <- enrichdt
    
    ggplots[[paste0('enrichAll_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
        geom_point() +
        facet_grid(rows="multiCall", scales="free_y") + 
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
        coord_flip() +
        labs(title=paste0('Enrichment (All, \n', dataset, ', ', et$name, ')'))
}

#+ enrichment_all, fig.height=8, fig.width=16
names2plot <- paste0('enrichAll_', sapply(enrich_table, "[[", "name"))
etl <- length(enrich_table)/2
order2plot <- rep(1:etl, each=2) + rep(c(0, etl), etl)
grid.arrange(ncol=2, grobs=ggplots[names2plot][order2plot])
# ggarrange(plotlist = ggplots[names2plot], align = "hv")

#'
#' # Enrichment (rare: MAF < 0.001)
MAF_LIMIT <- as.numeric(snakemake@config$MAF_LIMIT)
if(snakemake@wildcards$snptype %in% c("rareMMSplice", "rareSpliceAI")){
    MAF_LIMIT <- 0.001
}
mafsub <- copy(featuresSubset)
mafsub[MAF > MAF_LIMIT, simple_conseq:=NA]
for(et in enrich_table){
    enrichdt <- rbindlist(lapply(et$methods, function(x){
        dt <- calculateEnrichment(mafsub, cols=x, cutOff=et$cutOff,
                                  isZscore=et$isZscore)
        dt$dt[,.(cutoff, nRareEvent, total, fraction, nNA, Method=x,
                 enrichment=dt$enrichment, min.ci=dt$min.ci, max.ci=dt$max.ci)]
    }))
    enrich_final_ls[[paste0(dataset, ": ", et$name, " + MAF < ", MAF_LIMIT)]] <- enrichdt
    
    ggplots[[paste0('enrichMAF_', MAF_LIMIT, '_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
        geom_point() +
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
        coord_flip() +
        labs(title=paste0('Enrichment (MAF < ', MAF_LIMIT, ', \n', dataset, ', ', et$name, ')'))
}

#+ enrichment_maf, fig.height=10, fig.width=20
names2plot <- paste0('enrichMAF_', MAF_LIMIT, '_', sapply(enrich_table, "[[", "name"))
etl <- length(enrich_table)/2
order2plot <- rep(1:etl, each=2) + rep(c(0, etl), etl)
grid.arrange(ncol=2, grobs=ggplots[names2plot][order2plot])

#'
#' # Recall Rank plots all (HIGH/MODERATE)
#'
#+ create recall data 1
#' add random
methods2plot <- methodP
rrdt <- rbindlist(lapply(methods2plot, dt=mafsub, # BPPARAM=BPPARAM,
                           function(x, dt){
                               totalHits <- sum(!is.na(dt[,simple_conseq]))
                               dt <- calculateRecallRank(dt, x, grepl('_n?p$', x))
                               dt <- data.table(Method=x, 
                                                dt[,.(nTrueHits=get(paste0(x, "_recall")),
                                                      recall=get(paste0(x, "_recall"))/totalHits,
                                                      precision=get(paste0(x, "_recall")) / get(paste0(x, '_rank')),
                                                      rank=as.double(get(paste0(x, '_rank'))))])
                               dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 
                                                ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score'))]
                               dt[,Method:=gsub('_n?[pz]$', '', Method)]
                               
                               dt <- dt[rank <= 1e7 | rank == max(rank),]
                               return(dt)
                           }))
method_order <- c("LeafcutterMD", "SPOT", "FRASER", "FRASER2", "FRASER2_noBlacklist")
# method_order <- c("FRASER", "FRASER2")
rrdt[, Method:=factor(Method, levels=method_order)]
rrdt


#'
#' Merge with cutoffs from enrichment
#'
#+ get pval recall at pval cutoffs
dt_tmp <- rbindlist(lapply(names(enrich_final_ls[4:6]), function(x){
    Cutoff <- strsplit(x, " ")[[1]][4]
    enrich_final_ls[[x]][cutoff == TRUE, .(
        rank=total,
        Method=gsub("_n?[pz]$", "", Method),
        Type=ifelse(grepl('_p$', Method), 'P-value', 
                    ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score')),
        Cutoff=Cutoff)]
}))
dt4cutoffs <- merge(dt_tmp, rrdt)[order(Method, Type)]
dt4cutoffs[,Cutoff:=gsub("1e-", "10^-", Cutoff)]
dt4cutoffs[, Method:=factor(Method, levels=method_order)]
dt4cutoffs

#'
#' # Recall at 20 outliers/sample
#'
#+ create recall at 20 outlier/sample data
# n_samples <- length(unique(res_all[,subjectID]))
# rank_cutoff_20 <- n_samples*20
# recall_20_dt <- rrdt[rank == rank_cutoff_20, ]
# ggplots[['recall_at_20_outliers']] <- ggplot(recall_20_dt, aes(Method, recall)) +
#     geom_bar(stat="identity") + 
#     labs(x="pseudocount")
# ggplots[['recall_at_20_outliers']]
# rank_cutoff_10 <- n_samples*10
# recall_10_dt <- rrdt[rank == rank_cutoff_10, ]
# ggplots[['recall_at_10_outliers']] <- ggplot(recall_10_dt, aes(Method, recall)) +
#     geom_bar(stat="identity") + 
#     labs(x="pseudocount")
# ggplots[['recall_at_10_outliers']]

#+ save data for plotting
out_rds <- snakemake@output$rds
saveRDS(file=out_rds, object=list(variants=variants, 
                                  featuresSubset=featuresSubset,
                                  enrich_final_ls=enrich_final_ls,
                                  recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))
# recall_rds <- snakemake@output$recall_at_x_outliers_rds
# saveRDS(file=recall_rds, object=list(recallAt10Outliers=recall_10_dt,
#                                      recallAt20Outliers=recall_20_dt))

#'
#' The recall plots
#'
#+ recall plots 1, fig.height=6, fig.width=10
# snptype    <- "rareSplicingVariants"
snptype <- snakemake@wildcards$snptype
metric_labels <- method_order
# maxRank <- 20000
for(maxRank in c(1e5, 2e5, 5e5, 1e6, 5e6, 1e7)){
    
    ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
        rrdt, maxRank=maxRank, maxPoints=1e5) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        xlim(0, maxRank) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4", "violetred"),
                           labels=metric_labels) +
        # scale_color_manual(values=c("dodgerblue3", "purple4"),
        #                    labels=metric_labels) +
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Method", order = 1)) +
        theme_bw() + theme(text=element_text(size=12)) +
        xlim(0, maxRank) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) 
    ggplots[[paste0('recall_n=', maxRank)]]
    
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]] <- plotRecallRankForEnrichment(
        rrdt, maxRank=maxRank, maxPoints=1e5) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4", "violetred"),
                           labels=metric_labels) +
        # scale_color_manual(values=c("dodgerblue3", "purple4"),
        #                    labels=metric_labels) +
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Method", order = 1)) +
        theme_bw() + theme(text=element_text(size=12)) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) +
        scale_x_log10(limits=c(1, maxRank)) 
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]]
}

maxPlottedRank <- rrdt[, .SD[precision > 0.01, max(rank)], by="Method"][, max(V1)]
pr_dt <- rbind(rrdt[rank <= maxPlottedRank], 
               rrdt[rank > maxPlottedRank,][sample(rrdt[rank > maxPlottedRank,.N], 5e4*length(method_order))])
ggplots[["precision_recall"]] <- ggplot(pr_dt, aes(x=recall, y=precision, col=Method)) +
    geom_line() +
    geom_point(data=dt4cutoffs, aes(x=recall, y=precision, color=Method, shape=Cutoff), size=3) +
    labs(title=paste(dataset), 
         x="Recall", y="Precision") +
    grids(color="white") +
    scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4", "violetred"),
                       labels=metric_labels) +
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
           color=guide_legend(title="Method", order = 1)) +
    theme_bw() + 
    theme(text=element_text(size=12)) 
ggplots[["precision_recall"]]

# calc AUPRC
require(PRROC)
eval_dt <- rbindlist(lapply(methods2plot, function(x, dt){
    # get data for aucpr
    fg <- dt[!is.na(simple_conseq), -get(x)]
    fg[is.na(fg)] <- -1
    bg <- dt[is.na(simple_conseq), -get(x)]
    bg[is.na(bg)] <- -1
    
    # PR Curve
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, sorted=F, curve = F)
    # plot(pr)
    return(data.table(Method=gsub('_n?[pz]$', '', x),
                      AUPRC=pr$auc.integral))
    
}, dt=mafsub))
eval_dt
saveRDS(eval_dt, file=snakemake@output$auprc_rds)


#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
# out_rds <- snakemake@output$rds
# saveRDS(file=out_rds, object=list(variants=variants, 
#                                   featuresSubset=featuresSubset,
#                                   enrich_final_ls=enrich_final_ls,
#                                   recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))

