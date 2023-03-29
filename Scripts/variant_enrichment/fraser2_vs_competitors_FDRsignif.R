#'---
#' title: Compare FRASER2 (jaccard) to FRASER1 and other competitors on FDR significant results only
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/{dataset}/FRASER2_vs_competitors_fdrSignif_rv_recall_{snptype}.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/" + config["datasets"][wildcards.dataset]["vcf_group"] + "/variant_extraction/" + wildcards.snptype + "_filtered_VariantsTable.tsv.gz"
#'   - |
#'    def get_spot_tissue_clean(wildcards):
#'     t = wildcards.dataset
#'     tissue = t.replace('-_', '')
#'     tissue = tissue.replace('-', '_')
#'     return config["spot_results"] + tissue + "/spot__fullResults.tsv"
#'    def get_leafcutterMD_tissue_clean(wildcards):
#'     t = wildcards.dataset
#'     tissue = t.replace('-_', '')
#'     tissue = tissue.replace('-', '_')
#'     return config["leafcutterMD_results"] + tissue + "/leafcutterMD_testing/results_" + tissue + ".tsv"
#'   threads: 5
#'   resources:
#'     - mem_mb: 40000
#'   input:
#'     - FRASER2_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_25_minN10/PCA__pc0.1/optQ__newFilt/jaccard/pval_aberrant_matrix__delta0.1.tsv.gz"`'
#'     - IntronJaccardIndex_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/optQ__newFilt/jaccard/pval_aberrant_matrix__delta0.3.tsv.gz"`'
#'     - FRASER_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/FRASER_pval_aberrant_matrix__delta0.3.tsv.gz"`'
#'     - SPOT_p: '`sm get_spot_tissue_clean`'
#'     - LeafcutterMD_p: '`sm get_leafcutterMD_tissue_clean`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/{dataset_group}/variant_enrichment/{dataset}/FRASER2_vs_competitors_fdrSignif_rv_recall_{snptype}.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_rv_recall_plots_{snptype}.Rds"`'
#'     - rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_rv_recall_data_{snptype}.Rds"`'
#'     - auprc_rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_auprc_{snptype}.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

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
dataset <- snakemake@wildcards$dataset
input_files <- snakemake@input[!names(snakemake@input) %in% c("", "variant_table", "RScript")]

# blacklist regions
blacklist_regions <- snakemake@config$blacklist_regions_hg38
blacklist_gr <- rtracklayer::import(blacklist_regions, format = "BED")
gene_info <- fread(snakemake@config$datasets[[dataset]]$orgdb)
gene_gr <- makeGRangesFromDataFrame(gene_info, keep.extra.columns=TRUE)
findOverlaps(blacklist_gr, gene_gr, type="any")
blacklist_genes <- as.data.table(gene_gr[to(findOverlaps(blacklist_gr, gene_gr, type="any")),])
blacklist_genes <- blacklist_genes[!duplicated(blacklist_genes)]
blacklist_genes[, isBlacklist := TRUE]
nrow(blacklist_genes)
 
#+ read in variants
variants <- fread(variants_file)

#+ read input files
res_for_comparision <- lapply(names(input_files), function(method_name){
# res_for_comparision <- lapply(names(input_files), function(method_name){
    pval_mat <- fread(input_files[[method_name]])
    
    if(grepl("SPOT", method_name)){
        spot_res <- pval_mat
        spot_res[gene_fdr > snakemake@config$fdrCutoff, gene_p := NA] # only consider FDR significant results
        spot_res <- spot_res[, .(SAMPLE_ID, GENE_ID, gene_p)]
        setnames(spot_res, "SAMPLE_ID", "sampleID")
        setnames(spot_res, "GENE_ID", "geneID")
        setnames(spot_res, "gene_p", method_name)
        setkey(spot_res, geneID, sampleID)
        spot_res <- spot_res[geneID != "",]
        return(spot_res)
    }
    
    if(grepl("LeafcutterMD", method_name)){
        lfMD_res <- pval_mat
        lfMD_res[padj > snakemake@config$fdrCutoff, pvalue_gene := NA] # only consider FDR significant results
        lfMD_res <- lfMD_res[, .(sample, geneID, pvalue_gene)]
        setnames(lfMD_res, "sample", "sampleID")
        setnames(lfMD_res, "pvalue_gene", method_name)
        setkey(lfMD_res, geneID, sampleID)
        lfMD_res <- lfMD_res[geneID != "",]
        return(lfMD_res)
    }
    
    # FRASER(2) results
    dt <- melt(pval_mat, id.vars="geneID", value.name=method_name, variable.name="sampleID")
    # # mask outliers in genes with many outliers
    # newColName <- paste0(method_name, "__nrGeneOutliers")
    # dt[, (newColName) := .SD[get(method_name) < 1e-7, .N], by="geneID"]
    # dt[get(newColName) > 3, (method_name) := 1]
    setkey(dt, geneID, sampleID)
    return(dt)
})

#+ merge to one res table
res_all <- Reduce(function(...) merge(..., all = TRUE), res_for_comparision)
setnames(res_all, "sampleID", "subjectID")
res_all[, subjectID:= gsub("^", "GTEX-", gsub("-.*", "", gsub("GTEX-", "", res_all$subjectID)))]
res_all[, tissue:=dataset]
res_all

# merge with blacklist regions
res_all <- merge(res_all, 
      blacklist_genes[, .(gene_name, isBlacklist)],
      all.x=TRUE,
      by.x="geneID", by.y="gene_name")
res_all[is.na(isBlacklist), isBlacklist := FALSE]
res_all[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(FRASER2_p))), by="isBlacklist"]
res_all[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(FRASER_p))), by="isBlacklist"]
res_all[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(SPOT_p))), by="isBlacklist"]
res_all[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(LeafcutterMD_p))), by="isBlacklist"]


#' Dataset: `r dataset`  
#' Number of samples: `r res_all[,uniqueN(subjectID)]`  
#' Number of genes: `r res_all[,uniqueN(geneID)]`  

#+ comparison methods
methodP <- grep('_n?p$', names(res_all), value=TRUE)
total_test <- sapply(methodP, function(x) sum(!is.na(res_all[,get(x)])))
sort(total_test)

#+ set unified factors across tables
usamples <- unique(c(res_all$subjectID))
ugenes   <- unique(c(res_all$geneID))
res_all[ ,subjectID :=factor(subjectID, levels=usamples)]
res_all[ ,geneID    :=factor(geneID,    levels=ugenes)]

variants[,subjectID :=factor(as.character(subjectID, levels=usamples))]
variants[,geneID    :=factor(as.character(SYMBOL,    levels=ugenes))]
variants[,IMPACT    :=factor(IMPACT, levels=c("HIGH", "MODERATE", "LOW"))]

# #+ replace NAs with 1s
# for(i in methodP){
#     res_all[is.na(get(i)), c(i):=list(1)]
# }
###########################################
#'
#' # Merge variants with outlier calls
#'
###########################################

#'
#' ## Get tested genes for all methods (FRASER, SPOT, LEAFCUTTER)
#' 
# overlapGenesTested <- !rowAnyNAs(as.matrix(res_all[,c(methodP), with=FALSE]))
# overlapGenesTested <- !rowAlls(as.matrix(res_all[,c(methodP), with=FALSE]), value=NA)
# res_all <- res_all[overlapGenesTested]
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

featuresSubset[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(FRASER2_p)), nrRareVars=sum(!is.na(simple_conseq))), by="isBlacklist"]
featuresSubset[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(FRASER_p)), nrRareVars=sum(!is.na(simple_conseq))), by="isBlacklist"]
featuresSubset[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(SPOT_p)), nrRareVars=sum(!is.na(simple_conseq))), by="isBlacklist"]
featuresSubset[, .(nrGenes=uniqueN(geneID), nrSamples=uniqueN(subjectID), nrOutliers=sum(!is.na(LeafcutterMD_p)), nrRareVars=sum(!is.na(simple_conseq))), by="isBlacklist"]

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
rrdt <- rbindlist(bplapply(methods2plot, dt=mafsub, BPPARAM=BPPARAM,
                           function(x, dt){
                               totalHits <- sum(!is.na(dt[,simple_conseq]))
                               dt <- calculateRecallRank(dt[!is.na(get(x)),], x, grepl('_n?p$', x))
                               dt <- data.table(Method=x, 
                                                dt[,.(nTrueHits=get(paste0(x, "_recall")),
                                                      recall=get(paste0(x, "_recall"))/totalHits,
                                                      precision=get(paste0(x, "_recall")) / get(paste0(x, '_rank')),
                                                      rank=get(paste0(x, '_rank')))])
                               dt[,Type:="FDR"]
                               dt[,Method:=gsub('_n?[pz]$', '', Method)]
                           }))
method_order <- c("LeafcutterMD", "SPOT", "FRASER", "IntronJaccardIndex", "FRASER2")
rrdt[, Method:=factor(Method, levels=method_order)]
rrdt <- rbind(rrdt, data.table(Method="totalPossibleRank", nTrueHits=sum(!is.na(mafsub[,simple_conseq])), recall=0, precision=0, rank=mafsub[, .N], Type="FDR"))
rrdt


#'
#' Merge with cutoffs from enrichment
#'
dt_tmp <- data.table(Method = method_order, 
                     Type = "FDR", Cutoff = snakemake@config$fdrCutoff)
dt4cutoffs <- merge(dt_tmp, rrdt[, .SD[rank == max(rank)], by="Method"])

dt4cutoffs[, Method:=factor(Method, levels=method_order)]
dt4cutoffs[, Cutoff:=factor(Cutoff)]
dt4cutoffs



#+ save data for plotting
out_rds <- snakemake@output$rds
saveRDS(file=out_rds, object=list(variants=variants, 
                                  featuresSubset=featuresSubset,
                                  enrich_final_ls=enrich_final_ls,
                                  recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))

# dont have final FRASER2 in plots for first figure
# rrdt <- rrdt[Method != "FRASER2"]
# dt4cutoffs <- dt4cutoffs[Method != "FRASER2"]

getMetricLabels <- function(methods){
    vapply(methods, FUN = function(x) switch(x, 
                                             FRASER = "FRASER", 
                                             FRASER2 = "FRASER2",
                                             IntronJaccardIndex = "FRASER + Intron Jaccard Index",
                                             SPOT = "SPOT",
                                             LeafcutterMD = "LeafcutterMD")
    , FUN.VALUE = "a" )
}


#'
#' The recall plots
#'
#+ recall plots 1, fig.height=6, fig.width=10
# snptype    <- "rareSplicingVariants"
snptype <- snakemake@wildcards$snptype
metric_labels <- getMetricLabels(methods=sort(gsub("_p$", "", methods2plot)))
# maxRank <- 20000
for(maxRank in c(10000, 20000, 50000, 100000)){
    
    ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
        rrdt, maxRank=maxRank, maxPoints=1e4) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        xlim(0, maxRank) +
        ylim(0, rrdt[rank < maxRank, max(recall)]) +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        # geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
        geom_abline(intercept=0, slope=rrdt[Method == "totalPossibleRank",1/rank], col="firebrick", linetype="dashed") + 
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
                           labels=metric_labels) +
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
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        # geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
        # geom_abline(intercept=0, slope=rrdt[,1/mafsub[, .N]], col="firebrick", linetype="dashed") + 
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
                           labels=metric_labels) +
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
         x="Recall", y="Precision") +
    grids(color="white") +
    scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "darkblue", "purple4"),
                       labels=metric_labels) +
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
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

#average precision: AP=Σ(Rₙ-Rₙ₋₁)Pₙ
rrdt[, AP_i := (recall - data.table::shift(.SD[,recall], 1, type="lag") ) * precision, by="Method"]
ap_dt <- rrdt[, sum(.SD[2:.N, AP_i]), by="Method"]
setnames(ap_dt, "V1", "average_precision")

# combine eval metrics and save
eval_dt <- merge(eval_dt, ap_dt)[order(-AUPRC)]
saveRDS(eval_dt, file=snakemake@output$auprc_rds)

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
# out_rds <- snakemake@output$rds
# saveRDS(file=out_rds, object=list(variants=variants, 
#                                   featuresSubset=featuresSubset,
#                                   enrich_final_ls=enrich_final_ls,
#                                   recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))

