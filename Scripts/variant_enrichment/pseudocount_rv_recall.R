#'---
#' title: Compare different pseudocount values (FRASER2)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/{dataset}/FRASER2_pseudocount_rv_recall_rho{rho}_{snptype}.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/" + config["datasets"][wildcards.dataset]["vcf_group"] + "/variant_extraction/" + wildcards.snptype + "_filtered_VariantsTable.tsv.gz"
#'   threads: 5
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - pc_1_0_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_5_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.5/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_25_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.25/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_1_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.1/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_05_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.05/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_01_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.01/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_001_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.001/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - pc_0_00001_p: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/minK20_95_minN1/PCA__pc0.00001/optQ__newFilt/jaccard/pval_matrix__0.3__5__{rho}__FALSE.tsv.gz"`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/{dataset_group}/variant_enrichment/{dataset}/FRASER2_pseudocount_rv_recall_rho{rho}_{snptype}.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_rv_recall_plots_rho{rho}_{snptype}.Rds"`'
#'     - rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_rv_recall_data_rho{rho}_{snptype}.Rds"`'
#'     - recall_at_x_outliers_rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_recall_at_x_outliers_rho{rho}_{snptype}.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---


saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
# .libPaths("~/R/4.1/FRASER2")
.libPaths("~/R/4.1/FRASER2_BB_loss")
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

# # blacklist regions
# blacklist_regions <- snakemake@config$blacklist_regions
# blacklist_gr <- rtracklayer::import(blacklist_regions, format = "BED")
# 
#+ read in variants
variants <- fread(variants_file)
# if(snakemake@wildcards$snptype %in% c("rareSpliceAI", "rareMMSplice")){
#     sep <- ":"
# } else {
#     sep <- "_"
# }
# variants[, seqnames := sapply(variantID, function(x, sep="_"){
#     strsplit(x, sep, fixed=TRUE)[[1]][1]
# }, sep=sep)]
# variants[, pos := sapply(variantID, function(x, sep="_"){
#     strsplit(x, sep, fixed=TRUE)[[1]][2]
# }, sep=sep)] 
# 
# # remove vars in blacklist regions
# var_gr <- makeGRangesFromDataFrame(df=variants, keep.extra.columns=TRUE, 
#                                    ignore.strand=TRUE, start.field="pos", end.field="pos")
# hits <- findOverlaps(var_gr, blacklist_gr, minoverlap=1, type="any" )
# variants <- variants[-from(hits),] 

#+ read input files
res_for_comparision <- lapply(names(input_files), function(method_name){
    pval_mat <- fread(input_files[[method_name]])
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

#' Dataset: `r dataset`  
#' Number of samples: `r res_all[,uniqueN(subjectID)]`  
#' Number of genes: `r res_all[,uniqueN(geneID)]`  

#+ comparison methods
methodP <- grep('_n?p$', names(res_all), value=TRUE)
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
#' # Enrichment (rare: MAF < 0.01)
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
#' # Recall Rank plots all 
#'
#+ create recall data 1
#' add random
methods2plot <- methodP
rrdt <- rbindlist(bplapply(methods2plot, dt=mafsub, BPPARAM=BPPARAM,
                           function(x, dt){
                               totalHits <- sum(!is.na(dt[,simple_conseq]))
                               dt <- calculateRecallRank(dt, x, grepl('_n?p$', x))
                               dt <- data.table(Method=x, 
                                                dt[,.(nTrueHits=get(paste0(x, "_recall")),
                                                      recall=get(paste0(x, "_recall"))/totalHits,
                                                      rank=get(paste0(x, '_rank')))])
                               dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 
                                                ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score'))]
                               dt[,Method:=gsub('_n?[pz]$', '', Method)]
                           }))
#'
#' Merge with cutoffs from enrichment
#'
dt_tmp <- rbindlist(lapply(names(enrich_final_ls[4:6]), function(x){
    Cutoff <- strsplit(x, " ")[[1]][3]
    enrich_final_ls[[x]][cutoff == TRUE, .(
        rank=total,
        Method=gsub("_n?[pz]$", "", Method),
        Type=ifelse(grepl('_p$', Method), 'P-value', 
                    ifelse(grepl('_np$', Method), 'P-value norm', 'Z-score')),
        Cutoff=Cutoff)]
}))
dt4cutoffs <- merge(dt_tmp, rrdt)[order(Method, Type)]
dt4cutoffs[,Cutoff:=gsub("1e-", "10^-", Cutoff)]
all_methods <- sort(gsub("_", ".", gsub("pc_", "", gsub("_p$", "", methods2plot))))
rrdt[, Method:=gsub("_", ".", gsub("pc_", "", Method))]
rrdt[, Method:=factor(Method,  levels=all_methods)]
rrdt
dt4cutoffs[, Method:=gsub("_", ".", gsub("pc_", "", Method))]
dt4cutoffs[, Method:=factor(Method,  levels=all_methods)]
dt4cutoffs

#+ subset data for plotting
pc_vals_to_plot <- c("0.01", "0.1", "0.5", "1.0")
rrdt_plot <- rrdt[Method %in% pc_vals_to_plot,]
dt4cutoffs_plot <- dt4cutoffs[Method %in% pc_vals_to_plot,]
colors_all <- RColorBrewer::brewer.pal(length(methods2plot), "Blues")
colors_sub <- colors_all[all_methods %in% pc_vals_to_plot]

#'
#' The recall plots
#'
#+ recall plots 1, fig.height=6, fig.width=10
# snptype    <- "rareSplicingVariants"
snptype <- snakemake@wildcards$snptype
maxRank <- 20000
for(maxRank in c(10000, 20000, 25000, 30000, 50000, 100000)){
    ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
        rrdt_plot, maxRank=maxRank, maxPoints=1e4) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        xlim(0, maxRank) +
        ylim(0, rrdt_plot[rank < maxRank, max(recall)]) +
        geom_point(data=dt4cutoffs_plot, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        geom_abline(intercept=0, slope=rrdt_plot[,1/max(rank)], col="firebrick", linetype="dashed") + 
        # scale_color_brewer(palette="Blues") + 
        scale_color_manual(values=colors_sub) + 
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Pseudocount", order = 1)) +
        theme_bw() + #theme(text=element_text(size=16)) +
        xlim(0, maxRank) +
        ylim(0, rrdt_plot[rank < maxRank, max(recall)]) 
    ggplots[[paste0('recall_n=', maxRank)]]
    
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]] <- plotRecallRankForEnrichment(
        rrdt_plot, maxRank=maxRank, maxPoints=1e4) +
        labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
             x="Top N outliers", y="Recall") +
        grids(color="white") +
        geom_point(data=dt4cutoffs_plot, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
        geom_abline(intercept=0, slope=rrdt_plot[,1/max(rank)], col="firebrick", linetype="dashed") + 
        # scale_color_brewer(palette="Blues") + 
        scale_color_manual(values=colors_sub) + 
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Pseudocount", order = 1)) +
        theme_bw() + theme(text=element_text(size=16)) +
        ylim(0, rrdt_plot[rank < maxRank, max(recall)]) +
        scale_x_log10(limits=c(1, maxRank)) 
    ggplots[[paste0('recall_n=', maxRank, "_xlog")]]
}

#'
#' # Recall at 20 outliers/sample
#'
#+ create recall at 20 outlier/sample data
n_samples <- length(unique(res_all[,subjectID]))
rank_cutoff_20 <- n_samples*20
recall_20_dt <- rrdt[rank == rank_cutoff_20, ]
recall_20_dt[, pc:=as.numeric(as.character(Method))]
recall_20_dt[, rho:=as.numeric(snakemake@wildcards$rho)]
ggplots[['recall_at_20_outliers']] <- ggplot(recall_20_dt, aes(factor(pc), recall)) +
    geom_bar(stat="identity") + 
    labs(x="pseudocount")
ggplots[['recall_at_20_outliers']]
ggplots[['recall_at_20_outliers_xlog']] <- ggplot(recall_20_dt, aes(pc, recall)) +
    geom_bar(stat="identity") + 
    scale_x_log10() +
    labs(x="pseudocount")
ggplots[['recall_at_20_outliers_xlog']]
rank_cutoff_10 <- n_samples*10
recall_10_dt <- rrdt[rank == rank_cutoff_10, ]
recall_10_dt[, pc:=as.numeric(as.character(Method))]
recall_10_dt[, rho:=as.numeric(snakemake@wildcards$rho)]
ggplots[['recall_at_10_outliers']] <- ggplot(recall_10_dt, aes(factor(pc), recall)) +
    geom_bar(stat="identity") + 
    labs(x="pseudocount")
ggplots[['recall_at_10_outliers']]

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
out_rds <- snakemake@output$rds
saveRDS(file=out_rds, object=list(variants=variants, 
                                  featuresSubset=featuresSubset,
                                  enrich_final_ls=enrich_final_ls,
                                  recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))
recall_rds <- snakemake@output$recall_at_x_outliers_rds
saveRDS(file=recall_rds, object=list(recallAt10Outliers=recall_10_dt,
                                     recallAt20Outliers=recall_20_dt))
