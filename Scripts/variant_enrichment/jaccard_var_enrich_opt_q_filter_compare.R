#'---
#' title: Compare jaccard to old version, on optimal q
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset}/jaccard_var_enrich_optQ_venn.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/variant_extraction/" + config["datasets"][wildcards.dataset]["vcf_group"] + "/rareSplicing_filtered_VariantsTable.tsv.gz"
#'   threads: 5
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fraser2_delta03_newFilt: '`sm config["DATADIR"] + "/{dataset}/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_delta03_oldFilt: '`sm config["DATADIR"] + "/{dataset}/optQ/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_delta03_BB_decoder_newFilt: '`sm config["DATADIR"] + "/{dataset}/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser: '`sm config["DATADIR"] + "/{dataset}/optQ__newFilt/FRASER_pval_matrix__0.3__5__1.tsv.gz"`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/variant_enrichment/{dataset}/optQ_newFilt_variant_enrichment_jaccard_venn.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/{dataset}/plot_rds/optQ__newFilt/jaccard_enrichment_plots_venn.Rds"`'
#'     - rds: '`sm config["DATADIR"] + "/{dataset}/plot_rds/optQ__newFilt/jaccard_enrichment_data_venn.Rds"`'
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
source("src/R/enrichment_helper.R")

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ input
variants_file <- snakemake@input$variant_table
dataset <- snakemake@wildcards$dataset
input_files <- list()
input_files[["FRASER_p"]] <- snakemake@input$fraser
input_files[["FRASER2_delta>=0.3_newFiltering_p"]] <- snakemake@input$fraser2_delta03_newFilt
input_files[["FRASER2_delta>=0.3_oldFiltering_p"]] <- snakemake@input$fraser2_delta03_oldFilt
input_files[["FRASER2_delta>=0.3_BB-decoder_newFiltering_p"]] <- snakemake@input$fraser2_delta03_BB_decoder_newFilt

#+ read in variants
variants <- fread(variants_file)

#+ read input files
res_for_comparision <- lapply(names(input_files), function(method_name){
    pval_mat <- fread(input_files[[method_name]])
    dt <- melt(pval_mat, id.vars="geneID", value.name=method_name, variable.name="sampleID")
    setkey(dt, geneID, sampleID)
    return(dt)
})

#+ merge to one res table
res_all <- Reduce(function(...) merge(..., all = TRUE), res_for_comparision)
setnames(res_all, "sampleID", "subjectID")
res_all[, subjectID:= gsub("^", "GTEX-", gsub("-.*", "", gsub("GTEX-", "", res_all$subjectID)))]
res_all[, tissue:=dataset]
res_all

#+ get overlap between FRASER2 results based on different filterings
res_all[, `FRASER2_delta>=0.3_bothFilterings_p` := pmax(`FRASER2_delta>=0.3_newFiltering_p`, `FRASER2_delta>=0.3_oldFiltering_p`)]
res_all[(`FRASER2_delta>=0.3_oldFiltering_p` > 1e-3 | is.na(`FRASER2_delta>=0.3_oldFiltering_p`)) & 
            `FRASER2_delta>=0.3_newFiltering_p` < 1e-3, `FRASER2_delta>=0.3_newFiltering_only_p`:= `FRASER2_delta>=0.3_newFiltering_p`]
res_all[(`FRASER2_delta>=0.3_newFiltering_p` > 1e-3 | is.na(`FRASER2_delta>=0.3_newFiltering_p`)) & 
            `FRASER2_delta>=0.3_oldFiltering_p` < 1e-3, `FRASER2_delta>=0.3_oldFiltering_only_p`:= `FRASER2_delta>=0.3_oldFiltering_p`]

#+ get overlap between FRASER2 results based on different filterings (with BB decoder results)
res_all[, `FRASER2_delta>=0.3_bothFilterings_BB-decoder_p` := pmax(`FRASER2_delta>=0.3_BB-decoder_newFiltering_p`, `FRASER2_delta>=0.3_oldFiltering_p`)]
res_all[(`FRASER2_delta>=0.3_oldFiltering_p` > 1e-3 | is.na(`FRASER2_delta>=0.3_oldFiltering_p`)) & 
            `FRASER2_delta>=0.3_BB-decoder_newFiltering_p` < 1e-3, `FRASER2_delta>=0.3_BB-decoder_newFiltering_only_p`:= `FRASER2_delta>=0.3_BB-decoder_newFiltering_p`]
res_all[(`FRASER2_delta>=0.3_BB-decoder_newFiltering_p` > 1e-3 | is.na(`FRASER2_delta>=0.3_BB-decoder_newFiltering_p`)) & 
            `FRASER2_delta>=0.3_oldFiltering_p` < 1e-3, `FRASER2_delta>=0.3_oldFiltering_only(BB)_p`:= `FRASER2_delta>=0.3_oldFiltering_p`]

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
    list(name="Pval 1e-3", methods=methodP, isZscore=FALSE, cutOff=1e-3),
    #list(name="Pval 1e-4", methods=methodP, isZscore=FALSE, cutOff=1e-4),
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
        labs(title=paste0('Enrichment (All, ', dataset, ', ', et$name, ')'))
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
        labs(title=paste0('Enrichment (MAF < ', MAF_LIMIT, ', ', dataset, ', ', et$name, ')'))
}

#+ enrichment_maf, fig.height=8, fig.width=16
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
rrdt <- rbindlist(bplapply(methods2plot, dt=featuresSubset, BPPARAM=BPPARAM,
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
rrdt[, Method:=factor(Method, levels=c("FRASER", "FRASER2_delta>=0.3_oldFiltering", "FRASER2_delta>=0.3_newFiltering", "FRASER2_delta>=0.3_BB-decoder_newFiltering", 
                                       "FRASER2_delta>=0.3_newFiltering_only", "FRASER2_delta>=0.3_oldFiltering_only", "FRASER2_delta>=0.3_bothFilterings",
                                       "FRASER2_delta>=0.3_BB-decoder_newFiltering_only", "FRASER2_delta>=0.3_oldFiltering_only(BB)", "FRASER2_delta>=0.3_bothFilterings_BB-decoder"))]
rrdt

#'
#' Merge with cutoffs from enrichment
#'
dt_tmp <- rbindlist(lapply(names(enrich_final_ls[1:5]), function(x){
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
dt4cutoffs[, Method:=factor(Method, levels=c("FRASER", "FRASER2_delta>=0.3_oldFiltering", "FRASER2_delta>=0.3_newFiltering", "FRASER2_delta>=0.3_BB-decoder_newFiltering", 
                                             "FRASER2_delta>=0.3_newFiltering_only", "FRASER2_delta>=0.3_oldFiltering_only", "FRASER2_delta>=0.3_bothFilterings",
                                             "FRASER2_delta>=0.3_BB-decoder_newFiltering_only", "FRASER2_delta>=0.3_oldFiltering_only(BB)", "FRASER2_delta>=0.3_bothFilterings_BB-decoder"))]
dt4cutoffs

#'
#' The recall plots
#'
#+ recall plots 1, fig.height=6, fig.width=10
snptype    <- "rareSplicingVariants"
maxRank <- 50000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
    rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
         x="Rank", y="Recall") +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
    scale_color_brewer(palette="Paired") + 
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "P-value cutoff", "Cutoff")) ) +
    theme_bw() + theme(text=element_text(size=16))
ggplots[[paste0('recall_n=', maxRank)]]

#+ recall plots 2, fig.height=6, fig.width=10
maxRank <- 100000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
    rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
         x="Rank", y="Recall") +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
    scale_color_brewer(palette="Paired") + 
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "P-value cutoff", "Cutoff")) ) +
    theme_bw() + theme(text=element_text(size=16))
ggplots[[paste0('recall_n=', maxRank)]]

#+ recall plots 3, fig.height=6, fig.width=10
ggplots[[paste0('recall_n=', maxRank, "_xlog")]] <- plotRecallRankForEnrichment(
    rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank), 
         x="Rank", y="Recall") +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
    scale_color_brewer(palette="Paired") + 
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(dt4cutoffs[,Type == "P-value"]), "P-value cutoff", "Cutoff")) ) +
    theme_bw() + theme(text=element_text(size=16)) +
    scale_x_log10()
ggplots[[paste0('recall_n=', maxRank, "_xlog")]]

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)
out_rds <- snakemake@output$rds
saveRDS(file=out_rds, object=list(variants=variants, 
                                  featuresSubset=featuresSubset,
                                  enrich_final_ls=enrich_final_ls,
                                  recallData=rbind(rrdt, dt4cutoffs, fill=TRUE)))

