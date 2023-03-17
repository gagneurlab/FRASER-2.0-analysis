#'---
#' title: Compare jaccard to old version, on optimal q
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset}/jaccard_var_enrich_optQ_venn2.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/variant_extraction/" + config["datasets"][wildcards.dataset]["vcf_group"] + "/rareSplicing_filtered_VariantsTable.tsv.gz"
#'   threads: 5
#'   resources:
#'     - mem_mb: 45000
#'   input:
#'     - fraser2_BB_decoder_newFilt_k20_n1_q95: '`sm config["DATADIR"] + "/{dataset}/minK20_95_minN1/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_BB_decoder_newFilt_k20_n1_q25: '`sm config["DATADIR"] + "/{dataset}/minK20_25_minN1/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_PCA_newFilt_k20_n1_q95: '`sm config["DATADIR"] + "/{dataset}/minK20_95_minN1/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_PCA_newFilt_k20_n1_q25: '`sm config["DATADIR"] + "/{dataset}/minK20_25_minN1/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_BB_decoder_newFilt_k5_n10_q95: '`sm config["DATADIR"] + "/{dataset}/minK5_95_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_BB_decoder_newFilt_k5_n10_q25: '`sm config["DATADIR"] + "/{dataset}/minK5_25_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_PCA_newFilt_k5_n10_q95: '`sm config["DATADIR"] + "/{dataset}/minK5_95_minN10/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser2_PCA_newFilt_k5_n10_q25: '`sm config["DATADIR"] + "/{dataset}/minK5_25_minN10/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
#'     - fraser_newFilt: '`sm config["DATADIR"] + "/{dataset}/optQ__newFilt/FRASER_pval_matrix__0.3__5__1.tsv.gz"`'
#'     - fraser_oldFilt: '`sm config["DATADIR"] + "/{dataset}/optQ/FRASER_pval_matrix__0.3__5__1.tsv.gz"`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/variant_enrichment/{dataset}/optQ_newFilt_variant_enrichment_jaccard_venn2.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/{dataset}/plot_rds/optQ__newFilt/jaccard_enrichment_plots_venn2.Rds"`'
#'     - rds: '`sm config["DATADIR"] + "/{dataset}/plot_rds/optQ__newFilt/jaccard_enrichment_data_venn2.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

# #'     - fraser2_PCA_oldFilt: '`sm config["DATADIR"] + "/{dataset}/optQ/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_PCA_newFilt_k5_n10_q25: '`sm config["DATADIR"] + "/{dataset}/minK5_25_minN10/PCA/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k20_n1_q5: '`sm config["DATADIR"] + "/{dataset}/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n5_q25: '`sm config["DATADIR"] + "/{dataset}/minK5_25_minN5/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n10_q25: '`sm config["DATADIR"] + "/{dataset}/minK5_25_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n5_q40: '`sm config["DATADIR"] + "/{dataset}/minK5_40_minN5/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n10_q40: '`sm config["DATADIR"] + "/{dataset}/minK5_40_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n5_q20: '`sm config["DATADIR"] + "/{dataset}/minK5_20_minN5/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n10_q20: '`sm config["DATADIR"] + "/{dataset}/minK5_20_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n10_q10: '`sm config["DATADIR"] + "/{dataset}/minK5_10_minN10/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'
# #'     - fraser2_BB_decoder_newFilt_k5_n5_q10: '`sm config["DATADIR"] + "/{dataset}/minK5_10_minN5/PCA-BB-Decoder/optQ__newFilt/jaccard/pval_matrix__0.3__5__0.1.tsv.gz"`'

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
source("src/R/enrichment_helper.R")

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ input
variants_file <- snakemake@input$variant_table
dataset <- snakemake@wildcards$dataset
input_files <- list()
input_files[["FRASER_q95_p"]] <- snakemake@input$fraser_oldFilt
input_files[["FRASER_q5_p"]] <- snakemake@input$fraser_newFilt

input_files["FRASER2_PCA_k20_q95_n1_p"] <- snakemake@input$fraser2_PCA_newFilt_k20_n1_q95
input_files["FRASER2_PCA_k20_q25_n1_p"] <- snakemake@input$fraser2_PCA_newFilt_k20_n1_q25
input_files["FRASER2_BB_k20_q95_n1_p"] <- snakemake@input$fraser2_BB_decoder_newFilt_k20_n1_q95
input_files["FRASER2_BB_k20_q25_n1_p"] <- snakemake@input$fraser2_BB_decoder_newFilt_k20_n1_q25
input_files["FRASER2_PCA_k5_q95_n10_p"] <- snakemake@input$fraser2_PCA_newFilt_k5_n10_q95
input_files["FRASER2_PCA_k5_q25_n10_p"] <- snakemake@input$fraser2_PCA_newFilt_k5_n10_q25
input_files["FRASER2_BB_k5_q95_n10_p"] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q95
input_files["FRASER2_BB_k5_q25_n10_p"] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q25


#input_files[["FRASER2_PCA_k20_q95_n1_p"]] <- snakemake@input$fraser2_PCA_oldFilt
#input_files[["FRASER2_BB_k20_q5_n1_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k20_n1_q5
#input_files[["FRASER2_BB_k5_q25_n5_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n5_q25
#input_files[["FRASER2_BB_k5_q25_n10_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q25
#input_files[["FRASER2_BB_k5_q40_n5_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n5_q40
#input_files[["FRASER2_BB_k5_q40_n10_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q40
#input_files[["FRASER2_BB_k5_q20_n5_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n5_q20
#input_files[["FRASER2_BB_k5_q20_n10_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q20
## input_files[["FRASER2_BB_k5_q10_n5_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n5_q10
#input_files[["FRASER2_BB_k5_q10_n10_p"]] <- snakemake@input$fraser2_BB_decoder_newFilt_k5_n10_q10

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

# #+ get overlap between FRASER2 results based on different filterings
# res_all[, `FRASER2_PCA_bothFilterings_p` := pmax(`FRASER2_PCA_newFiltering_p`, `FRASER2_PCA_oldFiltering_p`)]
# res_all[, `FRASER2_bothFits_newFiltering_p` := pmax(`FRASER2_BB_Decoder_newFiltering_p`, `FRASER2_PCA_newFiltering_p`)]
# res_all[(`FRASER2_PCA_oldFiltering_p` > 1e-3 | is.na(`FRASER2_PCA_oldFiltering_p`)) &
#             (`FRASER2_BB_Decoder_newFiltering_p` > 1e-3 | is.na(`FRASER2_BB_Decoder_newFiltering_p`)) &
#             `FRASER2_PCA_newFiltering_p` < 1e-3, `FRASER2_PCA_newFiltering_only_p`:= `FRASER2_PCA_newFiltering_p`]
# res_all[(`FRASER2_PCA_oldFiltering_p` > 1e-3 | is.na(`FRASER2_PCA_oldFiltering_p`)) &
#             (`FRASER2_PCA_newFiltering_p` > 1e-3 | is.na(`FRASER2_PCA_newFiltering_p`)) &
#             `FRASER2_BB_Decoder_newFiltering_p` < 1e-3, `FRASER2_BB_Decoder_newFiltering_only_p`:= `FRASER2_BB_Decoder_newFiltering_p`]
# res_all[(`FRASER2_PCA_newFiltering_p` > 1e-3 | is.na(`FRASER2_PCA_newFiltering_p`)) &
#             (`FRASER2_BB_Decoder_newFiltering_p` > 1e-3 | is.na(`FRASER2_BB_Decoder_newFiltering_p`)) &
#             `FRASER2_PCA_oldFiltering_p` < 1e-3, `FRASER2_PCA_oldFiltering_only_p`:= `FRASER2_PCA_oldFiltering_p`]

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
rrdt[, Method:=factor(Method, levels=c("FRASER_q95", "FRASER_q5",
				       "FRASER2_PCA_k20_q95_n1", "FRASER2_PCA_k20_q25_n1",
				       "FRASER2_BB_k20_q95_n1", "FRASER2_BB_k20_q25_n1",
				       "FRASER2_PCA_k5_q95_n10", "FRASER2_PCA_k5_q25_n10",
				       "FRASER2_BB_k5_q95_n10", "FRASER2_BB_k5_q25_n10"))]
		      #c("FRASER_q95", "FRASER2_PCA_k20_q95_n1", 
                      #                 "FRASER_q5", "FRASER2_BB_k20_q5_n1", 
                      #                 "FRASER2_BB_k5_q25_n5", "FRASER2_BB_k5_q25_n10",
                      #                 "FRASER2_BB_k5_q10_n5", "FRASER2_BB_k5_q10_n10",
                      #                 "FRASER2_BB_k5_q20_n5", "FRASER2_BB_k5_q20_n10",
                      #                 "FRASER2_BB_k5_q40_n5", "FRASER2_BB_k5_q40_n10"))]
                      # levels=c("FRASER_oldFiltering", "FRASER2_PCA_oldFiltering", 
                      #                  "FRASER_newFiltering", "FRASER2_PCA_newFiltering", 
                      #                  "FRASER2_BB_Decoder_newFiltering", "FRASER2_BB_Decoder_newFiltering_only",
                      #                  "FRASER2_PCA_newFiltering_only", "FRASER2_PCA_bothFilterings"))] # "FRASER2_bothFits_newFiltering", "FRASER2_PCA_oldFiltering_only", 
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
dt4cutoffs[, Method:=factor(Method, c("FRASER_q95", "FRASER_q5",
                                      "FRASER2_PCA_k20_q95_n1", "FRASER2_PCA_k20_q25_n1",
                                      "FRASER2_BB_k20_q95_n1", "FRASER2_BB_k20_q25_n1",
                                      "FRASER2_PCA_k5_q95_n10", "FRASER2_PCA_k5_q25_n10",
				      "FRASER2_BB_k5_q95_n10", "FRASER2_BB_k5_q25_n10"))]
			    #c("FRASER_q95", "FRASER2_PCA_k20_q95_n1", 
                            #	    "FRASER_q5", "FRASER2_BB_k20_q5_n1", 
                            #          "FRASER2_BB_k5_q25_n5", "FRASER2_BB_k5_q25_n10",
                            #          "FRASER2_BB_k5_q10_n5", "FRASER2_BB_k5_q10_n10",
                            #          "FRASER2_BB_k5_q20_n5", "FRASER2_BB_k5_q20_n10",
                            #          "FRASER2_BB_k5_q40_n5", "FRASER2_BB_k5_q40_n10"))]
                            # levels=c("FRASER_oldFiltering", "FRASER2_PCA_oldFiltering", 
                            #                  "FRASER_newFiltering", "FRASER2_PCA_newFiltering", 
                            #                  "FRASER2_BB_Decoder_newFiltering", "FRASER2_BB_Decoder_newFiltering_only",
                            #                  "FRASER2_PCA_newFiltering_only", "FRASER2_PCA_bothFilterings"))] # "FRASER2_bothFits_newFiltering", "FRASER2_PCA_oldFiltering_only", 
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

