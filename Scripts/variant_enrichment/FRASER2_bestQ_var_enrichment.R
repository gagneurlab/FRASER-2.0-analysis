#'---
#' title: Compare different encoding dimension for fixed FRASER2 settings
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset_group}/{dataset}/minK{K}_{quantile}_minN{N}/{implementation}/{dPsi}__{minCoverage}__{rho}__{Act}/FRASER2_bestQ_var_enrich_{snptype}.Rds"`'
#'   py:
#'   - |
#'    def get_input_variants(wildcards):
#'     return config["DATADIR"] + "/" + config["datasets"][wildcards.dataset]["vcf_group"] + "/variant_extraction/" + wildcards.snptype + "_filtered_VariantsTable.tsv.gz"
#'   threads: 5
#'   resources:
#'     - mem_mb: 35000
#'   input:
#'     - fraser2_pvals: '`sm expand(config["DATADIR"] + "/{dataset_group}/{dataset}/minK{K}_{quantile}_minN{N}/{implementation}/fixedQ{q}/jaccard/pval_matrix__{dPsi}__{minCoverage}__{rho}__{Act}.tsv.gz", q=config["q"], allow_missing=True)`'
#'     - fds_optQ: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{K}_{quantile}_minN{N}/{implementation}/savedObjects/{dataset}__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'
#'     - variant_table: '`sm get_input_variants`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/{dataset_group}/{dataset}/variant_enrichment/minK{K}_{quantile}_minN{N}/{implementation}/jaccard/{dPsi}__{minCoverage}__{rho}__{Act}/all_q_variant_enrichment_{snptype}.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/{dataset_group}/{dataset}/plot_rds/minK{K}_{quantile}_minN{N}/{implementation}/jaccard/{dPsi}__{minCoverage}__{rho}__{Act}/all_q_enrichment_plots_{snptype}.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
.libPaths("~/R/4.1/FRASER2_BB_loss")
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(FRASER)
source("src/R/enrichment_helper.R")

#+ set nr of threads to use
BPPARAM <- MulticoreParam(snakemake@threads)

#+ input
variants_file <- snakemake@input$variant_table
dataset <- snakemake@wildcards$dataset
# q <- snakemake@wildcards$q
input_files <- snakemake@input$fraser2_pvals
implementation <- snakemake@wildcards$implementation

#+ read in variants
variants <- fread(variants_file)

#+ read input files
res_for_comparision <- lapply(input_files, function(file){
    message(file)
    pval_mat <- fread(file)
    # file_name <- strsplit(basename(file), ".", fixed=TRUE)[[1]][1]
    # filters <- strsplit(file_name, "__", fixed=TRUE)[[1]]
    # method_name <- paste0("FRASER_dPsi", filters[2], "_minN", filters[3], "_rho", filters[4], "_phi", filters[5], "_p")
    method_name <- paste0("FRASER2_", implementation, "_", 
                          gsub("fixedQ", "q", basename(dirname(dirname(file)))), 
                          "_p")
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

#' Dataset: `r dataset`  
#' Number of samples: `r res_all[,uniqueN(subjectID)]`  
#' Number of genes: `r res_all[,uniqueN(geneID)]`  
#' delta psi cutoff: `r snakemake@wildcards$dPsi`
#' min coverage cutoff: `r snakemake@wildcards$minCoverage`
#' rho cutoff: `r snakemake@wildcards$rho`
#' psi type: `r snakemake@wildcards$psiType`  

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
        # q_x <- as.integer(gsub("q", "", strsplit(x, "_", fixed=TRUE)[[1]][2]))
        q_x <- as.integer(gsub("q", "", stringr::str_extract(x, "q[0-9]+")))
        ans <- dt1$dt[, .(cutoff, nRareEvent, total, fraction, nNA, Method=x, 
                          q=q_x,
                          enrichment=dt1$enrichment, 
                          min.ci=dt1$min.ci, max.ci=dt1$max.ci,
                          multiCall=1)]
        ans
    }))
    enrich_final_ls[[paste0(dataset, ": ", et$name)]] <- enrichdt
    
    # ggplots[[paste0('enrichAll_', et$name)]] <-
    #     ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
    #     geom_point() +
    #     facet_grid(rows="multiCall", scales="free_y") + 
    #     geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
    #     coord_flip() +
    #     labs(title=paste0('Enrichment (All, ', dataset, ', ', et$name, ')'))
    
    ggplots[[paste0('enrichAll_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(q, enrichment)) +
        geom_point() + geom_smooth(method="loess", formula=y~x) +
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
        labs(title=paste0('Enrichment (All, ', dataset, ', ', et$name, ')'))
}

#+ enrichment_all, fig.height=8, fig.width=12
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
        # q_x <- as.integer(gsub("q", "", strsplit(x, "_", fixed=TRUE)[[1]][2]))
        q_x <- as.integer(gsub("q", "", stringr::str_extract(x, "q[0-9]+")))
        dt$dt[,.(cutoff, nRareEvent, total, fraction, nNA, Method=x, q=q_x,
                 enrichment=dt$enrichment, min.ci=dt$min.ci, max.ci=dt$max.ci)]
    }))
    enrich_final_ls[[paste0(dataset, ": ", et$name, " + MAF < ", MAF_LIMIT)]] <- enrichdt
    
    # ggplots[[paste0('enrichMAF_', MAF_LIMIT, '_', et$name)]] <-
    #     ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
    #     geom_point() +
    #     geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
    #     coord_flip() +
    #     labs(title=paste0('Enrichment (MAF < ', MAF_LIMIT, ', ', dataset, ', ', et$name, ')'))
    
    ggplots[[paste0('enrichMAF_', MAF_LIMIT, '_', et$name)]] <-
        ggplot(enrichdt[cutoff==TRUE], aes(q, enrichment)) +
        geom_point() + geom_smooth(method="loess", formula=y~x) +
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci), width=.2) +
        labs(title=paste0('Enrichment (MAF < ', MAF_LIMIT, ', ', dataset, ', ', et$name, ')'))
}

#+ enrichment_maf, fig.height=8, fig.width=12
names2plot <- paste0('enrichMAF_', MAF_LIMIT, '_', sapply(enrich_table, "[[", "name"))
etl <- length(enrich_table)/2
order2plot <- rep(1:etl, each=2) + rep(c(0, etl), etl)
grid.arrange(ncol=2, grobs=ggplots[names2plot][order2plot])

#'
#' # FRASER2 q hyper opt plot
#'
#+ plot enc dim search
fds_file <- snakemake@input$fds_optQ
fds <- loadFraserDataSet(file=fds_file)
plotEncDimSearch(fds, type="jaccard")

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
dt4cutoffs

#'
#' The recall plots
#'
#+ recall plots 1
# snptype    <- "rareSplicingVariants"
snptype <- snakemake@wildcards$snptype
maxRank <- 50000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
    rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank)) +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed")
ggplots[[paste0('recall_n=', maxRank)]]

maxRank <- 100000
ggplots[[paste0('recall_n=', maxRank)]] <- plotRecallRankForEnrichment(
    rrdt, maxRank=maxRank, maxPoints=1e4) +
    labs(title=paste(dataset, '\n', snptype, " rank < ", maxRank)) +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=dt4cutoffs, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed")
ggplots[[paste0('recall_n=', maxRank)]]

#+ output
out_file <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_file)
