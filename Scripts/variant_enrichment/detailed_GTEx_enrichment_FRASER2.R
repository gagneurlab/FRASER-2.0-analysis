#'---
#' title: Compare FRASER2 settings, on opt q, for all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_var_enrich_optQ_investigate_recallAt{x}.Rds"`'
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - jaccard_vs_fraser: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER_vs_jaccard_rv_recall_data_{snptype}.Rds", dataset=config["tissues_for_detailed_analysis"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'     - jaccard_vs_fraser_types: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER_types_vs_jaccard_rv_recall_data_{snptype}.Rds", dataset=config["tissues_for_detailed_analysis"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'     - pseudocount: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_pseudocount_recall_at_x_outliers_rho{rho}_{snptype}.Rds", dataset=config["tissues_for_detailed_analysis"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"], rho=[1, 0.5, 0.2, 0.1, 0.05, 0.025, 0.01, 0.001])`'
#'     - filtering: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_filtering_recall_at_x_outliers_pc0.1_rho1_delta{delta}_{snptype}.Rds", dataset=config["tissues_for_detailed_analysis"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"], delta=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7])`'
#'     - filtering_info: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/filtering_info.tsv"`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/all_rv_recallAt{x}.html"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/parameter_optimization_rv_recallAt{x}_ggplots.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

# #'     - goodness_of_fit: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_goodnessOfFit_recall_at_x_outliers_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice"])`'

# #'     - bestQ_PCA: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/minK20_95_minN1/PCA__pc01/jaccard/0.3__5__0.1__FALSE/all_q_enrichment_plots_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareMMSplice", "rareSpliceAI"])`'
# #'     - bestQ_BB: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/minK20_95_minN1/PCA-BB-Decoder/jaccard/0.3__5__0.1/all_q_enrichment_plots_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareMMSplice", "rareSpliceAI"])`'
# , "rareIntronic"
# #'     - all_done: '`sm config["htmlOutputPath"] + "/GTEx_v8/variant_enrichment/detailed_FRASER2_enrichment_investigation.done"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ load libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
recall_at <- as.numeric(snakemake@wildcards$x)

#+ create list of ggplots to save
ggplots <- list()

# #+ read in goodness of fit recall data
# goodnessOfFit_recall_dt <- rbindlist(lapply(snakemake@input$goodness_of_fit, function(file){
#     recall_dt <- readRDS(file)[["recallAt20Outliers"]]
#     t <- basename(dirname(dirname(dirname(file))))
#     recall_dt[, tissue:=t]
#     snptype <- grep("rare(.)+", 
#                     strsplit(basename(file), "_", fixed=T)[[1]], 
#                     value=T)
#     snptype <- gsub(".Rds", "", snptype)
#     recall_dt[, snptype:=snptype]
#     return(recall_dt)
# }))

#+ read in pseudocount recall data
recall_dt <- rbindlist(bplapply(snakemake@input$pseudocount, function(file){
    recall_dt <- readRDS(file)[[paste0("recallAt", recall_at, "Outliers")]]
    t <- basename(dirname(dirname(dirname(file))))
    recall_dt[, tissue:=t]
    snptype <- grep("rare(.)+",
                    strsplit(basename(file), "_", fixed=T)[[1]],
                    value=T)
    snptype <- gsub(".Rds", "", snptype)
    recall_dt[, snptype:=snptype]
    return(recall_dt)
}))
recall_dt[, rho:=factor(rho)]
recall_dt[, pc:=factor(pc)]

recall_dt[snptype == "rareSpliceAI", snptype:="rare SpliceAI"]
recall_dt[snptype == "rareMMSplice", snptype:="rare MMSplice"]
recall_dt[snptype == "rareAbSplice", snptype:="rare AbSplice"]
recall_dt[snptype == "rareSplicing", snptype:="rare splice\nsite vicinity\n(VEP)"]
recall_dt[, snptype:=factor(snptype, levels=c(
                                        "rare splice\nsite vicinity\n(VEP)",
                                        "rare MMSplice",
                                        "rare SpliceAI",
                                        "rare AbSplice"))]

#+ opt_pc_for_goodnessOfFit, fig.width=21, fig.height=14
g_joined_pc_rho_xPc <- ggplot(recall_dt, aes(pc, recall, fill=pc)) +
    facet_grid(snptype ~ rho, labeller=label_both) +
    geom_boxplot() +
    labs(x="Pseudocount", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Blues") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size)) 
g_joined_pc_rho_xPc

#+ opt_pc_for_fixed_rho1, fig.width=21, fig.height=14
g_pc <- ggplot(recall_dt[rho == 1], aes(pc, recall, fill=pc)) +
    facet_grid(~ snptype) +
    geom_boxplot() +
    labs(x="Pseudocount", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Blues") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))  
g_pc

#+ opt_rho_for_pc, fig.width=21, fig.height=14
g_joined_pc_rho_xRho <- ggplot(recall_dt, aes(rho, recall, fill=rho)) +
    facet_grid(snptype ~ pc, labeller=label_both) +
    geom_boxplot() +
    labs(x="rho cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) +  
    scale_fill_brewer(palette="Greens") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_joined_pc_rho_xRho

#+ opt_rho_for_fixed_pc0.1, fig.width=21, fig.height=14
g_rho <- ggplot(recall_dt[pc == 0.1], aes(rho, recall, fill=rho)) +
    facet_grid(~ snptype) +
    geom_boxplot() +
    labs(x="rho cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) +  
    scale_fill_brewer(palette="Greens") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_rho

ggplots[["pc_opt_fixed_rho"]] <- g_pc
ggplots[["rho_opt_fixed_pc"]] <- g_rho
ggplots[["joined_pc_rho_opt_xPc"]] <- g_joined_pc_rho_xPc
ggplots[["joined_pc_rho_opt_xRho"]] <- g_joined_pc_rho_xRho

#+ optimal_combination_pc_rho, fig.width=12, fig.height=5
# ggplot(recall_dt[, .SD[which.max(recall), (pc, rho)], by="tissue,snptype"],
ggplot(recall_dt[, .SD[which.max(recall), paste0("pc", pc, "_rho", rho)], by="tissue,snptype"],
       aes(V1)) +
    geom_bar() +
    facet_wrap(~snptype) +
    labs(x="Pseudocount and rho cutoff combination", y="number of tissues in which cutoff combination is optimal") + 
    theme(axis.text = element_text(angle=45, hjust=1, size=font_size),
          axis.title=element_text(size=font_size),
          strip.text = element_text(size = font_size)) 


# #+ goodnessOfFit_recall_across_tissues, fig.width=12, fig.height=5
# ggplot(goodnessOfFit_recall_dt, aes(factor(rho), recall)) +
#     facet_wrap(~ snptype) +
#     geom_boxplot() +
#     labs(x="rho cutoff", y="Recall at 20 outliers/sample")
# 
# #+ goodnessOfFit_recall_best_rho_per_tissue, fig.width=12, fig.height=5
# ggplot(goodnessOfFit_recall_dt[, .SD[which.max(recall),rho], by="tissue,snptype"],
#        aes(factor(V1))) + 
#     geom_bar() + 
#     facet_wrap(~snptype) + 
#     labs(x="rho cutoff", y="number of tissues in which rho cutoff is optimal")
# 
# ggscatter(dcast(goodnessOfFit_recall_dt[, .(Method, recall, tissue, snptype)], tissue + snptype ~ Method, value.var="recall"),
#           x="rho_1_0", y="rho_0_1", 
#           xlab="Recall at 20/outliers per sample \nfor rho_cutoff=1", 
#           ylab="Recall at 20/outliers per sample \nfor rho_cutoff=0.1") +
#     facet_wrap(~ snptype) +
#     geom_abline(intercept=0, slope=1, linetype="dotted") 

# #+ pseudocount_recall_across_tissues, fig.width=12, fig.height=5
# ggplot(pseudocount_recall_dt, aes(factor(pc), recall)) +
#     facet_wrap(~ snptype) +
#     geom_boxplot() +
#     labs(x="Pseudocount", y="Recall at 20 outliers/sample")

#+ read in filtering recall data
filtering_dt <- rbindlist(bplapply(snakemake@input$filtering, function(file){
    recall_dt <- readRDS(file)[[paste0("recallAt", recall_at, "Outliers")]]
    t <- basename(dirname(dirname(dirname(file))))
    recall_dt[, tissue:=t]
    snptype <- grep("rare(.)+",
                    strsplit(basename(file), "_", fixed=T)[[1]],
                    value=T)
    snptype <- gsub(".Rds", "", snptype)
    recall_dt[, snptype:=snptype]
    return(recall_dt)
}))
filtering_dt[, rho:=factor(rho)]
filtering_dt[, pc:=factor(pc)]
filtering_dt[, delta:=factor(delta)]
filtering_dt[, filter_settings:=paste0("minK", k, "_q", q, "_N", n)]
filtering_dt[, filter_settings:=factor(filter_settings)]
filtering_dt[, k:=factor(k)]
filtering_dt[, q:=factor(q)]
filtering_dt[, n:=factor(n)]
filtering_dt[snptype == "rareSpliceAI", snptype:="rare SpliceAI"]
filtering_dt[snptype == "rareAbSplice", snptype:="rare AbSplice"]
filtering_dt[snptype == "rareMMSplice", snptype:="rare MMSplice"]
filtering_dt[snptype == "rareSplicing", snptype:="rare splice\nsite vicinity\n(VEP)"]
filtering_dt[, snptype:=factor(snptype, levels=c("rare splice\nsite vicinity\n(VEP)",
                                                "rare MMSplice",
                                                "rare SpliceAI", 
                                                "rare AbSplice"
                                              ))]


#+ joined_delta_filtering, fig.width=21, fig.height=14
g_joined_delta_filter <- ggplot(filtering_dt, aes(delta, recall, fill=delta)) +
    facet_grid(snptype ~ filter_settings, labeller=label_both) +
    geom_boxplot() +
    labs(x="Delta jaccard cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Oranges") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_joined_delta_filter

#+ opt_delta_for_fixed_filtering_n10q25, fig.width=21, fig.height=14
g_delta_opt_filter <- ggplot(filtering_dt[k == 20 & n == 10 & q == 25,], aes(delta, recall, fill=delta)) +
    facet_grid(~ snptype) +
    geom_boxplot() +
    labs(x="Delta jaccard cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Oranges") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_delta_opt_filter

#+ opt_delta_for_fixed_filtering_n1q95, fig.width=21, fig.height=14
g_delta_fixed_filter_n1q95 <- ggplot(filtering_dt[k == 20 & n == 1 & q == 95,], aes(delta, recall, fill=delta)) +
    facet_grid(~ snptype) +
    geom_boxplot() +
    labs(x="Delta jaccard cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Oranges") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_delta_fixed_filter_n1q95

#+ joined_delta_filtering_rareSpliceAI, fig.width=21, fig.height=14
g_joined_delta_filter_spliceAI <- ggplot(filtering_dt[snptype == "rare SpliceAI" & k == 20,], 
                                         aes(delta, recall, fill=delta)) +
    facet_grid(q ~ n, labeller=label_both) +
    geom_boxplot() +
    labs(x="Delta jaccard cutoff", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Oranges") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(font_size),
          strip.text = element_text(size = font_size))
g_joined_delta_filter_spliceAI

#+ joined_filtering_delta_spliceAI, fig.width=21, fig.height=14
g_joined_filter_delta_spliceAI <- ggplot(filtering_dt[snptype == "rare SpliceAI" & k == 20 & delta %in% c(0.1, 0.2, 0.3),], 
                                         aes(Method, recall, fill=Method)) +
    facet_grid(~delta, labeller=label_both) +
    geom_boxplot() +
    labs(x="Filtering setting", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), 
          axis.text.y=element_text(size=font_size),
          axis.text.x=element_text(size=font_size, angle=45, hjust=0.5, vjust=0.5),
          strip.text = element_text(size = font_size))
    
g_joined_filter_delta_spliceAI

#+ filtering_delta0.3, fig.width=21, fig.height=14
g_filter_delta0.3 <- ggplot(filtering_dt[delta == 0.3 & k == 20,], aes(n, recall, fill=n)) +
    facet_grid(snptype ~ q, labeller=label_both) +
    geom_boxplot() +
    labs(x="Minimal N in at least q% samples", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_filter_delta0.3

#+ filtering_delta0.2, fig.width=21, fig.height=14
g_filter_delta0.2 <- ggplot(filtering_dt[delta == 0.2 & k == 20,], aes(n, recall, fill=n)) +
    facet_grid(snptype ~ q, labeller=label_both) +
    geom_boxplot() +
    labs(x="Minimal N in at least q% samples", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_filter_delta0.2

#+ filtering_delta0.3_rareSpliceAI, fig.width=15, fig.height=5
g_filter_delta0.3_spliceAI <- ggplot(filtering_dt[snptype == "rare SpliceAI" & delta == 0.3 & k == 20,], aes(n, recall, fill=n)) +
    facet_grid( ~ q, labeller=label_both) +
    geom_boxplot() +
    labs(x="Minimal N in at least q% samples", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_filter_delta0.3_spliceAI

#+ filtering_delta0.2_rareSpliceAI, fig.width=15, fig.height=5
g_filter_delta0.2_spliceAI <- ggplot(filtering_dt[snptype == "rare SpliceAI" & delta == 0.2 & k == 20,], aes(n, recall, fill=n)) +
    facet_grid( ~ q, labeller=label_both) +
    geom_boxplot() +
    labs(x="Minimal N in at least q% samples", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_filter_delta0.2_spliceAI

#+ filtering_delta0.1_rareSpliceAI, fig.width=15, fig.height=5
g_filter_delta0.1_spliceAI <- ggplot(filtering_dt[snptype == "rare SpliceAI" & delta == 0.1 & k == 20,], aes(n, recall, fill=n)) +
    facet_grid( ~ q, labeller=label_both) +
    geom_boxplot() +
    labs(x="Minimal N in at least q% samples", y=paste0("Recall at ", recall_at, " outliers/sample")) + 
    scale_fill_brewer(palette="Purples") + 
    guides(fill="none") +
    theme_bw() +
    theme(axis.title=element_text(size=font_size), axis.text=element_text(size=font_size),
          strip.text = element_text(size = font_size))
g_filter_delta0.1_spliceAI

ggplots[["filter_opt_fixed_delta0.3"]] <- g_filter_delta0.3
ggplots[["filter_opt_fixed_delta0.2"]] <- g_filter_delta0.2
ggplots[["filter_opt_fixed_delta0.3+spliceAIvars"]] <- g_filter_delta0.3_spliceAI
ggplots[["filter_opt_fixed_delta0.2+spliceAIvars"]] <- g_filter_delta0.2_spliceAI
ggplots[["filter_opt_fixed_delta0.1+spliceAIvars"]] <- g_filter_delta0.1_spliceAI
ggplots[["joined_delta_filter_opt"]] <- g_joined_delta_filter
ggplots[["joined_delta_filter_opt_spliceAIvars"]] <- g_joined_delta_filter_spliceAI
ggplots[["joined_filter_delta_opt_spliceAIvars"]] <- g_joined_filter_delta_spliceAI
ggplots[["delta_filter_opt"]] <- g_delta_opt_filter
ggplots[["delta_filter_fixed_n1q95"]] <- g_delta_fixed_filter_n1q95

#+ filtering_nrJunc, fig.width=12, fig.height=5
filter_info_dt <- fread(snakemake@input$filtering_info)
setnames(filter_info_dt, "njuncs", "Introns")
setnames(filter_info_dt, "ngenes", "Genes")
g_filter_stats <- ggplot(melt(filter_info_dt, id.vars=c("tissue", "k", "q", "n", "nsamples"), 
            variable.name="type", value.name="nr"), 
       aes(factor(n), nr, fill=factor(n))) +
    facet_grid(type ~ q, labeller=labeller(q=label_both), scales="free_y") +
    geom_boxplot() +
    scale_fill_brewer(palette = "Purples") +
    guides(fill="none") +
    # scale_y_log10(limits=c(1, max(filter_info_dt$Introns + 1e4))) +
    # scale_y_log10() +
    # annotation_logticks(sides="l") +
    labs(x="Minimal N in at least q% samples", y="Number after filtering") +
    theme_bw() +
    theme(axis.title=element_text(size=12), axis.text=element_text(size=12),
          strip.text = element_text(size = 12))
g_filter_stats

ggplots[["filter_stats"]] <- g_filter_stats

#+ output
out_gg_rds <- snakemake@output$gg_rds
saveRDS(ggplots, file=out_gg_rds)