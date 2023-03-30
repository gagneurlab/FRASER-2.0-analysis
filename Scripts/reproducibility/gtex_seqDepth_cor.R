#'---
#' title: seq depth vs nr of outliers (FRASER1 vs FRASER2)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/gtex_seqDepthCor_delta{delta}.Rds"`'
#'   threads: 5
#'   resources:
#'    - mem_mb: 32000
#'   input:
#'    - bam_coverage_tsv: '`sm expand(config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode34/outrider/{tissue}/bam_coverage.tsv", tissue=config["tissues_for_reproducibility"])`'
#'    - res_FRASER2: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK20_25_minN10/PCA__pc0.1/{dataset}/optQ__newFilt/delta{{delta}}/results_gene.tsv", dataset=config["tissues_for_reproducibility"])`'
#'    - res_FRASER1: '`sm expand("config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{tissue}_old_filter/results.tsv", tissue=config["tissues_for_reproducibility"])`'
#'    - res_SPOT_LeafcutterMD: '`sm expand(config["DATADIR"] + "/GTEx_v8/fraser2_improvements/{tissue}/seq_depth_cor_spot_leafcutterMD.tsv", tissue=config["tissues_for_reproducibility"])`'
#'   output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/FRASER1_vs_FRASER2/gtex_seqDepth_cor_delta{delta}.html"`'
#'    - gg_single_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/seq_depth_cor_single_ggplot.Rds"`'
#'    - gg_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/seq_depth_cor_all_ggplot.Rds"`'
#'    - skin_seqDepth_cor: '`sm config["htmlOutputPath"] + "/GTEx_v8/FRASER1_vs_FRASER2/delta{delta}/skin_seqDepth_vs_nrOultiers.png"`'
#'    - skin_seqDepth_cor_log: '`sm config["htmlOutputPath"] + "/GTEx_v8/FRASER1_vs_FRASER2/delta{delta}/skin_seqDepth_vs_nrOultiers_log.png"`'
#'    - skin_seqDepth_cor_noOutlier: '`sm config["htmlOutputPath"] + "/GTEx_v8/FRASER1_vs_FRASER2/delta{delta}/skin_seqDepth_vs_nrOultiers_noOutlier.png"`'
#'    - all_seqDepth_cor: '`sm config["htmlOutputPath"] + "/GTEx_v8/FRASER1_vs_FRASER2/delta{delta}/all_seqDepth_vs_nrOultiers.png"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
register(MulticoreParam(snakemake@threads))
library(data.table)
library(ggplot2)
library(ggpubr)

#+ read in result tables and coverage information
coverage_dt <- rbindlist(bpmapply(snakemake@input$bam_coverage_tsv, 
                                  snakemake@input$res_FRASER1,
                                  # res_files_f1,
                                  snakemake@input$res_FRASER2, 
                                  # snakemake@input$res_FRASER1_newFilter, snakemake@input$res_FRASER2_newFilter, 
    FUN=function(cov, res_f1, res_f2){ # , res_f1_newFilt, res_f2_newFilt){
        # get tissue name
        t <- basename(dirname(cov))
        # for each sample, get seq depth
        bam_coverage <- fread(cov)
        # get nr of aberrant events (psi3/5/theta)
        res_FRASER1 <- fread(res_f1)
        res_FRASER1[, nrOutPerSample:=.N, by="sampleID"]
        # res_FRASER1_newFilt <- fread(res_f1_newFilt)
        # res_FRASER1_newFilt[, nrOutPerSample:=.N, by="sampleID"]
        # get nr of aberrant events (jaccard metric)
        res_FRASER2 <- fread(res_f2)
        res_FRASER2[, nrOutPerSample:=.N, by="sampleID"]
        # res_FRASER2_newFilt <- fread(res_f2_newFilt)
        # res_FRASER2_newFilt[, nrOutPerSample:=.N, by="sampleID"]
        
        # merge coverage with nr outliers
        coverage_dt <- merge(res_FRASER1[, .N, by="sampleID"], bam_coverage, by="sampleID")
        setnames(coverage_dt, "N", "FRASER")
        coverage_dt <- merge(coverage_dt, res_FRASER2[, .N, by="sampleID"])
        setnames(coverage_dt, "N", "FRASER2")
        # coverage_dt <- merge(coverage_dt, res_FRASER1_newFilt[, .N, by="sampleID"])
        # setnames(coverage_dt, "N", "FRASER_newFilt")
        # coverage_dt <- merge(coverage_dt, res_FRASER2_newFilt[, .N, by="sampleID"])
        # setnames(coverage_dt, "N", "FRASER2_newFilt")
        coverage_dt <- melt(coverage_dt, id.vars=c("sampleID", "record_count"), value.name="nr_outliers", variable.name="method")    
        # coverage_dt[grepl("new", method), filtering:="new"]
        # coverage_dt[!grepl("new", method), filtering:="old"]
        # coverage_dt[grepl("FRASER2", method), method:="FRASER2"]
        # coverage_dt[!grepl("FRASER2", method), method:="FRASER"]
        coverage_dt[, tissue:=t]
        return(coverage_dt)
}, SIMPLIFY=FALSE))
coverage_dt[, nr_outliers_plus_one := nr_outliers+1]

# read in results of spot and leafcutterMD
cov_spot_lmd <- rbindlist(lapply(snakemake@input$res_SPOT_LeafcutterMD, fread))
coverage_dt <- rbind(coverage_dt, cov_spot_lmd)

DT::datatable(
    coverage_dt,
    caption = 'Correlation with seq depth for GTEx tissues',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

# scatterplot seq depth correlation for one tissue (Skin)
g_cor_skin1 <- ggscatter(coverage_dt[tissue == "Skin_-_Not_Sun_Exposed_Suprapubic"], 
                        x="record_count", y="nr_outliers_plus_one", 
                        size=1.5,
                        add="reg.line", 
                        conf.int=TRUE, 
                        add.params=list(color="blue", fill="lightgray"), 
                        cor.coef=TRUE, 
                        cor.method="spearman",
                        cor.coef.size=3) + 
    facet_wrap(~method, nrow=2) +
    # facet_grid(filtering~method, labeller=label_both) + 
    labs(x="Sequencing depth", y="Splicing outliers + 1") + 
    theme(text=element_text(size=14))
g_cor_skin1

# show instead in log scale
g_cor_skin2 <- ggscatter(coverage_dt[tissue == "Skin_-_Not_Sun_Exposed_Suprapubic"], 
                         x="record_count", y="nr_outliers_plus_one", 
                         size=1.5,
                         add="reg.line", 
                         conf.int=TRUE, 
                         add.params=list(color="blue", fill="lightgray"), 
                         cor.coef=TRUE, 
                         cor.method="spearman",
                         cor.coef.size=3) + 
    facet_wrap(~method, nrow=2, scales="free_y") +
    # facet_grid(filtering~method, labeller=label_both, scales="free_y") + 
    scale_y_log10(limits=c(1, 10000)) + 
    labs(x="Sequencing depth", y="Splicing outliers + 1") + 
    theme(text=element_text(size=14))
g_cor_skin2

# remove outlier in scatterplot seq depth correlation
samples_many_outliers <- coverage_dt[tissue == "Skin_-_Not_Sun_Exposed_Suprapubic" & method == "FRASER2" & nr_outliers > 500, sampleID]
g_cor_skin3 <- ggscatter(coverage_dt[tissue == "Skin_-_Not_Sun_Exposed_Suprapubic" & !sampleID %in% samples_many_outliers], x="record_count", y="nr_outliers_plus_one", 
                         size=1.5,
                         add="reg.line", 
                   conf.int=TRUE, 
                   add.params=list(color="blue", fill="lightgray"), 
                   cor.coef=TRUE, 
                   cor.method="spearman",
                   cor.coef.size=3) + 
    facet_wrap(~method, nrow=2, scales="free_y") +
    # facet_grid(filtering~method, labeller=label_both, scales="free_y") + 
    labs(x="Sequencing depth", y="Splicing outliers + 1") + 
    theme(text=element_text(size=14))
g_cor_skin3


# save plots
ggsave(g_cor_skin1, filename=snakemake@output$skin_seqDepth_cor,
       width=8, height=8)
ggsave(g_cor_skin2, filename=snakemake@output$skin_seqDepth_cor_log,
       width=8, height=8)
ggsave(g_cor_skin3, filename=snakemake@output$skin_seqDepth_cor_noOutlier,
       width=8, height=8)

#+ save ggplot (single tissue)
saveRDS(g_cor_skin2, file=snakemake@output$gg_single_rds)

# scatterplot seq depth correlation over tissues
cor_dt <- coverage_dt[, cor(record_count, nr_outliers, method="spearman"), by="tissue,method"] #,filtering"]
setnames(cor_dt, "V1", "cor_coef")
cor_dt <- dcast(cor_dt, "tissue ~ method", value.var="cor_coef")
# cor_dt <- dcast(cor_dt, "tissue + filtering ~ method", value.var="cor_coef")
g_cor_all <- ggscatter(cor_dt, x="FRASER", y="FRASER2", size=1.5) + 
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    labs(x="FRASER spearman correlation coefficient", y="FRASER2 spearman correlation coefficient") + 
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    # facet_wrap(~filtering, labeller=label_both) +
    theme(text=element_text(size=14))
g_cor_all

# save plot
ggsave(g_cor_all, filename=snakemake@output$all_seqDepth_cor,
       width=6, height=6)

#+ save ggplot
saveRDS(g_cor_all, file=snakemake@output$gg_all_rds)

