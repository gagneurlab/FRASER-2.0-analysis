#'---
#' title: seq depth vs nr of outliers (old vs jaccard version)
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/10_seqDepthCor_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 1
#'  resources:
#'   - mem_mb: 10000
#'  input:
#'   - bam_coverage_tsv: '`sm config["mito_bam_coverage_tsv"]`'
#'   - resultTableJunction_old: '`sm config["mito_fraser1_results"] + "/results/" + config["mito_annotation"] + "/fraser/fib/results_per_junction.tsv"`'
#'   - resultTableJunction_jaccard: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_per_junction_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - filtered_fds_FRASER: '/s/project/prokisch/processed_results/aberrant_splicing/datasets/savedObjects/fib--gencode34/fds-object.RDS'
#'   - filtered_fds_FRASER2: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                "/fds-object.RDS"`'
#'  output:
#'   - bam_coverage: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + "/seqDepth_vs_nrOultiers.png"`'
#'   - total_junc_cov: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + "/totalJunctionCoverage_vs_nrOutliers.png"`'
#'   - total_junc_cov_filtered: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + "/totalJunctionCoverage_vs_nrOutliers_filtered.png"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
library(ggplot2)
library(ggpubr)

# for each sample, get seq depth
bam_coverage <- fread(snakemake@input$bam_coverage_tsv)
hist(bam_coverage$record_count)

# get nr of aberrant events (psi3/5/theta)
res_old <- fread(snakemake@input$resultTableJunction_old)
res_old[, nrOutPerSample:=.N, by="sampleID"]

# get nr of aberrant events (jaccard metric)
res_jacc <- fread(snakemake@input$resultTableJunction_jaccard)
res_jacc[, nrOutPerSample:=.N, by="sampleID"]

# merge coverage with nr outliers
coverage_dt <- merge(res_old[, .N, by="sampleID"], bam_coverage, by="sampleID")
setnames(coverage_dt, "N", "FRASER")
coverage_dt <- merge(coverage_dt, res_jacc[, .N, by="sampleID"])
setnames(coverage_dt, "N", "FRASER2")
coverage_dt <- melt(coverage_dt, id.vars=c("sampleID", "record_count"), value.name="nr_outliers", variable.name="metric")

# scatterplot
g_cor <- ggscatter(coverage_dt, x="record_count", y="nr_outliers", 
                  add="reg.line", 
                  conf.int=TRUE, 
                  add.params=list(color="blue", fill="lightgray"), 
                  cor.coef=TRUE, 
                  cor.method="spearman",
                  cor.coef.size=6) + 
    facet_wrap(~metric, nrow=2) + 
    labs(x="Sequencing depth", y="Splicing outliers") + 
    theme(text=element_text(size=14))
# g_cor

# save plot
ggsave(g_cor, filename=snakemake@output$bam_coverage,
       width=8, height=8)

# load fds
fds_raw <- loadFraserDataSet(file=snakemake@config$raw_mito_fds)
fds_old <- loadFraserDataSet(file=snakemake@input$filtered_fds_FRASER)
fds_jacc <- loadFraserDataSet(file=snakemake@input$filtered_fds_FRASER2)

# get seq depth as total junction coverage (sum of K(fds) per sample)
# raw fds
total_j_cov_raw <- colSums(K(fds_raw, "psi5"))
total_j_cov_raw <- data.table(sampleID=names(total_j_cov_raw), total_junc_cov=total_j_cov_raw)

# filtered fds
total_j_cov_old <- colSums(K(fds_old, "psi5"))
total_j_cov <- data.table(sampleID=names(total_j_cov_old), FRASER=total_j_cov_old)
total_j_cov_jacc <- colSums(K(fds_jacc, "psi5"))
total_j_cov <- merge(total_j_cov, 
                     data.table(sampleID=names(total_j_cov_jacc), Jaccard=total_j_cov_jacc),
                     by="sampleID")
total_j_cov <- melt(total_j_cov, id.vars="sampleID", variable.name="metric", value.name="total_junc_cov")

# merge raw total junction coverage with nr outliers
total_cov_dt_raw <- merge(res_old[, .N, by="sampleID"], total_j_cov_raw, by="sampleID")
setnames(total_cov_dt_raw, "N", "FRASER")
total_cov_dt_raw <- merge(total_cov_dt_raw, res_jacc[, .N, by="sampleID"])
setnames(total_cov_dt_raw, "N", "FRASER2")
total_cov_dt_raw <- melt(total_cov_dt_raw, id.vars=c("sampleID", "total_junc_cov"), value.name="nr_outliers", variable.name="metric")

g_cor_junc_cov <- ggscatter(total_cov_dt_raw, 
                            x="total_junc_cov", y="nr_outliers", 
                            add="reg.line", 
                            conf.int=TRUE, 
                            add.params=list(color="blue", fill="lightgray")) + 
    facet_wrap(~metric) + 
    labs(x="total junction coverage (unfiltered fds)", y="number of outliers") + 
    stat_cor(method="spearman") 
# g_cor_junc_cov
# save plot
ggsave(g_cor_junc_cov, filename=snakemake@output$total_junc_cov,
       width=16, height=8)


# merge total_junction_coverage with nr outliers
total_cov_dt <- merge(res_old[, .N, by="sampleID"], total_j_cov[metric=="FRASER",], by="sampleID")
total_cov_dt <- rbind(total_cov_dt, merge(total_j_cov[metric=="Jaccard",], res_jacc[, .N, by="sampleID"], by="sampleID"))
setnames(total_cov_dt, "N", "nr_outliers")

g_cor_junc_cov_filtered <- ggscatter(total_cov_dt, 
                                    x="total_junc_cov", y="nr_outliers", 
                                    add="reg.line", 
                                    conf.int=TRUE, 
                                    add.params=list(color="blue", fill="lightgray")) + 
    facet_wrap(~metric) + 
    labs(x="total junction coverage (filtered fds)", y="number of outliers") + 
    stat_cor(method="spearman") 
# g_cor_junc_cov_filtered

# save plot
ggsave(g_cor_junc_cov_filtered, filename=snakemake@output$total_junc_cov_filtered,
       width=16, height=8)
