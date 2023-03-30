#'---
#' title: Number of outliers for FRASER1
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/{dataset}/minK20_95_minN1/PCA/FRASER1_nrOutliers.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - fds_fraser: '`sm config["general_data_dir"] + "/gtex_genetic_diagnosis/v8/processed_data/aberrant_splicing/datasets/savedObjects/{dataset}_old_filter/fds-object.RDS"`'
#'     - res_fraser: '`sm  config["general_data_dir"] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv"`'
#'     - res_fraser_junction: '`sm  config["general_data_dir"] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results_per_junction.tsv"`'
#'   output:
#'     - nrOutliers_table: '`sm config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_95_minN1/PCA/optQ/FRASER1_nrOutliers.tsv"`'
#'   type: script
#'---



saveRDS(snakemake, snakemake@log$snakemake)

.libPaths("~/R/4.1/FRASER2_BB_loss")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

tissue <- snakemake@wildcards$dataset

# FRASER1: filtering quantile: 95, n:1
fds <- loadFraserDataSet(file=snakemake@input$fds_fraser)
res_junc <- fread(snakemake@input$res_fraser_junction)
res_junc[,tmp:=paste(seqnames, start, end, strand, sep="_")]
ab_junc <- res_junc[, uniqueN(tmp), by="sampleID,type"]
ab_junc <- merge(data.table(sampleID=rep(samples(fds), each=3), type=rep(c("psi5", "psi3", "theta"), ncol(fds))), ab_junc, by=c("sampleID", "type"), all.x=TRUE, sort=FALSE)
ab_junc[is.na(V1), V1:=0]
ab_junc <- dcast(ab_junc, formula="sampleID ~ type", value.var="V1")
res_gene <- fread(snakemake@input$res_fraser)
ab_old <- res_gene[, uniqueN(hgncSymbol), by="sampleID"]
ab_old <- merge(data.table(sampleID=samples(fds)), ab_old, by="sampleID", all.x=TRUE, sort=FALSE)
ab_old[is.na(V1), V1:=0]
dt <- data.table(sampleID=ab_junc$sampleID, FRASER_nrOutlierJunctions_psi5=ab_junc$psi5, 
                 FRASER_nrOutlierJunctions_psi3=ab_junc$psi3, FRASER_nrOutlierJunctions_theta=ab_junc$theta, 
                 FRASER_nrOutlierGenes=ab_old$V1)
dt[,tissue:=tissue]
dt[,implementation:="PCA"]
dt[,k:=20]
dt[,q:=95]
dt[,n:=1]

#+ save table
fwrite(dt, file=snakemake@output$nrOutliers_table)
