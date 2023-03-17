#'---
#' title: Number of outliers for FRASER2
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/{dataset}/minK{k}_{q}_minN{n}/{implementation}/delta{delta}/FRASER2_nrOutliers.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - fds_fraser2:  '`sm config["DATADIR"] + "/GTEx_v8/fds/minK{k}_{q}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/fds-object.RDS"`'
#'     - res_fraser2: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_{q}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv"`'
#'     - res_fraser2_junction: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_{q}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_junction.tsv"`'
#'   output:
#'     - nrOutliers_table: '`sm config["DATADIR"] + "/GTEx_v8/{dataset}/minK{k}_{q}_minN{n}/{implementation}/optQ/delta{delta}/FRASER2_nrOutliers.tsv"`'
#'   type: script
#'---


saveRDS(snakemake, snakemake@log$snakemake)


.libPaths("~/R/4.1/FRASER2_BB_loss")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

tissue <- snakemake@wildcards$dataset
implementation <- snakemake@wildcards$implementation
k <- snakemake@wildcards$k
q <- snakemake@wildcards$q
n <- snakemake@wildcards$n

message(date(), ": loading fds object for tissue = ", tissue, " ...")
fds <- loadFraserDataSet(file=snakemake@input$fds_fraser2)
# ab <- aberrant(fds, type="jaccard", by="sample", aggregate=TRUE)
# ab_junc <- aberrant(fds, type="jaccard", by="sample", aggregate=FALSE)
res_junc <- fread(snakemake@input$res_fraser2_junction)
res_junc[,tmp:=paste(seqnames, start, end, strand, sep="_")]
ab_junc <- res_junc[, uniqueN(tmp), by="sampleID"]
ab_junc <- merge(data.table(sampleID=samples(fds)), ab_junc, by=c("sampleID"), all.x=TRUE, sort=FALSE)
ab_junc[is.na(V1), V1:=0]
res_gene <- fread(snakemake@input$res_fraser2)
ab <- res_gene[, uniqueN(hgncSymbol), by="sampleID"]
ab <- merge(data.table(sampleID=samples(fds)), ab, by="sampleID", all.x=TRUE, sort=FALSE)
ab[is.na(V1), V1:=0]
dt <- data.table(sampleID=ab$sampleID, FRASER2_nrOutlierGenes=ab$V1, FRASER2_nrOutlierJunctions_jaccard=ab_junc$V1)
dt[,tissue:=tissue]
dt[,implementation:=implementation]
dt[,k:=k]
dt[,q:=q]
dt[,n:=n]

#+ save table
fwrite(dt, file=snakemake@output$nrOutliers_table)
