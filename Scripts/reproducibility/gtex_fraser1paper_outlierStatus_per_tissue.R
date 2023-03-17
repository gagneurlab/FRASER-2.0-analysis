#'---
#' title: GTEx Rare variant enrichtment Outlier extraction (per tissue)
#' author: Christian Mertes
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/{dataset}/outlierStatus__{deltaPsi}.Rds"`'   
#'   py:
#'    - |
#'     def get_spot_tissue_clean(wildcards):
#'      t = wildcards.dataset
#'      tissue = t.replace('-_', '')
#'      tissue = tissue.replace('-', '_')
#'      return config["spot_results"] + tissue + "/spot__fullResults.tsv"
#'     def get_leafcutterMD_tissue_clean(wildcards):
#'      t = wildcards.dataset
#'      tissue = t.replace('-_', '')
#'      tissue = tissue.replace('-', '_')
#'      return config["leafcutterMD_results"] + tissue + "/leafcutterMD_testing/results_" + tissue + ".tsv"
#'   threads: 9
#'   resources:
#'     - mem_mb: 64000
#'   input:
#'     - annotation:    '`sm config["gtex_sample_anno"]`'
#'     - leafcutterMD: '`sm get_leafcutterMD_tissue_clean`'
#'     - spot:    '`sm get_spot_tissue_clean`'
#'     - fraser:  '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}__old_filter--gencode34/padjBetaBinomial_theta.h5"`'
#'     - results_fraser:   '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv"`'
#'     - fraser2:  '`sm config["DATADIR"] + "/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho1_jaccard.h5"`'
#'     - results_fraser2:   '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK20_25_minN10/PCA__pc0.1/{dataset}/optQ__newFilt/delta{deltaPsi}/results_gene.tsv"`'
#'   output:
#'     - table:         '`sm config["DATADIR"] + "/GTEx_v8/reproducibility_fraser1paper/{dataset}__outlierStatus__{deltaPsi}.tsv.gz"`'
#'     - wBhtml:        '`sm config["htmlOutputPath"] + "/GTEx_v8/reproducibility_fraser1paper/{dataset}__outlierStatus__{deltaPsi}.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# source config
# source("./src/r/config.R")
source("src/R/variant_enrichment_helper.R")
.libPaths("~/R/4.1/FRASER2/")
library(FRASER)
library(BBmisc)
library(parallel)
library(data.table)
library(tidyr)

curAEVersion <- "PCA"
tissue     <- snakemake@wildcards$dataset
deltaPsi   <- as.numeric(snakemake@wildcards$deltaPsi)
lcmd_file  <- snakemake@input$leafcutterMD
spot_file  <- snakemake@input$spot
fds_file_fraser  <- snakemake@input$fraser
res_file_fraser <- snakemake@input$results_fraser
fds_file_fraser2  <- snakemake@input$fraser2
res_file_fraser2 <- snakemake@input$results_fraser2
anno_file  <- snakemake@input$annotation
threads    <- snakemake@threads

outtable   <- snakemake@output$table

MIN_DELTA_PSI <- deltaPsi
MIN_COVERAGE  <- 10
BPPARAM       <- MulticoreParam(min(4, threads))
register(BPPARAM)

#' 
#' # User input data
#' 
anno_file
lcmd_file
spot_file
res_file_fraser
res_file_fraser2
outtable
BPPARAM

###########################################
#'
#' # Read in outlier calls
#'
###########################################

#+ anno readin
anno <- fread(anno_file)

#'
#' ## Read LeafcutterMD calls
#' 
#+ leafcutterMD readin
tmp_dt <- fread(lcmd_file)
tmp_dt[abs(maxEffect) < MIN_DELTA_PSI, c("pvalue_gene", "padj"):=list(1, 1)]
tmp_dt <- merge(tmp_dt, anno[, .(sample=RNA_ID, subjectID=INDIVIDUAL_ID)], by="sample", all.x=TRUE)
enrich_lcmd <- tmp_dt[,.(subjectID, geneID, tissue=tissue,
                         LeafcutterMD_p=pvalue_gene, LeafcutterMD_fdr=padj, Leafcutter_effect=maxEffect)]
tibble(enrich_lcmd)

###########################################
#'
#' ## Read SPOT calls
#' 
#+ spot readin
tmp_dt <- fread(spot_file)
tmp_dt <- merge(tmp_dt, anno[, .(SAMPLE_ID=RNA_ID, subjectID=INDIVIDUAL_ID)], by="SAMPLE_ID", all.x=TRUE)
enrich_spot <- tmp_dt[,.(subjectID, geneID=GENE_ID, tissue=tissue,
                         SPOT_p=gene_p, SPOT_fdr=gene_fdr)]
# reduce it to uniq rows per subject/gene
enrich_spot <- enrich_spot[!is.na(SPOT_fdr)]
tibble(enrich_spot)


###########################################
#'
#' ## Read FRASER calls
#'
#+ fds readin
loadFds_obj <- FALSE
ncpus <- 1
enrich_fds_obj <- load_fraser_enrichment_data(file=fds_file_fraser, 
                    mc.cores=threads, internCPUs=ncpus, minCoverage=MIN_COVERAGE,
                    minDeltaPsi=MIN_DELTA_PSI, debug=loadFds_obj)
enrich_fds <- enrich_fds_obj[["enrich_obj"]]
# fds_ls     <- enrich_fds_obj[["fds"]]
# res_ls     <- enrich_fds_obj[["res"]]
# names(enrich_fds) <- gsub("__old_filter--gencode34", "", basename(dirname(fds_file_fraser)))
# names(fds_ls) <- gsub("__old_filter--gencode34", "", basename(dirname(fds_file_fraser)))
# names(res_ls) <- gsub("__old_filter--gencode34", "", basename(dirname(fds_file_fraser)))

###########################################
#'
#' ## Read FRASER2 calls
#'
#+ fds fraser2 readin
enrich_fds_obj_fraser2 <- load_fraser2_enrichment_data(file=fds_file_fraser2, 
                    internCPUs=threads, minCoverage=MIN_COVERAGE,
                    minDeltaPsi=MIN_DELTA_PSI, debug=loadFds_obj)
enrich_fds_fraser2 <- enrich_fds_obj_fraser2[["enrich_obj"]]

#'
#' ### use subject ID aka individual identifier instead of run id
#' 
setnames(enrich_fds, "subjectID", "run")
enrich_fds <- merge(enrich_fds, anno[,.(run=RNA_ID, INDIVIDUAL_ID)], by="run")
setnames(enrich_fds, "INDIVIDUAL_ID", "subjectID")
enrich_fds[,run:=NULL]
enrich_fds[,tissue:=tissue]
setcolorder(enrich_fds, unique(c("subjectID", "geneID", "tissue", colnames(enrich_fds))))

setnames(enrich_fds_fraser2, "subjectID", "run")
enrich_fds_fraser2 <- merge(enrich_fds_fraser2, anno[,.(run=RNA_ID, INDIVIDUAL_ID)], by="run")
setnames(enrich_fds_fraser2, "INDIVIDUAL_ID", "subjectID")
enrich_fds_fraser2[,run:=NULL]
enrich_fds_fraser2[,tissue:=tissue]
setcolorder(enrich_fds_fraser2, unique(c("subjectID", "geneID", "tissue", colnames(enrich_fds_fraser2))))


#' 
#' Merge all data sets
#' 
list2merge <- list(lcmd=enrich_lcmd, spot=enrich_spot, fraser=enrich_fds, fraser2=enrich_fds_fraser2)
merged_dt <- Reduce(x=list2merge, f=function(x, y){ 
    merge(x, y, by=c('subjectID', 'geneID', "tissue"), all=TRUE) })


#' 
#' write out table
#' 
fwrite(merged_dt, outtable)
