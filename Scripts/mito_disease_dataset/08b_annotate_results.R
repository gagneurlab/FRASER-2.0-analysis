#'---
#' title: Annotate blacklist regions to results
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/08b_annotate_blacklist_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 10
#'  resources:
#'   - mem_mb: 30000
#'  input:
#'   - fdsin: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - resultTableJunc: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                          "/deltaJaccard{delta}/results_per_junction.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableGene: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - blacklist: '`sm config["blacklist_regions_hg19"]`'
#'   - txdb: '`sm config["mito_txdb"]`'
#'  output:
#'   - resultTableJunc_blacklist: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                          "/deltaJaccard{delta}/results_per_junction_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableGene_blacklist: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load libraries
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(AnnotationDbi)
library(data.table)

register(MulticoreParam(snakemake@threads))

# input
fds_file <- snakemake@input$fdsin

# load fds
fds <- loadFraserDataSet(file=fds_file)

# read in junction and gene level results table
res_junc_dt <- fread(snakemake@input$resultTableJunc)
res_gene_dt <- fread(snakemake@input$resultTableGene)

# load reference annotation
txdb <- loadDb(snakemake@input$txdb)

# annotate the type of splice event and UTR overlap
res_junc_dt <- annotateSpliceEventType(result=res_junc_dt, txdb=txdb, fds=fds)
res_gene_dt <- annotateSpliceEventType(result=res_gene_dt, txdb=txdb, fds=fds)

# load blacklist regions
blacklist_input <- snakemake@input$blacklist
# blacklist_gr <- import(blacklist_input, format = "BED")

# annotate overlap with blacklist regions
res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                    blacklist_regions=blacklist_input)
res_gene_dt <- flagBlacklistRegions(result=res_gene_dt, 
                                    blacklist_regions=blacklist_input)


# save annotated results
fwrite(res_junc_dt, file=snakemake@output$resultTableJunc_blacklist)
fwrite(res_gene_dt, file=snakemake@output$resultTableGene_blacklist)
