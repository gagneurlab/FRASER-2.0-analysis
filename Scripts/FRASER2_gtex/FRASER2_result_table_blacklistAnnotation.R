#'---
#' title: Write FRASER2 results table
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds_fitted: '`sm expand(config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho{rho}_jaccard.h5", rho=config["rhoCutoff"], allow_missing=True)`'
#'     - res_table_junction: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_junction.tsv"`'
#'     - res_table_gene: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv"`'
#'   output:
#'     - blacklist_res_table_junction: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_junction_w_blacklist.tsv"`'
#'     - blacklist_res_table_gene: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene_w_blacklist.tsv"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ read in res tables
res_junc_dt <- fread(snakemake@input$res_table_junction)
res_genes_dt <- fread(snakemake@input$res_table_gene)

#+ read in fds object
fds <- loadFraserDataSet(file=snakemake@input$fds_fitted)

# Annotate results with spliceEventType and blacklist region overlap
# load reference annotation
library(AnnotationDbi)
txdb <- loadDb(snakemake@config$datasets[[snakemake@wildcards$dataset]]$txdb)

#+ annotate the type of splice event and UTR overlap
res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)

#+ set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- snakemake@config$datasets[[snakemake@wildcards$dataset]]$genomeAssembly
if(grepl("grch37", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg19"
}
if(grepl("grch38", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg38"
}

# annotate overlap with blacklist regions
if(assemblyVersion %in% c("hg19", "hg38")){
    res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                        assemblyVersion=assemblyVersion)
    res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, 
                                         assemblyVersion=assemblyVersion)
} else{
    message(date(), ": cannot annotate blacklist regions as no blacklist region\n", 
            "BED file is available for genome assembly version ", assemblyVersion, 
            " as part of FRASER.")
}

#+ save results tables
fwrite(res_junc_dt, file=snakemake@output$blacklist_res_table_junction, sep="\t")
fwrite(res_genes_dt, file=snakemake@output$blacklist_res_table_gene, sep="\t")
