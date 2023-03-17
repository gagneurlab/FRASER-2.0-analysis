#'---
#' title: Results of FRASER analysis
#' author: Ines Scheller, Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/08_results_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_results"] + "/datasets/{implementation}/"`'
#'   - padjCutoff: '`sm config["fdrCutoff"]`'
#'   - zScoreCutoff: '`sm config["zScoreCutoff"]`'
#'   - deltaJaccardCutoff: '`sm config["deltaCutoff"]`'
#'  threads: 10
#'  resources:
#'   - mem_mb: 30000
#'  input:
#'   - fdsin: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - genePvals: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                "/genePvals.done"`'
#'  output:
#'   - resultTableJunc: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                          "/deltaJaccard{delta}/results_per_junction.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableGene: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)
library(AnnotationDbi)
library(readr)

annotation    <- snakemake@config$annotation
dataset    <- snakemake@config$dataset_name
fds_name   <- basename(dirname(snakemake@input$fdsin))
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
outputDir  <- snakemake@params$outputDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=fds_name) # paste(fds_name, annotation, sep="--"))

# Extract results per junction
res_junc <- results(fds, psiType="jaccard",
                    padjCutoff=snakemake@params$padjCutoff,
                    zScoreCutoff=snakemake@params$zScoreCutoff,
                    deltaPsiCutoff=snakemake@params$deltaJaccardCutoff)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')

# Add features
if(nrow(res_junc_dt) > 0){
    
    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
    
    # add colData to the results
    res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
    res_junc_dt[, c("bamFile", "pairedEnd") := NULL]
    setnames(res_junc_dt, 'STRAND', 'STRAND_SPECIFIC')  # otherwise it's confusing with the 'strand' column from the junction
} else{
    warning("The aberrant splicing pipeline gave 0 results for the ", fds_name, " fds_name")
}

# Aggregate results by gene
if(length(res_junc) > 0){
    # Extract results per gene
    res_gene <- results(fds, psiType="jaccard",
                        aggregate=TRUE, collapse=FALSE,
                        padjCutoff=snakemake@params$padjCutoff,
                        zScoreCutoff=snakemake@params$zScoreCutoff,
                        deltaPsiCutoff=snakemake@params$deltaJaccardCutoff)
    res_genes_dt   <- as.data.table(res_gene)
    print('Results per gene extracted')
    res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
    res_genes_dt[, c("bamFile", "pairedEnd") := NULL]
    setnames(res_genes_dt, 'STRAND', 'STRAND_SPECIFIC')
    
} else res_genes_dt <- data.table()

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)