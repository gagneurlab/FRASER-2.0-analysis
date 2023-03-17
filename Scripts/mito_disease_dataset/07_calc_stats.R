#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/07_stats_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_results"] + "/datasets/{implementation}/"`'
#'  threads: 20
#'  resources:
#'   - mem_mb: 75000
#'  input:
#'   - fdsin:  '`sm config["mito_processed_results"] + 
#'                 "/datasets/{implementation}/savedObjects/" + 
#'                 config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                 "/fds-object.RDS"`'
#'  output:
#'   - fdsout: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                  "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - genePvals: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                "/genePvals.done"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

dataset    <- snakemake@config$dataset_name
fds_name   <- basename(dirname(snakemake@input$fdsin))
annotation <- snakemake@config$annotation
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load Zscores data
fds <- loadFraserDataSet(file=fdsFile)

# Calculate stats
for (type in c("jaccard")) {
    # Zscores
    fds <- calculateZscore(fds, type=type)
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type)
}

fds <- saveFraserDataSet(fds)
file.create(snakemake@output$genePvals)
