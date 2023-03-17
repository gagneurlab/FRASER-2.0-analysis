#'---
#' title: Filter and clean dataset
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/03_filter_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/"`'
#'  threads: 3
#'  resources:
#'   - mem_mb: 30000
#'  input:
#'   - intron_jaccard: '`sm config["mito_processed_data"] + 
#'                    "/datasets/savedObjects/raw-" + 
#'                    config["mito_dataset_name"] + "/jaccard.h5"`'
#'  output:
#'   - fds: '`sm config["mito_processed_data"] +
#'                "/datasets/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/fds-object.RDS"`'
#'   - done: '`sm config["mito_processed_data"] + 
#'                "/datasets/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/filter.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

# input
dataset    <- snakemake@config$dataset_name
workingDir <- snakemake@params$workingDir

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Apply filter
minExpressionInOneSample <- as.integer(snakemake@wildcards$minK)
quantile <- as.numeric(snakemake@wildcards$quant)
quantileMinExpression <- as.integer(snakemake@wildcards$minN)
minFilterDelta <- snakemake@config$minFilterDelta

# filter introns with low read support and corresponding splice sites
fds <- FRASER:::filterExpression_jaccard(fds,
                        minExpressionInOneSample=minExpressionInOneSample,
                        quantile=(1-quantile),
                        quantileMinExpression=quantileMinExpression,
                        filter=FALSE)

# filter introns that are not variable across samples
fds <- FRASER:::filterVariability_jaccard(fds, minDelta=minFilterDelta, filter=FALSE)
message(date(), ": Filtering done!")

# devNull <- saveFraserDataSet(fds)

# Keep junctions that pass filter
name(fds) <- paste0(dataset, "-minExpr", minExpressionInOneSample, 
                    "-quantile", quantile, "-quantCoverage", quantileMinExpression)
filtered <- mcols(fds, type="j")[,"passed"]
fds <- fds[filtered,]
message(paste("filtered to", nrow(fds), "junctions"))

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$done)