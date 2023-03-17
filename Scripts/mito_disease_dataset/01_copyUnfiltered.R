#'---
#' title: Copy raw unfiltered fds to new folder
#' author: Ines Scheller
#' wb:
#'  log:
#'   - snakemake: '`sm config["log_dir"] + "mito/01_copyUnfiltered.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/"`'
#'  resources:
#'   - mem_mb: 10000
#'  threads: 10
#'  input:
#'   - raw_fds: '`sm config["raw_mito_fds"]`'
#'  output:
#'  - theta:     '`sm config["mito_processed_data"] +
#'                    "/datasets/savedObjects/raw-" + 
#'                    config["mito_dataset_name"] + "/theta.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

raw_input    <- snakemake@input$raw_fds
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# load and copy to new dir
fds <- loadFraserDataSet(file=raw_input)
workingDir(fds) <- workingDir
name(fds) <- paste0("raw-", snakemake@config$dataset_name)
fds <- saveFraserDataSet(fds, rewrite=TRUE)
