#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/05_fit_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/{implementation}/"`'
#'  threads: 20
#'  resources:
#'   - mem_mb: 125000
#'  input:
#'   - hyper: '`sm config["mito_processed_data"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/hyper.done" `'
#'  output:
#'   - fdsout: '`sm config["mito_processed_data"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/predictedMeans_jaccard.h5"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

dataset    <- snakemake@config$dataset_name
fds_name   <- basename(dirname(snakemake@input$hyper))
workingDir <- snakemake@params$workingDir
implementation <- snakemake@wildcards$implementation

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# set pseudocount if requested by implementation
if(grepl("pc", implementation)){
    split_impl <- strsplit(implementation, "__pc", fixed=T)
    implementation <- split_impl[[1]][1]
    pc <- as.numeric(split_impl[[1]][2])
    pseudocount(pc)
}

fds <- loadFraserDataSet(dir=workingDir, name=fds_name)

# Fit autoencoder
# run it for every type
for(type in c("jaccard")){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
    fds <- saveFraserDataSet(fds)
}