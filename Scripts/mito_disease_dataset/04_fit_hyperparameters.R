#'---
#' title: Hyper parameter optimization
#' author: Ines Scheller, Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/04_hyper_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/"`'
#'  threads: 30
#'  resources:
#'   - mem_mb: 150000
#'  input:
#'   - filter: '`sm config["mito_processed_data"] + 
#'                "/datasets/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/filter.done" `'
#'  output:
#'   - hyper: '`sm config["mito_processed_data"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/hyper.done" `'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

if ("random_seed" %in% names(snakemake@config)){
    rseed <- snakemake@config$random_seed
    if(isTRUE(rseed)){
        set.seed(42)
    } else if (is.numeric(rseed)){
        set.seed(as.integer(rseed))
    }
}

#+ input
dataset        <- snakemake@config$dataset_name
fds_name       <- basename(dirname(snakemake@input$filter))
workingDir     <- snakemake@params$workingDir
implementation <- snakemake@wildcards$implementation
fds_out_dir    <- file.path(workingDir, implementation)

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

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=fds_name)

# copy to new dir according to used implementation
workingDir(fds) <- fds_out_dir
fds <- saveFraserDataSet(fds)

# Run hyper parameter optimization
# Get range for latent space dimension
mp <- snakemake@config$maxTestedDimensionProportion
a <- 2 
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
    maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- unique(round(exp(seq(log(a),log(b),length.out = Nsteps))))

for(type in c("jaccard")){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type,
                            implementation=implementation,
                            q_param=pars_q,
                            setSubset=30000,
                            plot = FALSE)
    # optData <- data.table(
    #     q=15,
    #     noise=0,
    #     nsubset=50000,
    #     eval=1.7,
    #     aroc=0.5)
    # FRASER:::hyperParams(fds, type=type) <- optData
    
    # fds <- saveFraserDataSet(fds)
}
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$hyper)