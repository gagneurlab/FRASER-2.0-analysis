#'---
#' title: Run FRASER fit for specified q
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 200000
#'   input:
#'     - fds_file: '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/jaccard.h5"`'
#'     - ss_update_done: '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/ss_update.done"`'
#'   output:
#'     - fds_out: '`sm config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/predictedMeans_jaccard.h5"`'
#'   type: script
#'---


#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)

#+ input
fds_file <- snakemake@input$fds_file
implementation <- snakemake@wildcards$implementation # "PCA-BB-Decoder"
txdb <- AnnotationDbi::loadDb(snakemake@config$datasets[[snakemake@wildcards$dataset]]$txdb)
orgdb <- fread(snakemake@config$datasets[[snakemake@wildcards$dataset]]$orgdb)
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ output
out_file <- snakemake@output$fds_out[1]
out_dir <- dirname(dirname(dirname(out_file)))
out_name <- basename(dirname(out_file))

# set pseudocount if requested by implementation
if(grepl("pc", implementation)){
    split_impl <- strsplit(implementation, "__pc", fixed=T)
    implementation <- split_impl[[1]][1]
    pc <- as.numeric(split_impl[[1]][2])
    pseudocount(pc)
}

#+ load fds
fds <- loadFraserDataSet(file=fds_file)
message("Loaded fds.")

#+ set names and dir to new ones
message("Setting new  name of fds to ", out_name)
name(fds) <- as.character(out_name)
message("Setting workingDir of fds to ", out_dir)
workingDir(fds) <- out_dir

#+ filter junctions based on jaccard index
minExpressionInOneSample <- as.integer(snakemake@wildcards$k)
quantile <- 1 - (as.integer(snakemake@wildcards$quantile)/100)
quantileMinExpression <- as.integer(snakemake@wildcards$n)
minFilterDelta <- snakemake@config$minFilterDelta

#+ filter introns with low read support and corresponding splice sites
fds <- FRASER:::filterExpression_jaccard(fds,
                                         minExpressionInOneSample=minExpressionInOneSample,
                                         quantile=quantile,
                                         quantileMinExpression=quantileMinExpression,
                                         filter=FALSE)

#+ filter introns that are not variable across samples
fds <- FRASER:::filterVariability_jaccard(fds, minDelta=minFilterDelta, filter=FALSE)
message(date(), ": Filtering done!")

# Keep junctions that pass filter
filtered <- mcols(fds, type="j")[,"passed"]
fds <- fds[filtered,]
message(paste("filtered to", nrow(fds), "junctions"))

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds)

# Get range for latent space dimension
# mp <- snakemake@config$maxTestedDimensionProportion
mp <- 2
a <- 2 
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp
maxSteps <- 12
if(mp < 6){
    maxSteps <- 15
}
Nsteps <- min(maxSteps, b)
pars_q <- unique(round(exp(seq(log(a),log(b),length.out = Nsteps))))

#+ run FRASER fit with hyper param opt for q
for(pt in c("jaccard")){
    message(date(), ": starting hyper param opt for type ", pt, " ...")
    fds <- optimHyperParams(fds, type=pt,
                            implementation=implementation,
                            q_param=pars_q,
                            setSubset=30000,
                            plot = FALSE)
    gc()
    message(date(), ": starting ", implementation, " fit for type ", pt, " ...")
    fds <- fit(fds, type=pt, q=bestQ(fds, type=pt), 
               implementation=implementation, rhoRange=c(-30, 30))
    gc()
}

#+ write output
fds <- saveFraserDataSet(fds)

