#'---
#' title: Time FRASER 2.0 run (fit + hyper param search, pval, results table)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/runtimes/{dataset_group}/fraser_v{version}/{implementation}/runtimes_{dataset}.Rds"`'
#'   py:
#'   - |
#'    def get_mem_mb(wildcards, attempt):
#'     return attempt * 200000
#'   threads: 15
#'   resources:
#'     - mem_mb: '`sm get_mem_mb`'
#'   input:
#'     - fds_file:           '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/jaccard.h5"`'
#'     - ss_update_done:     '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/ss_update.done"`'
#'   output:
#'     - fds_out:            '`sm config["DATADIR"] + "/{dataset_group}/runtimes/fraser_v{version}/{implementation}/fds/savedObjects/{dataset}/fds-object.RDS"`'
#'     - res_table_junction: '`sm config["DATADIR"] + "/{dataset_group}/runtimes/fraser_v{version}/{implementation}/results/{dataset}/results_junction.tsv"`'
#'     - res_table_gene:     '`sm config["DATADIR"] + "/{dataset_group}/runtimes/fraser_v{version}/{implementation}/results/{dataset}/results_gene.tsv"`'
#'     - timings:            '`sm config["DATADIR"] + "/{dataset_group}/runtimes/fraser_v{version}/{implementation}/timings/{dataset}_runtime.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load FRASER
# .libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)

#+ input
fds_file <- snakemake@input$fds_file
fraser_version <- snakemake@wildcards$version
implementation <- snakemake@wildcards$implementation # "PCA__pc0.1"
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
} else{
    pseudocount(1)
}

#+ prepare object to store runtime timings
timings <- list()
timings[["method"]] <- switch(fraser_version,
                              "2" = "FRASER 2.0 (recommmended filtering settings for FRASER 2.0)", 
                              "2_oldFilt" = "FRASER 2.0 (filtering settings as in FRASER)", 
                              "1" = "FRASER")
timings[["tissue"]] <- snakemake@wildcards$dataset
timings[["begin"]] <- Sys.time()

#+ load fds
timings[["begin_load_fds"]] <- Sys.time()
fds <- loadFraserDataSet(file=fds_file)
timings[["end_load_fds"]] <- Sys.time()
message("Loaded fds.")

timings[["sample_size"]] <- ncol(fds)
timings[["nr_introns"]] <- nrow(fds)

#+ set names and dir to new ones
message("Setting new  name of fds to ", out_name)
name(fds) <- as.character(out_name)
message("Setting workingDir of fds to ", out_dir)
workingDir(fds) <- out_dir

#+ set splice metrics to fit
if(fraser_version == "2" | fraser_version == "2_oldFilt"){
    fit_metrics <- "jaccard" 
} else{
    fit_metrics <- c("psi5", "psi3", "theta")
    fds <- calculatePSIValues(fds, types=fit_metrics)
}
fitMetrics(fds) <- fit_metrics

#+ filter junctions based on jaccard index
minExpressionInOneSample <- 20
minFilterDelta <- 0
if(fraser_version == "2"){
    quantile <- 0.75
    quantileMinExpression <- 10
    filterJaccard <- TRUE
} else if(fraser_version == "2_oldFilt"){
    quantile <- 0.05
    quantileMinExpression <- 1
    filterJaccard <- TRUE
} else{
    quantile <- 0.05
    quantileMinExpression <- 1
    filterJaccard <- FALSE
}
    
#+ filter introns with low read support and corresponding splice sites
timings[["begin_filtering"]] <- Sys.time()
fds <- filterExpressionAndVariability(fds, filterOnJaccard=filterJaccard,
                                         minExpressionInOneSample=minExpressionInOneSample,
                                         quantile=quantile,
                                         quantileMinExpression=quantileMinExpression,
                                         minDeltaPsi=minFilterDelta,
                                         filter=FALSE)
timings[["end_filter_calc"]] <- Sys.time()
message(date(), ": Filtering done!")

# Keep junctions that pass filter
filtered <- mcols(fds, type="j")[,"passed"]
fds <- fds[filtered,]
timings[["end_filtering"]] <- Sys.time()
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
timings[["begin_hyperpar"]] <- Sys.time()
for(pt in fit_metrics){
    message(date(), ": starting hyper param opt for type ", pt, " ...")
    fds <- optimHyperParams(fds, type=pt,
                            implementation=implementation,
                            q_param=pars_q,
                            setSubset=30000,
                            plot = FALSE)
}
timings[["end_hyperpar"]] <- Sys.time()
gc()

timings[["begin_fit"]] <- Sys.time()
for(pt in fit_metrics){
    message(date(), ": starting ", implementation, " fit for type ", pt, " ...")
    fds <- fit(fds, type=pt, q=bestQ(fds, type=pt), 
               implementation=implementation, rhoRange=c(-30, 30))
}
timings[["end_fit"]] <- Sys.time()
gc()

#+ annotate genes
seqlevels_fds <- seqlevelsStyle(fds)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgdb, 
                              # filter=list('gene_type'='protein_coding'),
                              feature = 'gene_name', 
                              featureName = 'hgnc_symbol', 
                              keytype = 'gene_id')
timings[["nr_genes"]] <- length(FRASER:::getGeneIDs(fds, type="jaccard", unique=TRUE))

timings[["begin_pvals"]] <- Sys.time()
for(pt in fit_metrics){
    message(date(), ": starting BB pval calculation for type ", pt, " ...")
    fds <- calculatePvalues(fds, type=pt, implementation=implementation)
    fds <- calculatePadjValues(fds, type=pt, rhoCutoff=NA)
}
timings[["end_pvals"]] <-  Sys.time()
gc()

#+ save fds
fds <- saveFraserDataSet(fds)

#+ set delta cutoff
deltaCutoff <- 0.1
fdrCutoff <- as.numeric(snakemake@config$fdrCutoff)
minCountCutoff <- as.numeric(snakemake@config$minCountCutoff)

#+ compute junction level results table
timings[["begin_results"]] <-  Sys.time()
timings[["begin_res_intronlevel"]] <-  Sys.time()
res_junction <- results(fds, psiType=fit_metrics, aggregate=FALSE, 
                        padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                        rhoCutoff=NA, minCount=minCountCutoff)
res_junction <- as.data.table(res_junction)
timings[["end_res_intronlevel"]] <- Sys.time()

#+ compute gene level results table
timings[["begin_res_genelevel"]] <- Sys.time()
res_gene <- results(fds, psiType=fit_metrics, aggregate=TRUE, 
                    padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                    rhoCutoff=NA, minCount=minCountCutoff)
res_gene <- as.data.table(res_gene)
timings[["end_res_genelevel"]] <- Sys.time()
timings[["end_results"]] <-  Sys.time()

#+ save results tables
fwrite(res_junction, file=snakemake@output$res_table_junction, sep="\t")
fwrite(res_gene, file=snakemake@output$res_table_gene, sep="\t")

#+ get last time point and save
timings[["end"]] <- Sys.time()

timings_outfile <- snakemake@output$timings
if(!dir.exists(dirname(timings_outfile))){
    dir.create(dirname(timings_outfile), recursive=TRUE)
}
saveRDS(timings, file=timings_outfile)


