#'---
#' title: Results of FRASER analysis
#' author: Ines Scheller, Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/06_annotateResults_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/{implementation}/"`'
#'   - outputDir: '`sm config["mito_processed_results"] + "/datasets/{implementation}/"`'
#'  threads: 20
#'  resources:
#'   - mem_mb: 50000
#'  input:
#'   - fdsin:  '`sm config["mito_processed_data"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" +
#'                "/predictedMeans_jaccard.h5"`'
#'   - txdb: '`sm config["mito_txdb"]`'
#'   - gene_name_mapping: '`sm config["mito_gene_name_mapping"]`'
#'  output:
#'   - annotated_fds: '`sm config["mito_processed_results"] + 
#'                 "/datasets/{implementation}/savedObjects/" + 
#'                 config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["mito_annotation"] + 
#'                 "/fds-object.RDS"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)
library(AnnotationDbi)

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
fds_input <- loadFraserDataSet(dir=workingDir, name=fds_name)

# Read annotation and match the chr style
txdb <- loadDb(snakemake@input$txdb)
orgdb <- fread(snakemake@input$gene_name_mapping)

seqlevels_fds <- seqlevelsStyle(fds_input)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds

# Annotate the fds with gene names and save it as a new object
fds_input <- annotateRangesWithTxDb(fds_input, txdb = txdb, orgDb = orgdb, 
                                    feature = 'gene_name', featureName = 'hgnc_symbol', keytype = 'gene_id')

# add basic annotations for overlap with the reference annotation
# run this function before creating the results table
fds_input <- annotateIntronReferenceOverlap(fds_input, txdb)

# save fds
fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(fds_name, annotation, sep = '--'), rewrite = TRUE)
