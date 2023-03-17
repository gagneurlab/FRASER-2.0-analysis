#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes
#' wb:
#'   py:
#'   - |
#'    def get_input_fds(wildcards):
#'     return config["datasets"][wildcards.dataset]["fds_file"]
#'   threads: 15
#'   resources:
#'     - mem_mb: 100000
#'   params:
#'    - workingDir: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_splicing/datasets/"`'
#'    - outputDir: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/datasets/"`'
#'    - padjCutoff: '`sm 0.1`'
#'    - zScoreCutoff: '`sm 0`'
#'    - deltaPsiCutoff: '`sm 0.3`'
#'   input:
#'     - fdsin: '`sm get_input_fds`'
#'     - txdb: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode34/txdb.db"`'
#'     - gene_name_mapping: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode34/gene_name_mapping_gencode34.tsv"`'
#'   output:
#'    - resultTableJunc: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results" +
#'                          "/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results_per_junction.tsv"`'
#'    - resultTableGene: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results" +
#'                          "/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv"`'
#'    - fds: '`sm "/s/project/gtex_genetic_diagnosis/v8/processed_results" +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}_old_filter--gencode34/fds-object.RDS"`'
#'   type: script
#'---

library(FRASER)
library(data.table)
library(AnnotationDbi)
library(BiocParallel)
library(magrittr)

write_tsv <- function(x, file, row.names = FALSE, ...){
    write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}

annotation    <- "gencode34"
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
outputDir  <- snakemake@params$outputDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
DelayedArray::setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds_input <- loadFraserDataSet(dir=workingDir, name=dataset)

# Read annotation and match the chr style
txdb <- loadDb(snakemake@input$txdb)
orgdb <- fread(snakemake@input$gene_name_mapping)

seqlevels_fds <- seqlevelsStyle(fds_input)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds

# Annotate the fds with gene names and save it as a new object
fds_input <- annotateRangesWithTxDb(fds_input, txdb = txdb, orgDb = orgdb, 
                                    feature = 'gene_name', featureName = 'hgnc_symbol', keytype = 'gene_id')

fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(paste0(dataset, "_old_filter"), annotation, sep = '--'), rewrite = TRUE)

# Extract results per junction
res_junc <- results(fds,
                    padjCutoff=snakemake@params$padjCutoff,
                    zScoreCutoff=snakemake@params$zScoreCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')
saveFraserDataSet(fds, dir=outputDir)


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
    warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}

# Aggregate results by gene
if(length(res_junc) > 0){
    res_genes_dt <- resultsByGenes(res_junc) %>% as.data.table
    res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
    res_genes_dt[, c("bamFile", "pairedEnd") := NULL]
    setnames(res_genes_dt, 'STRAND', 'STRAND_SPECIFIC')
    
} else res_genes_dt <- data.table()

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)
