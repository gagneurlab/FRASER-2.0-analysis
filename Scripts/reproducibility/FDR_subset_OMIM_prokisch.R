#'---
#' title: Write FRASER2 results table
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/FDR_subset_OMIM_results_delta{delta}.Rds"`'
#'   threads: 15
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds_fitted: '`sm config["DATADIR"] + "/mito/processed_results/datasets/PCA__pc0.1/savedObjects/fib-minExpr20-quantile0.25-quantCoverage10--gencode34__FDR_sub_MAF0.001_no_utr/fds-object.RDS"`'
#'     - omim_genes: '`sm config["omim_genes"]`'
#'   output:
#'     - omim_res_table_junction: '`sm config["DATADIR"] + "/mito/processed_results/results/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/deltaJaccard{delta}/results_junction_FDRomim.tsv"`'
#'     - omim_res_table_gene: '`sm config["DATADIR"] + "/mito/processed_results/results/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/deltaJaccard{delta}/results_gene_FDRomim.tsv"`'
#'   type: script
#'---
saveRDS(snakemake, snakemake@log$snakemake)

#+ load FRASER
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ set cutoffs
deltaCutoff <- as.numeric(snakemake@wildcards$delta)
fdrCutoff <- as.numeric(snakemake@config$fdrCutoff)
rhoCutoff <- NA # 0.1
minCountCutoff <- as.numeric(snakemake@config$minCountCutoff)

#+ read in fds object
fds <- loadFraserDataSet(file=snakemake@input$fds_fitted)

#+ get list of OMIM genes
omim_genes <- fread(snakemake@input$omim_genes)[,unique(gene_v29)]
omim_subset <- rep(list(omim_genes), ncol(fds))
names(omim_subset) <- samples(fds)

#+ get list of OMIM genes with rare variants
omim_plus_rv_genes <- metadata(fds)$genes_rare_omim_var

# calc padj values on subsets
fds <- calculatePadjValues(fds, type="jaccard", 
                           subsets=list("OMIM"=omim_subset,
                                        "OMIM + RV"=omim_plus_rv_genes) 
)

#+ get results for FDR on OMIM genes only
res_omim_junc <- results(fds, 
                         returnTranscriptomewideResults=FALSE,
                         psiType="jaccard", aggregate=FALSE, 
                         padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                         rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_omim_junc <- as.data.table(res_omim_junc)

#+ compute gene level results table
res_omim_gene <- results(fds, 
                         returnTranscriptomewideResults=FALSE,
                         psiType="jaccard", aggregate=TRUE, 
                         padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                         rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_omim_gene <- as.data.table(res_omim_gene)


#+ save results tables
fwrite(res_omim_junc, file=snakemake@output$omim_res_table_junction, sep="\t")
fwrite(res_omim_gene, file=snakemake@output$omim_res_table_gene, sep="\t")
