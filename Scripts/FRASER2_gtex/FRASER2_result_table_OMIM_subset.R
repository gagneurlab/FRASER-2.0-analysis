#'---
#' title: Write FRASER2 results table
#' author: Ines Scheller
#' wb:
#'   threads: 15
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds_fitted: '`sm expand(config["DATADIR"] + "/{dataset_group}/fds/minK{k}_{quantile}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/padjBetaBinomial_rho{rho}_jaccard.h5", rho=config["rhoCutoff"], allow_missing=True)`'
#'     - omim_genes: '`sm config["omim_genes"]`'
#'   output:
#'     - omim_res_table_junction: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_junction_FDRomim.tsv"`'
#'     - omim_res_table_gene: '`sm config["DATADIR"] + "/{dataset_group}/FRASER2_results/minK{k}_{quantile}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene_FDRomim.tsv"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
register(MulticoreParam(snakemake@threads))

#+ set cutoffs
deltaCutoff <- as.numeric(snakemake@wildcards$delta)
fdrCutoff <- as.numeric(snakemake@config$fdrCutoff)
rhoCutoff <- as.numeric(snakemake@config$rhoCutoff)
minCountCutoff <- as.numeric(snakemake@config$minCountCutoff)

#+ read in fds object
fds <- loadFraserDataSet(file=snakemake@input$fds_fitted)

#+ get list of OMIM genes
omim_genes <- fread(snakemake@input$omim_genes)[,unique(gene_v29)]
omim_subset <- rep(list(omim_genes), ncol(fds))
names(omim_subset) <- samples(fds)

#+ get results for FDR on OMIM genes only
res_omim_junc <- results(fds, subsets=list("OMIM"=omim_subset),
                        returnTranscriptomewideResults=FALSE,
                        psiType="jaccard", aggregate=FALSE, 
                        padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                        rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_omim_junc <- as.data.table(res_omim_junc)

#+ compute gene level results table
res_omim_gene <- results(fds, subsets=list("OMIM"=omim_subset),
                    returnTranscriptomewideResults=FALSE,
                    psiType="jaccard", aggregate=TRUE, 
                    padjCutoff=fdrCutoff, deltaPsiCutoff=deltaCutoff, 
                    rhoCutoff=rhoCutoff, minCount=minCountCutoff)
res_omim_gene <- as.data.table(res_omim_gene)


#+ save results tables
fwrite(res_omim_junc, file=snakemake@output$omim_res_table_junction, sep="\t")
fwrite(res_omim_gene, file=snakemake@output$omim_res_table_gene, sep="\t")
