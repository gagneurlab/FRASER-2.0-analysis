#'---
#' title: Nr of outliers per sample (rare disease cohorts)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/mito/nrOutliers_rareDiseaseCohorts.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - prokisch_fdr_comparison: '`sm config["mito_processed_results"] + "/ggplot_rds/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/FDR_subset_MAF0.001_no_utr_ggplots_F1_old_filter.Rds"`'
#'     - prokisch_omim_res_table: '`sm config["mito_processed_results"] + "/results/PCA__pc0.1/gencode34/fib-minExpr20-quantile0.25-quantCoverage10/deltaJaccard0.1/results_gene_FDRomim.tsv"`'
#'     - udn_nr_outliers_comparison: '`sm config["DATADIR"] + "/udn/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/nrOutliers_comparison_ggplot.Rds"`'
#'   output:
#'    - nr_out_combined: '`sm config["DATADIR"] + "/mito/nrOutliers_rareDiseaseCohorts.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

library(data.table)

#+ read in data for each cohort
fdr_comparison_plots <- readRDS(snakemake@input$prokisch_fdr_comparison)
udn_plots <- readRDS(snakemake@input$udn_nr_outliers_comparison)

#+ mito (yepez et al) dataset
mito_dt <- fdr_comparison_plots[["boxplots_numOut_gene"]]$data
mito_dt[, dataset_group:="Yépez et al. dataset"] # "Mitochondrial disease cohort"]
# mito_dt[, dataset_label:=paste0("Fibroblasts (N=", uniqueN(sampleID),")")]
mito_dt[, dataset_label:="Fibroblasts"]
mito_dt

#+ mito (yepez et al) dataset (OMIM genes only)
mito_dt_omim <- fread(snakemake@input$prokisch_omim_res_table)
mito_dt_omim <- mito_dt_omim[!duplicated(mito_dt_omim, by=c("sampleID","seqnames","start","end","strand")), ]
# get nr of outliers per sample (gene-level)
mito_dt_omim[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
mito_dt_omim[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
# add name of dataset
mito_dt_omim[, dataset_group:="Yépez et al. dataset"] # "Mitochondrial disease cohort"]
# mito_dt_omim[, dataset_label:=paste0("Fibroblasts (N=", uniqueN(mito_dt$sampleID),")")]
mito_dt_omim[, dataset_label:="Fibroblasts"]
mito_dt_omim[, method:="FRASER2 (OMIM)"]
mito_dt_omim <- mito_dt_omim[, .(num_out = unique(numOutliersPerSample)), by="sampleID,method,dataset_label,dataset_group"]
mito_dt_omim

#+ all mito scenarios
mito_dt <- rbind(mito_dt, mito_dt_omim)
dcast_dt <- dcast(mito_dt, sampleID + dataset_label + dataset_group ~ method, value.var="num_out")
dcast_dt[is.na(`FRASER2 (OMIM)`), `FRASER2 (OMIM)` := 0]
mito_dt <- melt(dcast_dt, id.vars=c("sampleID", "dataset_label","dataset_group"), value.name="num_out", variable.name="method")
mito_dt

#+ udn
udn_dt <- udn_plots[["boxplots_numOut_gene"]]$data
udn_dt[, dataset_group := "Undiagnosed Disease Network (UDN)"]
udn_dt[, dataset_label := NULL]
setnames(udn_dt, "dataset_name", "dataset_label")

#+ combine the different cohorts
num_out_dt <- rbind(mito_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)], 
                    udn_dt[,.(sampleID, dataset_label, dataset_group, method, num_out)])
tissue_order <- num_out_dt[dataset_group == "Undiagnosed Disease Network (UDN)", median(num_out),by="dataset_label,method"][method == "FRASER"][order(-V1), dataset_label]
num_out_dt[, dataset_label := factor(dataset_label, levels=tissue_order)]
num_out_dt[method == "FRASER2", method := "FRASER 2.0"]
num_out_dt[method == "FRASER2 (OMIM)", method := "FRASER 2.0\n(OMIM)"]
num_out_dt[method == "FRASER2\n(OMIM + RV)", method := "FRASER 2.0\n(OMIM + RV)"]
num_out_dt[, method := factor(method, levels=c("FRASER", "FRASER 2.0", "FRASER 2.0\n(OMIM)", "FRASER 2.0\n(OMIM + RV)"))]
num_out_dt[, dataset_group := factor(dataset_group, levels=c("Yépez et al. dataset", "Undiagnosed Disease Network (UDN)"))]

#+ save
fwrite(num_out_dt, file=snakemake@output$nr_out_combined, sep="\t")