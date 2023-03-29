#'---
#' title: Paper numbers
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/paper_numbers.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 32000
#'   input:
#'     - filtering_info: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/filtering_info.tsv"`'
#'     - variant_recall_all_tissues: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/plot_rds/FRASER2_vs_others_allTissues_rv_recall_plots_{varType}.Rds", varType=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'     - variant_enrich_fdr_all: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/all_final_prec_rec_ggplots.Rds"`'
#'     - variant_enrich_all: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/all_final_rv_recall_ggplots.Rds"`'
#'     - variant_pr_fdr_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_rv_recall_plots_{snptype}.Rds", dataset=config["tissues_for_reproducibility"], snptype=["rareSplicing", "rareSpliceAI", "rareMMSplice", "rareAbSplice"])`'
#'     - nr_outliers_per_sample_gtex:  '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta0.1/nrOutliers_comparison_ggplot.Rds"`'
#'     - nr_outliers_per_sample_rd_cohorts: '`sm config["DATADIR"] + "/mito/nrOutliers_rareDiseaseCohorts.tsv"`'
#'     - reproducibility_table:  '`sm config["DATADIR"] + "/GTEx_v8/reproducibility_fraser1paper/rareSplicing__0.0_reproducibility.tsv.gz"`'
#'     - fraser2_FDR_sub: '`sm config["mito_processed_results"] + "/datasets/PCA__pc0.1/savedObjects/" + 
#'                "fib-minExpr20-quantile0.25-quantCoverage10--gencode34__FDR_sub_MAF0.001_no_utr" +
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'     - results_mito_tw: '`sm config["mito_processed_results"] + 
#'                          "/results/PCA__pc0.1/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr20-quantile0.25-quantCoverage10" + 
#'                          "/deltaJaccard0.1/results_blacklist.tsv"`'
#'     - patho_sample_anno: '`sm config["mito_sample_anno"]`'
#'     - power_analysis_full_res_patho: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_full.tsv"`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/paper/paper_numbers.html"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load packages
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)

#' 
#' # GTEx analysis
#' 

filter_info_dt <- fread(snakemake@input$filtering_info)
#' 
#' ### Nr of genes and introns for each filter setting across `r filter_info_dt[, uniqueN(tissue)]` tissues 
#' 
out_dt <- filter_info_dt[, .(median_nr_genes = median(ngenes), 
                          median_nr_introns = median(njuncs)),
                      by="k,q,n"]
DT::datatable(
    out_dt,
    caption = 'Nr of genes and introns for each filter setting',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)
#' 
#' ### Precision and recall on splice-disrupting candidate variants at FDR cutoff 0.1
#' 
# pr_dt_tissue <- readRDS(snakemake@input$variant_pr_fdr_all)$precision_recall$layers[[2]]$data
pr_rec_data_tissue <- rbindlist(lapply(snakemake@input$variant_pr_fdr_all, function(infile){
    var_recall_plots <- readRDS(infile)
    pr_rec_plot <- var_recall_plots$precision_recall
    pr_rec_data <- pr_rec_plot$layers[[2]]$data[Method %in% c("FRASER", "FRASER2")]
    vartype <- strsplit(gsub(".Rds", "", basename(infile)), "_", fixed=TRUE)[[1]][8]
    pr_rec_data[, snptype := vartype]
    dataset <- basename(dirname(dirname(dirname(infile))))
    pr_rec_data[, tissue := dataset]
    return(pr_rec_data)
}))
comp_pr_dt_tissue <- dcast(melt(pr_rec_data_tissue, id.vars=c("Method", "snptype", "tissue", "Type", "Cutoff"), 
                                measure.vars=c("precision", "recall"), 
                                variable.name="eval_metric"), 
                           value.var="value", snptype + tissue + Type + Cutoff + eval_metric ~ Method)
comp_pr_dt_tissue[, fc_FRASER2_vs_FRASER := FRASER2 / FRASER]
comp_pr_dt_tissue
DT::datatable(
    comp_pr_dt_tissue[, .(median_fc = median(fc_FRASER2_vs_FRASER),
                          median_FRASER = median(FRASER),
                          median_FRASER2 = median(FRASER2)), 
                      by="snptype,eval_metric"],
    caption = 'Comparison of Precision and Recall of FRASER2 and FRASER on splice-disrupting variants at FDR cutoff of 0.1 across tissues',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)


#' 
#' ### Precision and recall of splice-disrupting candidate variants (jointly across tissues)
#' 
pr_rec_data <- rbindlist(lapply(snakemake@input$variant_recall_all_tissues, function(infile){
    var_recall_plots <- readRDS(infile)
    pr_rec_plot <- var_recall_plots$precision_recall
    pr_rec_data <- pr_rec_plot$layers[[2]]$data[Method %in% c("FRASER", "FRASER2")]
    vartype <- strsplit(gsub(".Rds", "", basename(infile)), "_", fixed=TRUE)[[1]][8]
    pr_rec_data[, snptype := vartype]
    return(pr_rec_data)
}))
comp_dt <- dcast(melt(pr_rec_data, id.vars=c("Method", "Cutoff", "snptype"), measure.vars=c("recall", "precision"), variable.name="eval_metric"), 
                 value.var="value", snptype + Cutoff ~ Method + eval_metric)
comp_dt[, fc_precision := FRASER2_precision / FRASER_precision]
comp_dt[, fc_recall := FRASER2_recall / FRASER_recall]
DT::datatable(
    comp_dt,
    caption = 'Comparison of FRASER2 to FRASER on splice-disrupting variants jointly on all tissues',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#' 
#' ### Recall of splice-disrupting candidate variants at rank of 10/20 outliers per sample
#' 
var_enrich_all <- readRDS(snakemake@input$variant_enrich_all)
recall_dt_tissue <- var_enrich_all$recall_across_tissues$data[Method %in% c("FRASER", "FRASER2")]
comp_dt_tissue <- dcast(recall_dt_tissue, value.var="recall", snptype + mean_outliers_per_sample + tissue ~ Method)
comp_dt_tissue[, fc_recall := FRASER2 / FRASER]
comp_dt_tissue
DT::datatable(
    comp_dt_tissue[, .(median_fc_recall = median(fc_recall),
                       median_recall_FRASER = median(FRASER),
                       median_recall_FRASER2 = median(FRASER2)), 
                   by="snptype,mean_outliers_per_sample"],
    caption = 'Comparison of FRASER2 to FRASER on splice-disrupting variants on median fc across tissues',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#' 
#' ### AUPRC on splice-disrupting candidate variants at FDR cutoff 0.1
#' 
var_enrich_fdr_all <- readRDS(snakemake@input$variant_enrich_fdr_all)
auprc_dt_tissue <- var_enrich_fdr_all$auprc_across_tissues$data[Method %in% c("FRASER", "FRASER2")]
comp_auprc_dt_tissue <- dcast(auprc_dt_tissue, value.var="AUPRC", snptype + tissue ~ Method)
comp_auprc_dt_tissue[, fc_auprc := FRASER2 / FRASER]
comp_auprc_dt_tissue
DT::datatable(
    comp_auprc_dt_tissue[, .(median_fc_auprc = median(fc_auprc),
                       median_auprc_FRASER = median(FRASER),
                       median_auprc_FRASER2 = median(FRASER2)), 
                   by="snptype"],
    caption = 'Comparison of AUPRC of FRASER2 and FRASER on splice-disrupting variants on median across tissues',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)



nr_out_gtex <- readRDS(snakemake@input$nr_outliers_per_sample_gtex)$data
#' 
#' ### Number of outliers per sample (FRASER vs FRASER2) 
#' Total nr of samples across all tissues: `r nr_out_gtex[, .N, by="method"][,unique(N)]`
#' #' 
#' FRASER: median nr of outliers across all GTEx samples: `r nr_out_gtex[method=="FRASER", median(nr_outliers_per_sample)]`
#' #' 
#' FRASER 2.0: median nr of outliers across all GTEx samples: `r nr_out_gtex[method=="FRASER2", median(nr_outliers_per_sample)]`
#' #' 
#' Improvement of FRASER 2.0 compared to FRASER:
fc_dt <- dcast(nr_out_gtex[, median(nr_outliers_per_sample), by="tissue,method"], tissue ~ method, value.var="V1")[, .(tissue, FRASER, FRASER2, fc_f1_over_f2=FRASER/FRASER2, fc_f2_over_f1=FRASER2/FRASER)]
fc_dt
fc_dt[, .(median_fc=median(fc_f1_over_f2), mean_fc=mean(fc_f1_over_f2), sd_fc=sd(fc_f1_over_f2))]
fc_dt[, .(median_fc=median(fc_f2_over_f1), mean_fc=mean(fc_f2_over_f1), sd_fc=sd(fc_f2_over_f1))]

#'
#' # Reproducibility analysis stats
#'
required_tissues_per_gene <- snakemake@config$required_tissues_per_gene
required_tissues_per_individual <- snakemake@config$required_tissues_per_individual
#' Minimum nr of tissues per gene: `r required_tissues_per_gene`
#'
#' Minimum nr of available tissues per individual: `r required_tissues_per_individual`
#'
repro_res <- fread(snakemake@input$reproducibility_table)
DT::datatable(
    repro_res[totalTested >= snakemake@config$required_tissues_per_gene, .(tested_genes = uniqueN(geneID)), by="Method"][, Method := gsub("_p$", "", Method)],
    caption = 'Number of genes tested for each method',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)
DT::datatable(
    repro_res[, .(tested_individuals = uniqueN(subjectID)), by="Method"][, Method := gsub("_p$", "", Method)],
    caption = 'Number of genes tested for each method',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#' 
#' # Yepez et al. dataset (mito)
#' 
nr_out_rd <- fread(snakemake@input$nr_outliers_per_sample_rd_cohorts)
#'  
#' Number of samples in each dataset:
#'  
nr_out_rd[, .(nr_of_samples = uniqueN(sampleID)), by="dataset_group"]
nr_out_rd[, .(nr_of_samples = uniqueN(sampleID)), by="dataset_group,dataset_label"]
#' Median nr of outliers per sample:
#' 
fc_dt <- dcast(nr_out_rd, dataset_label + dataset_group ~ method, 
               value.var="num_out", 
               fun.aggregate=function(x) as.numeric(median(x)))[, .(dataset_label, dataset_group, FRASER, `FRASER 2.0`, `FRASER 2.0\n(OMIM)`, `FRASER 2.0\n(OMIM + RV)`, fc_f1_over_f2 = FRASER / `FRASER 2.0`, fc_f2_over_f1 = `FRASER 2.0`/ FRASER)]
fc_dt


#+ number of tested genes and introns
fdsFile <- snakemake@input$fraser2_FDR_sub
fds <- loadFraserDataSet(file=fdsFile) 
#' Median tested genes per sample (transcriptome-wide): `r nrow(pVals(fds, type="jaccard", level="gene", filters=list(rho=0.1)))`
#' 
#' Median tested introns per sample (transcriptome-wide): `r nrow(fds)`
#' 

#+ omim and rv only
fdr_dt <- metadata(fds)$FDR_rare_omim
introns_tested_omim_rv <- fdr_dt[, uniqueN(jidx), by="sampleID"]
setnames(introns_tested_omim_rv, "V1", "Introns")
genes_tested_omim_rv <- fdr_dt[, uniqueN(gene_symbol), by="sampleID"]
setnames(genes_tested_omim_rv, "V1", "Genes")
dt <- merge(genes_tested_omim_rv, introns_tested_omim_rv, by="sampleID")
#' Median tested genes per sample (OMIM + RV): `r dt[, median(Genes)]`
#' 
#' Median tested introns per sample (OMIM + RV): `r dt[, median(Introns)]`
#' 

#'
#' ### Pathogenic cases (mito)
#' 
# pathogenic cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj) & !grepl("deletion|cnv", VARIANT_EFFECT),]
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]
#'
#' Number of pathogenic cases: `r patho_sa[, .N]`
#' 
DT::datatable(
    patho_sa,
    caption = 'Pathogenic cases (Yepez et al.)',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#'
#' ### FRASER 2.0 results of pathogenic cases
#'
#+ read in transcriptome-wide (gene-level) FRASER 2.0 results on mito
res_mito_tw <- fread(snakemake@input$results_mito_tw)
res_mito_tw[, tmp := paste(sampleID, hgncSymbol, sep="_")]
#'
#' Number of pathogenic cases found with FRASER 2.0: `r table(patho_tmp %in% res_mito_tw[,tmp])`
#' 
found <- patho_tmp[patho_tmp %in% res_mito_tw[,tmp]]
not_found <- patho_tmp[!patho_tmp %in% res_mito_tw[,tmp]]
DT::datatable(
    res_mito_tw[tmp %in% patho_tmp],
    caption = 'Pathogenic cases (Yepez et al.) found with FRASER 2.0',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)
#'
#' ### Pathogenic cases NOT found with FRASER 2.0: `r not_found`
#' 
fdr_dt[, tmp := paste(sampleID, gene_symbol, sep="_")]
DT::datatable(
    fdr_dt[tmp %in% not_found][, .SD[which.min(pval)] ,by="tmp"],
    caption = 'Pathogenic cases (Yepez et al.) found with FRASER 2.0',
    options=list(scrollX=TRUE),
    escape=FALSE,
    filter = 'top'
)

#'
#' ### Power analysis (mito)
#' 
res_tps <- fread(snakemake@input$power_analysis_full_res_patho)
res_power_analysis <- res_tps[padjustGene <= 0.1, .(Noutliers=.N, prop=.N/length(patho_tmp)), by=.(size, sim, FDR_set)]
#' at sample size 50:
#' Recovered pathogenic cases: `r res_power_analysis[size == 50, max(Noutliers)]`
#' 
#' Proportion of recovered cases: `r res_power_analysis[size == 50, max(prop)]`
#' 
