#'---
#' title: Get venn diagrams for different filtering settings and AE implementations
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 2000
#'   input:
#'     - qq_plot_table_f2: '`sm expand(config["mito_processed_data"] + 
#'                "/qqPlot_data/minK{minK}_{quant}_minN{n}/{implementation}/fraser2_qqPlot_data.tsv.gz",
#'                 minK=config["minExpressionInOneSample"], quant=[25], n=config["quantileMinExpression"], implementation=config["AE_implementations"])`'
#'     - qq_plot_table_f1: '`sm config["mito_processed_data"] + "/qqPlot_data/fraser1_qqPlot_data.tsv.gz"`'
#'     - numOut: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/numOutGene.png", 
#'                delta=config["deltaCutoff"], implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"])`'
#'     - venns_FDR_sub: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/venn_qopt_deltaJaccard{delta}_FDR_subset_MAF{maf}_{varSubset}_F1_old_filter.png", 
#'                delta=config["deltaCutoff"], implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"], varSubset=["all", "no_intergenic", "no_utr"], maf=[0.001])`'
#'     - sampleRanks: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/nominal_sampleRank_pathogenicEvents.png", 
#'                implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"])`'
#'     - outlier_annotation: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/outlier_annotation.png", 
#'                implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"])`'
#'     - bam_coverage: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + "/seqDepth_vs_nrOultiers.png", 
#'                implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"])`'
#'     - delta_comparison_png: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/compare_delta_pathogenicEvents.png", 
#'                implementation=config["AE_implementations"],
#'                minK=config["minExpressionInOneSample"], quant=config["quantile"], minN=config["quantileMinExpression"])`'
#'   output:
#'     - all_done: '`sm config["figdir"] + "/FRASER_vs_FRASER2/all_plots_done.txt"`'
#'   type: script
#'---

#+ create dummy file
out <- file.create(snakemake@output$all_done)

