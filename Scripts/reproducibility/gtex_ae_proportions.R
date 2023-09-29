#'---
#' title: Proportion of splice outliers that are expression outliers
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/AE_proportion_minK{k}_{q}_minN{n}_{implementation}_delta{delta}.Rds"`'
#'   threads: 5
#'   resources:
#'     - mem_mb: 32000
#'   py:
#'   - |
#'    def get_spot_tissue_clean(wildcards):
#'     t = config["tissues_for_reproducibility"]
#'     tissues = [sub.replace('-_', '') for sub in t]
#'     tissues = [sub.replace('-', '_') for sub in tissues]
#'     return [config["spot_results"] + tissue + "/spot__fullResults.tsv" for tissue in tissues]
#'   - |
#'    def get_leafcutterMD_tissue_clean(wildcards):
#'     t = config["tissues_for_reproducibility"]
#'     tissues = [sub.replace('-_', '') for sub in t]
#'     tissues = [sub.replace('-', '_') for sub in tissues]
#'     return [config["leafcutterMD_results"] + tissue + "/leafcutterMD_testing/results_" + tissue + ".tsv" for tissue in tissues]
#'   input:
#'     - res_fraser2: '`sm expand(config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_{q}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'     - res_fraser1: '`sm  expand(config["general_data_dir"] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'     - res_outrider: '`sm expand(config["general_data_dir"] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_expression/gencode34/outrider/{dataset}/OUTRIDER_results.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'     - res_spot: '`sm get_spot_tissue_clean`'
#'     - res_leafcutterMD: '`sm get_leafcutterMD_tissue_clean`'
#'     - outliers_venn_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK{k}_{q}_minN{n}/{implementation}/optQ/delta{delta}/combined_outliers_venn.Rds"`'
#'   output:
#'     - ae_proportions_table: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK{k}_{q}_minN{n}/{implementation}/optQ/delta{delta}/AE_proportions.tsv"`'
#'   type: script
#'---


saveRDS(snakemake, snakemake@log$snakemake)

library(data.table)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in res tables
plot_dt <- rbindlist(bpmapply(snakemake@input$res_fraser1,
                              snakemake@input$res_fraser2,
                              snakemake@input$res_outrider, 
                              snakemake@input$res_spot,
                              snakemake@input$res_leafcutterMD,
                              FUN=function(file_f1, file_f2, file_ae, file_spot, file_lmd){
    # read in result tables
    res_f1 <- fread(file_f1)
    res_f2 <- fread(file_f2)
    res_ae <- fread(file_ae)
    res_spot <- fread(file_spot)[gene_fdr <= snakemake@config$fdrCutoff,]
    res_lmd <- fread(file_lmd)[padj <= snakemake@config$fdrCutoff,]
    ts <- basename(dirname(file_f1))
    
    # get proportion of splicing outliers that are also AE
    res_ae[, sampleGene := paste(sampleID, hgncSymbol, sep="_")]
    res_f1[, sampleGene := paste(sampleID, hgncSymbol, sep="_")]
    res_f2[, sampleGene := paste(sampleID, hgncSymbol, sep="_")]
    res_spot[, sampleGene := paste(SAMPLE_ID, GENE_ID, sep="_")]
    res_lmd[, sampleGene := paste(sample, geneID, sep="_")]
    # split between up/down AE outliers
    nr_AE_down_f1 <- res_f1[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N]
    nr_AE_down_f2 <- res_f2[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N]
    nr_AE_down_spot <- res_spot[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N]
    nr_AE_down_lmd <- res_lmd[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N]
    nr_AE_down_total <- res_ae[l2fc < 0, .N]
    nr_AE_up_f1 <- res_f1[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N]
    nr_AE_up_f2 <- res_f2[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N]
    nr_AE_up_spot <- res_spot[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N]
    nr_AE_up_lmd <- res_lmd[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N]
    nr_AE_up_total <- res_ae[l2fc > 0, .N]
    
    prop_down_f1 <- res_f1[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N] / res_f1[, .N]
    prop_down_f2 <- res_f2[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N] / res_f2[, .N]
    prop_down_spot <- res_spot[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N] / res_spot[, .N]
    prop_down_lmd <- res_lmd[sampleGene %in% res_ae[l2fc < 0]$sampleGene, .N] / res_lmd[, .N]
    prop_up_f1 <- res_f1[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N] / res_f1[, .N]
    prop_up_f2 <- res_f2[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N] / res_f2[, .N]
    prop_up_spot <- res_spot[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N] / res_spot[, .N]
    prop_up_lmd <- res_lmd[sampleGene %in% res_ae[l2fc > 0]$sampleGene, .N] / res_lmd[, .N]
        
    dt <- data.table(prop=c(prop_down_f1, prop_down_f2, prop_down_spot, prop_down_lmd, 1, 
                            prop_up_f1, prop_up_f2, prop_up_spot, prop_up_lmd, 1), 
                     nr=c(nr_AE_down_f1, nr_AE_down_f2, nr_AE_down_spot, nr_AE_down_lmd, nr_AE_down_total, 
                          nr_AE_up_f1, nr_AE_up_f2, nr_AE_up_spot, nr_AE_up_lmd, nr_AE_up_total), 
                     ae_type=rep(c("Underexpression outliers", "Overexpression outliers"), each=5),
                     method=rep(c("FRASER", "FRASER 2.0", "SPOT", "LeafCutterMD", "OUTRIDER"), 2), 
                     tissue=ts)
    return(dt)
}, SIMPLIFY=FALSE))
plot_dt[, method:=factor(method, levels=c("LeafCutterMD", "SPOT", "FRASER", "FRASER 2.0", "OUTRIDER"))]
plot_dt[, ae_type:=factor(ae_type, levels=c("Underexpression outliers", "Overexpression outliers"))]

#+ write to file
fwrite(plot_dt, file=snakemake@output$ae_proportions_table)




