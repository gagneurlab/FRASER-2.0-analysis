#'---
#' title: seq depth vs nr of outliers (SPOT, LeafcutterMD)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/{dataset}/gtex_seqDepthCor_spot_leafcutterMD.Rds"`'
#'   threads: 5
#'   resources:
#'    - mem_mb: 32000
#'   py:
#'   - |
#'    def get_spot_tissue_clean(wildcards):
#'     t = wildcards.dataset
#'     tissue = t.replace('-_', '')
#'     tissue = tissue.replace('-', '_')
#'     return config["spot_results"] + tissue + "/spot__fullResults.tsv"
#'    def get_leafcutterMD_tissue_clean(wildcards):
#'     t = wildcards.dataset
#'     tissue = t.replace('-_', '')
#'     tissue = tissue.replace('-', '_')
#'     return config["leafcutterMD_results"] + tissue + "/leafcutterMD_testing/results_" + tissue + ".tsv"
#'   input:
#'    - bam_coverage_tsv: '`sm "config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode29/outrider/{dataset}/bam_coverage.tsv"`'
#'    - res_SPOT: '`sm get_spot_tissue_clean`'
#'    - res_LeafcutterMD: '`sm get_leafcutterMD_tissue_clean`'
#'   output:
#'    - combined_res: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/{dataset}/seq_depth_cor_spot_leafcutterMD.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

library(data.table)

# get tissue name
cov <- snakemake@input$bam_coverage_tsv
t <- basename(dirname(cov))
# for each sample, get seq depth
bam_coverage <- fread(cov)
# get nr of aberrant events 
res_spot <- fread(snakemake@input$res_SPOT)
res_spot[, nrOutPerSample:=.SD[gene_fdr <= snakemake@config$fdrCutoff,.N], by="SAMPLE_ID"]
setnames(res_spot, "SAMPLE_ID", "sampleID")
# get nr of aberrant events (jaccard metric)
res_lmd <- fread(snakemake@input$res_LeafcutterMD)
res_lmd[, nrOutPerSample:=.SD[padj <= snakemake@config$fdrCutoff,.N], by="sample"]
setnames(res_lmd, "sample", "sampleID")

# merge coverage with nr outliers
coverage_dt <- merge(res_spot[, unique(nrOutPerSample), by="sampleID"], bam_coverage, by="sampleID")
setnames(coverage_dt, "V1", "SPOT")
coverage_dt <- merge(coverage_dt, res_lmd[, unique(nrOutPerSample), by="sampleID"])
setnames(coverage_dt, "V1", "LeafcutterMD")
coverage_dt <- melt(coverage_dt, id.vars=c("sampleID", "record_count"), value.name="nr_outliers", variable.name="method")    
coverage_dt[, tissue:=t]
                                  
coverage_dt[, nr_outliers_plus_one := nr_outliers+1]

# save
fwrite(coverage_dt, file=snakemake@output$combined_res, sep="\t")

