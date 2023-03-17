#'---
#' title: Compile SPOT pvalues across all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/SPOT_allTissues_pvals_for_rv.Rds"`'
#'   py:
#'   - |
#'    def get_spot_tissue_clean(wildcards):
#'     t = config["tissues_for_reproducibility"]
#'     tissues = [sub.replace('-_', '') for sub in t]
#'     tissues = [sub.replace('-', '_') for sub in tissues]
#'     return [config["spot_results"] + tissue + "/spot__fullResults.tsv" for tissue in tissues]
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - SPOT_p: '`sm get_spot_tissue_clean`'
#'   output:
#'     - spot_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/SPOT_allTissues_pvals_for_rv.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
library(data.table)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read and save combined pvals (SPOT)
all_res <- rbindlist(bplapply(snakemake@input$SPOT_p, function(in_file, method_name){
    # SPOT results
    spot_res <- fread(in_file)
    t <- basename(dirname(in_file))
    spot_res <- spot_res[, .(SAMPLE_ID, GENE_ID, gene_p)]
    setnames(spot_res, "SAMPLE_ID", "sampleID")
    setnames(spot_res, "GENE_ID", "geneID")
    setnames(spot_res, "gene_p", method_name)
    spot_res[, tissue:=t]
    setkey(spot_res, geneID, sampleID, tissue)
    spot_res <- spot_res[geneID != "",]
    return(spot_res)
    
}, method_name="SPOT_p"))
saveRDS(all_res, file=snakemake@output$spot_all_rds)

