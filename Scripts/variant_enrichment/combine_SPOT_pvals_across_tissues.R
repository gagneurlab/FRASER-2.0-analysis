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
#'     - spot_all_fdrSignif_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/SPOT_allTissues_pvals_aberrant_for_rv.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
library(data.table)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ function to combine pval matrices across tissues to one data.table
combine_pvals_across_tissues <- function(input_files, method_name="SPOT_p", FDRsignif=FALSE){
    combined_res <- rbindlist(bplapply(input_files, function(in_file, method_name, FDRsignif=FALSE){
    # SPOT results
    spot_res <- fread(in_file)
    t <- basename(dirname(in_file))
    
    if(isTRUE(FDRsignif)){
        spot_res[gene_fdr > snakemake@config$fdrCutoff, gene_p := NA] # only consider FDR significant results
    }
    
    spot_res <- spot_res[, .(SAMPLE_ID, GENE_ID, gene_p)]
    setnames(spot_res, "SAMPLE_ID", "sampleID")
    setnames(spot_res, "GENE_ID", "geneID")
    setnames(spot_res, "gene_p", method_name)
    spot_res[, tissue:=t]
    setkey(spot_res, geneID, sampleID, tissue)
    spot_res <- spot_res[geneID != "",]
    return(spot_res)
    
    }, method_name=method_name, FDRsignif=FDRsignif))
    return(combined_res)
}

#+ read and save combined pvals (SPOT)
all_res <- combine_pvals_across_tissues(input_files=snakemake@input$SPOT_p, method_name="SPOT_p", FDRsignif=FALSE)
saveRDS(all_res, file=snakemake@output$spot_all_rds)


#+ read and save combined pvals (SPOT FDR signficant)
all_res_fdrSignif <- combine_pvals_across_tissues(input_files=snakemake@input$SPOT_p, method_name="SPOT_p", FDRsignif=TRUE)
saveRDS(all_res_fdrSignif, file=snakemake@output$spot_all_fdrSignif_rds)

