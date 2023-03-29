#'---
#' title: Compile LeafcutterMD pvalues across all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/LeafcutterMD_allTissues_pvals_for_rv.Rds"`'
#'   py:
#'   - |
#'    def get_leafcutterMD_tissue_clean(wildcards):
#'     t = config["tissues_for_reproducibility"]
#'     tissues = [sub.replace('-_', '') for sub in t]
#'     tissues = [sub.replace('-', '_') for sub in tissues]
#'     return [config["leafcutterMD_results"] + tissue + "/leafcutterMD_testing/results_" + tissue + ".tsv" for tissue in tissues]
#'   threads: 10
#'   resources:
#'     - mem_mb: 75000
#'   input:
#'     - LeafcutterMD_p: '`sm get_leafcutterMD_tissue_clean`'
#'   output:
#'     - lmd_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/LeafcutterMD_allTissues_pvals_for_rv.Rds"`'
#'     - lmd_all_FDRsignif_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/LeafcutterMD_allTissues_pvals_aberrant_for_rv.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
library(data.table)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ function to combine pval matrices across tissues to one data.table
combine_pvals_across_tissues <- function(input_files, method_name="LeafcutterMD_p", FDRsignif=FALSE){
    comb_res <- rbindlist(bplapply(input_files, function(in_file, method_name, FDRsignif){
        # SPOT results
        lfMD_res <- fread(in_file)
        t <- basename(dirname(in_file))
        
        if(isTRUE(FDRsignif)){
            lfMD_res[padj > snakemake@config$fdrCutoff, pvalue_gene := NA] # only consider FDR significant results
        }
        
        lfMD_res <- lfMD_res[, .(sample, geneID, pvalue_gene)]
        setnames(lfMD_res, "sample", "sampleID")
        setnames(lfMD_res, "pvalue_gene", method_name)
        lfMD_res[, tissue:=t]
        setkey(lfMD_res, geneID, sampleID, tissue)
        lfMD_res <- lfMD_res[geneID != "",]
        return(lfMD_res)
        
    }, method_name=method_name, FDRsignif=FDRsignif))
    return(comb_res)
}
#+ read and save combined pvals (LeafcutterMD)
all_res <- combine_pvals_across_tissues(input_files=snakemake@input$LeafcutterMD_p, method_name="LeafcutterMD_p", FDRsignif=FALSE)
saveRDS(all_res, file=snakemake@output$lmd_all_rds)

#+ read and save combined pvals (LeafcutterMD, only FDR significant)
all_res_fdrSignif <- combine_pvals_across_tissues(input_files=snakemake@input$LeafcutterMD_p, method_name="LeafcutterMD_p", FDRsignif=TRUE)
saveRDS(all_res_fdrSignif, file=snakemake@output$lmd_all_FDRsignif_rds)
