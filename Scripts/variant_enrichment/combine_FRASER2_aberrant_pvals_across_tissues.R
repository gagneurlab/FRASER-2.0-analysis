#'---
#' title: Compile FRASER2 pvalues across all tissues (at FDR <= 0.1)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/FRASER2_allTissues_pvals_aberrant_for_rv.Rds"`'
#'   py:
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - FRASER_p: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/FRASER_pval_aberrant_matrix__delta0.3.tsv.gz", dataset=config["tissues_for_reproducibility"])`'
#'     - FRASER2_p: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_25_minN10/PCA__pc0.1/optQ__newFilt/jaccard/pval_aberrant_matrix__delta0.1.tsv.gz", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'     - fraser_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER_allTissues_pvals_aberrant_for_rv.Rds"`'
#'     - fraser2_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER2_allTissues_pvals_aberrant_for_rv.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
library(data.table)
library(GenomicRanges)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ function to combine pval matrices across tissues to one data.table
combine_pvals_across_tissues <- function(input_files, method_name){
    combined_res <- rbindlist(bplapply(input_files, function(in_file, method_name){
        pval_mat <- fread(in_file)
        if(grepl("FRASER2", method_name)){
            t <- basename(dirname(dirname(dirname(dirname(dirname((in_file)))))))
        } else{
            t <- basename(dirname(dirname(dirname(dirname((in_file))))))
        }
        
        # FRASER2 results
        dt <- melt(pval_mat, id.vars="geneID", value.name=method_name, variable.name="sampleID")
        dt[, tissue:=t]
        
        
        setkey(dt, geneID, sampleID, tissue)
        return(dt)
    }, method_name=method_name))
    return(combined_res)
}

#+ read and save combined pvals (FRASER2)
all_res_fraser2 <- combine_pvals_across_tissues(snakemake@input$FRASER2_p, "FRASER2_p")
saveRDS(all_res_fraser2, file=snakemake@output$fraser2_all_rds)

#+ read and save combined pvals (FRASER)
all_res_fraser <- combine_pvals_across_tissues(snakemake@input$FRASER_p, "FRASER_p")
saveRDS(all_res_fraser, file=snakemake@output$fraser_all_rds)
