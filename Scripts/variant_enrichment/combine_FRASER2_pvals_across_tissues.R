#'---
#' title: Compile FRASER2 pvalues across all tissues
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/FRASER2_allTissues_pvals_for_rv.Rds"`'
#'   py:
#'   threads: 10
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - FRASER_p: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_95_minN1/PCA/FRASER1_optQ/FRASER_pval_matrix__0.3__5__1__FALSE.tsv.gz", dataset=config["tissues_for_reproducibility"])`'
#'     - FRASER2_p: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_25_minN10/PCA__pc0.1/optQ__newFilt/jaccard/pval_matrix__0.1__5__1__FALSE.tsv.gz", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'     - fraser_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER_allTissues_pvals_for_rv.Rds"`'
#'     - fraser2_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER2_allTissues_pvals_for_rv.Rds"`'
#'     - fraser2_noBlacklist_all_rds: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_enrichment/FRASER2_noBlacklist_allTissues_pvals_for_rv.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load needed packages
library(data.table)
library(GenomicRanges)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

# blacklist regions
blacklist_regions <- snakemake@config$blacklist_regions
blacklist_gr <- rtracklayer::import(blacklist_regions, format = "BED")
gene_anno <- fread(snakemake@config$datasets$Whole_Blood$orgdb)
gene_anno_gr <- makeGRangesFromDataFrame(gene_anno, keep.extra.columns=TRUE)
gene_hits <- findOverlaps(blacklist_gr, gene_anno_gr, type="any")
blacklist_genes <- gene_anno[to(gene_hits), gene_name]
# blacklist_genes

#+ function to combine pval matrices across tissues to one data.table
combine_pvals_across_tissues <- function(input_files, method_name, removeBlacklist=FALSE){
    combined_res <- rbindlist(bplapply(input_files, function(in_file, method_name, removeBlacklist){
        pval_mat <- fread(in_file)
        if(grepl("FRASER2", method_name)){
            t <- basename(dirname(dirname(dirname(dirname(dirname((in_file)))))))
        } else{
            t <- basename(dirname(dirname(dirname(dirname((in_file))))))
        }
        
        # FRASER2 results
        dt <- melt(pval_mat, id.vars="geneID", value.name=method_name, variable.name="sampleID")
        dt[, tissue:=t]
        
        if(isTRUE(removeBlacklist)){
            # set pval to 1 to rank last if gene is in blacklist region
            dt[geneID %in% blacklist_genes, (method_name) := 1]
        }
        
        setkey(dt, geneID, sampleID, tissue)
        return(dt)
    }, method_name=method_name, removeBlacklist=removeBlacklist))
    return(combined_res)
}

#+ read and save combined pvals (FRASER2)
all_res_fraser2 <- combine_pvals_across_tissues(snakemake@input$FRASER2_p, "FRASER2_p", removeBlacklist=FALSE)
saveRDS(all_res_fraser2, file=snakemake@output$fraser2_all_rds)

#+ read and save combined pvals (FRASER2)
all_res_fraser2 <- combine_pvals_across_tissues(snakemake@input$FRASER2_p, "FRASER2_noBlacklist_p", removeBlacklist=TRUE)
saveRDS(all_res_fraser2, file=snakemake@output$fraser2_noBlacklist_all_rds)

#+ read and save combined pvals (FRASER)
all_res_fraser <- combine_pvals_across_tissues(snakemake@input$FRASER_p, "FRASER_p", removeBlacklist=FALSE)
saveRDS(all_res_fraser, file=snakemake@output$fraser_all_rds)
