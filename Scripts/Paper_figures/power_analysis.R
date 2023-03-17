#'---
#' title: Power analysis
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/power_analysis_mito_{FDR_set}.Rds"`'
#'   threads: 8
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'   - res_simulations: '`sm expand(config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/results/" + 
#'                "v29/fraser/{group}/results.tsv", group=config["power_analysis_sim"])`'
#'   - patho_sample_anno: '`sm config["mito_sample_anno"]`' 
#'   output:
#'    - comb_results: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/combined_results_{FDR_set}.tsv"`'
#'    - patho_results: '`sm config["DATADIR"] + "/power_analysis/mito/processed_results/aberrant_splicing/patho_results_{FDR_set}.tsv"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
# .libPaths("~/R/4.1/FRASER2")
# library(FRASER)
library(data.table)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ get FDR set
fdr_set <- snakemake@wildcards$FDR_set
if(fdr_set == "full"){
    fdr_set <- "transcriptome-wide"
}

#+ read in res table and cases to check and extract p values of solved cases
patho_sa <- fread(snakemake@input$patho_sample_anno)
patho_sa <- patho_sa[!is.na(FRASER_padj) & !grepl("deletion|cnv", VARIANT_EFFECT),]
# patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] 
patho_sa[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "C19orf70"]
patho_tmp <- patho_sa[, paste(RNA_ID, KNOWN_MUTATION, sep="_")]

res_all <- rbindlist(bplapply(snakemake@input$res_simulations, FUN=function(res_file){
    res <- fread(res_file)
    res <- res[FDR_set == fdr_set,]
    res[, tmp := paste(sampleID, hgncSymbol, sep="_")]
    res <- res[, .(sampleID, hgncSymbol, pValue, padjust, pValueGene, padjustGene, psiValue, deltaPsi, FDR_set, tmp)]
    group <- basename(dirname(res_file))
    fsplit <- strsplit(group, "_", fixed=TRUE)[[1]]
    res[, group := group]
    res[, size := as.numeric(fsplit[3])]
    res[, sim  := gsub("sim", "run", fsplit[4])]
    return(res)
}))
setorder(res_all, size, sim)

# set name of FDR set for plotting
res_all[FDR_set == "OMIM_RV", FDR_set := "OMIM + RV"]

# create plotting table
res_tps <- res_all[tmp %in% patho_tmp, ]

#+ write combined results table of solved cases
fwrite(res_all, file=snakemake@output$comb_results, sep="\t")
fwrite(res_tps, file=snakemake@output$patho_results, sep="\t")
