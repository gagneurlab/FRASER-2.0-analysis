#'---
#' title: FDR only on genes with rare splice affecting variant
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/12c_FDR_genes_with_variants_MAF{maf}_{varSubset}_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 5
#'  resources:
#'   - mem_mb: 30000
#'  input:
#'   - fraser2_fds: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + config["annotation"] + 
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - sample_anno: '/s/project/genetic_diagnosis/rna_paper/sa_solved.tsv'
#'   - full_sample_anno: '/s/project/prokisch/sample_annotation.tsv'
#'   - rare_omim_variants: '`sm config["mito_processed_data"] + "/variants/rare_omim_variants_MAF{maf}_{varSubset}.Rds"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/fdr_rare_omim_variants_subset_MAF{maf}_{varSubset}.html"`'
#'   - fraser2_FDR_sub: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + 
#'                config["mito_annotation"] + "__FDR_sub_MAF{maf}_{varSubset}" +
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'  type: noindex
#' output:
#'   html_document
#'---

# #'   - gene_name_mapping: '`sm config["gene_name_mapping"]`'
# #'   - pred_gene_level_dir: '`sm expand("/s/project/absplice/data/results/prokisch_WES/splicing_predictions/pred_gene_level/{pred_type}", pred_type=["spliceai", "mmsplice_baseline"])`'
# #'   - raw_pred_dir: '`sm "/s/project/absplice/data/results/prokisch_WES/splicing_predictions/raw_pred/{pred_type}"`'

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
# library(arrow)
library(ggplot2)
library(ggrepel)
library(ggpubr)


fds2File    <- snakemake@input$fraser2_fds
fds1File    <- snakemake@input$fraser1_fds

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations

# Load fds and create a new one
fds2 <- loadFraserDataSet(file=fds2File)

# read in sample anno of subset and known pathogenic variants
sample_anno <- fread(snakemake@input$sample_anno)
sample_anno[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] # symbol was updated
pathogenic_vars <- sample_anno[!is.na(FRASER_padj),]
pathogenic_vars[, sampleGene:=paste(RNA_ID, KNOWN_MUTATION, sep="__")]

# read full proksich sample anno
full_sa <- fread(snakemake@input$full_sample_anno)

# subset to only the 303 samples to consider
maf <- snakemake@wildcards$maf
name(fds2) <- paste0(name(fds2), "__FDR_sub_MAF", maf, "_", snakemake@wildcards$varSubset)
fds2 <- fds2[, sample_anno$RNA_ID]

# # read in gene_name_mapping
# gene_name_mapping <- fread(snakemake@input$gene_name_mapping)

# read in variants
MAF_LIMIT <- 0.001 # not needed as input is already filtered to MAF<=0.001
variants <- readRDS(snakemake@input$rare_omim_variants)

# pred_gene_level_dir <- snakemake@input$pred_gene_level_dir
# pred_gene_level <- rbindlist(lapply(pred_gene_level_dir, function(dir){
#     pred_type <- basename(dir)
#     
#     if(pred_type == "spliceai"){
#         PRED_SCORE_LIMIT <- 0.2
#         score_column <- "delta_score"
#     } else{
#         PRED_SCORE_LIMIT <- 2
#         score_column <- "delta_logit_psi"
#     }
#     # read in gene level predictions
#     pred_gene_level_files <- file.path(dir, list.files(dir, recursive=TRUE, pattern=".csv"))
#     pred_gene_level <- rbindlist(lapply(pred_gene_level_files, function(in_file){
#         dt <- fread(in_file)
#         dt <- dt[abs(get(score_column)) >= PRED_SCORE_LIMIT,]
#         return(dt)
#     }))
#     return(pred_gene_level)
# }), fill=TRUE)
# pred_gene_level[, uniqueN(sample)]
# # pred_gene_level[, .N, by="sample"][, hist(log10(N+1), breaks=30)]

# # read in raw predictions
# raw_pred_dir <- snakemake@input$raw_pred_dir
# raw_pred_files <- file.path(raw_pred_dir, list.files(raw_pred_dir, recursive=TRUE, pattern=ifelse(pred_type == "spliceai", ".parquet", ".csv")))
# raw_pred <- rbindlist(lapply(raw_pred_files, function(in_file){
#     if(grepl(".parquet", in_file)){
#         dt <- as.data.table(arrow::read_parquet(in_file))
#     } else{
#         dt <- fread(in_file)
#     }
#     dt <- dt[abs(get(score_column)) > PRED_SCORE_LIMIT,]
#     return(dt)
# }))
# raw_pred[, .N]
# raw_pred <- raw_pred[!duplicated(raw_pred), ]
# raw_pred[, .N]

# get list of geneIDs with rare predicted splice affecting variant per sample
all_rna_ids <- sample_anno$RNA_ID
genes_with_var <- bplapply(all_rna_ids, function(sid){
    message(date(), ": ", sid)
    dna_id <- full_sa[RNA_ID == sid, DNA_ID]
    # genes <- pred_gene_level[sample == dna_id, unique(gene_id)]
    genes <- variants[sample == dna_id, unique(hgncid)]
    # map gene names to annotation used in fds object
    # sapply(genes, function(gid) gene_name_mapping[grepl(gid, gene_id), gene_name])
})
names(genes_with_var) <- all_rna_ids
head(lapply(genes_with_var, length))

# overlap with genes present in fds object
fds_genes_with_var <- lapply(genes_with_var, function(genes){
    # fds_genes <- genes[sapply(genes, function(g) any(grepl(g, mcols(fds2, type="j")$hgnc_symbol)) )]
    fds_genes <- genes[genes %in% rownames(metadata(fds2)$pvaluesBetaBinomial_gene_rho0.1_jaccard)]
    return(fds_genes)
})
head(fds_genes_with_var)
dt_nr_genes_to_test <- melt(rbind(data.table(type="all genes", as.data.table(lapply(genes_with_var, length))), 
                                  data.table(type="genes with junctions in fds", as.data.table(lapply(fds_genes_with_var, length))) ),
                            id.vars="type", value.name="nr_genes_to_test", variable.name="sampleID")

rnas_no_variants <- dt_nr_genes_to_test[type == "all genes" & nr_genes_to_test == 0, sampleID]
table(sample_anno$RNA_ID %in% rnas_no_variants)

# as histogram
ggplot(dt_nr_genes_to_test, aes(nr_genes_to_test+1)) + 
    facet_wrap(~ type) +
    geom_histogram(bins=50) + 
    scale_x_log10() +
    annotation_logticks(sides="b") +
    # labs(x=paste0("log10(nr of genes with rare ", paste(basename(pred_gene_level_dir), collapse=" or "), " variant +1)")) +
    labs(x="log10(nr of omim genes with rare variants + 1)") +
    theme_classic() 
# only for pathogenic cases
ggplot(dt_nr_genes_to_test[sampleID %in% pathogenic_vars$RNA_ID], aes(nr_genes_to_test+1)) + 
    facet_wrap(~ type) +
    geom_histogram(bins=10) +
    # geom_density() + 
    scale_x_log10() +
    annotation_logticks(sides="b") +
    # labs(x=paste0("log10(nr of genes with rare ", paste(basename(pred_gene_level_dir), collapse=" or "), " variant +1)")) +
    labs(x="log10(nr of omim genes with rare variants + 1)") +
    theme_classic() 
# scatterplot overall nr of genes with var vs genes with var present in fds
ggscatter(dcast(dt_nr_genes_to_test, formula="sampleID ~ type", value.var="nr_genes_to_test"), 
                  'all genes', 'genes with junctions in fds') +
    geom_abline(intercept=0,slope=1, linetype="dotted") + 
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl")

#+ add genes with rare omim var to fds2 object
metadata(fds2)$genes_rare_omim_var <- fds_genes_with_var

#+ compute FDR on subset of genes with rare var & OMIM
FDR_subset <- rbindlist(mapply(all_rna_ids, fds_genes_with_var, FUN=function(sample_id, genes_to_test_sample){
    message(date(), ": sample = ", sample_id)
    # get idx of junctions corresponding to genes with var
    jidx <- unlist(lapply(genes_to_test_sample, function(gene){
        idx <- which(grepl(gene, mcols(fds2, type="j")$hgnc_symbol))
        names(idx) <- rep(gene, length(idx))
        return(idx)
    }))
    jidx <- sort(jidx[!duplicated(jidx)])
    # retrieve pvalues of junctions to test
    p <- pVals(fds2, type="jaccard", level="junction")[jidx, sample_id]
    # FDR correction
    pa <- p.adjust(p, method="BY")
    # gene level pvals
    dt <- data.table(sampleID=sample_id, pval=p, FDR_sub=pa, gene_symbol=names(jidx), jidx=jidx)
    dt[, pval_gene := min(p.adjust(pval, method="holm")), by=gene_symbol]
    # gene level FDR
    dt2 <- dt[, unique(pval_gene), by="gene_symbol"]
    dt2[, FDR_sub_gene := p.adjust(V1, method="BY")]
    dt <- merge(dt, dt2[, .(gene_symbol, FDR_sub_gene)], by="gene_symbol", all.x=TRUE)
    # return new FDR
    return(dt)
}, SIMPLIFY=FALSE))

#+ add FDR subset info to fds object and save
metadata(fds2)$FDR_rare_omim <- FDR_subset
saveFraserDataSet(fds2, rewrite=TRUE)
