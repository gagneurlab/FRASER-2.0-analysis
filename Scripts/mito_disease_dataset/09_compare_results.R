#'---
#' title: Compare with old FRASER results
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/09_compareVenn_f1_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 1
#'  resources:
#'   - mem_mb: 10000
#'  input:
#'   - resultTableGene_jaccard: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableGene_f1_old: '`sm config["mito_fraser1_results"] + "/results/" + config["mito_annotation"] + "/fraser/fib/results.tsv"`'
#'   - resultTableJunction_jaccard: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_per_junction_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableJunction_f1_old:  '`sm config["mito_fraser1_results"] + "/results/" + config["mito_annotation"] + "/fraser/fib/results_per_junction.tsv"`'
#'   - sample_anno: '`sm config["mito_sample_anno"]`'
#'  output:
#'   - ggplots: '`sm config["mito_processed_results"] + 
#'                  "/ggplot_rds/{implementation}/" + config["mito_annotation"] + "/" +
#'                  config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                  "/FRASER2_vs_FRASER1_ggplots.Rds"`'
#'   - patho_sampleRank: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/sampleRank_pathogenicEvents.png"`'
#'   - patho_pvals: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/pvalues_pathogenicEvents.png"`'
#'   - numOutGene: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/numOutGene.png"`'
#'   - numOutJunc: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/numOutJunc.png"`'
#'  type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)
# snakemake <- readRDS("logs/mito/09_compareVenn_PCA_minExpr20-quantile0.95-quantCoverage1.Rds")

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
# library(VennDiagram)
library(ggvenn)
library(ggplot2)

# read in gene-level and junction-level results tables
res_dt_jaccard <- fread(snakemake@input$resultTableGene_jaccard)
res_junction_jaccard <- fread(snakemake@input$resultTableJunction_jaccard)

# old FRASER1 results
res_dt_old <- fread(snakemake@input$resultTableGene_f1_old)
res_junction_old <- fread(snakemake@input$resultTableJunction_f1_old)

# read in sample anno of subset and known pathogenic variants
sample_anno <- fread(snakemake@input$sample_anno)
sample_anno[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] # symbol was updated
pathogenic_vars <- sample_anno[!is.na(FRASER_padj),]
pathogenic_vars[, sampleGene:=paste(RNA_ID, KNOWN_MUTATION, sep="__")]
pathogenic_vars <- pathogenic_vars[!grepl("deletion|cnv", VARIANT_EFFECT)]

# subset and get sample-gene pairs
res_dt_jaccard <- res_dt_jaccard[sampleID %in% sample_anno$RNA_ID,]
res_dt_jaccard[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]
res_dt_old <- res_dt_old[sampleID %in% sample_anno$RNA_ID,]
res_dt_old[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]

# remove duplicate entries for same outlier when overlapping several genes
# res_dt_jaccard[, inOldRes:=(paste(sampleID,seqnames,start,end,strand, sep="_") %in% res_dt_old[,paste(sampleID,seqnames,start,end,strand, sep="_")])]
# res_dt_jaccard[, numGenesPerOutlier:=uniqueN(hgncSymbol), by="sampleID,seqnames,start,end,strand"]
res_dt_jaccard <- res_dt_jaccard[!duplicated(res_dt_jaccard, by=c("sampleID","seqnames","start","end","strand")), ]

# get information about which cases are also expression outliers
ae_cases <- pathogenic_vars[AE_signif == TRUE, ]

# draw Venn diagram of results overlap
res_all <- list("FRASER"=res_dt_old[,sampleGene], 
                "FRASER2"=res_dt_jaccard[,sampleGene], 
                "pathogenicSplicing"=pathogenic_vars[,sampleGene],
                "aberrantExpression"=ae_cases[,sampleGene])
# res_all$pathogenicSplicing[!res_all$pathogenicSplicing %in% res_all$FRASER2]

# venn.plot <- venn.diagram(res_all, 
#                           # filename="Output/venn_test.png",
#                           filename=snakemake@output$venn,
#                           imagetype="png",
#                           disable.logging=TRUE,
#                           width=7500, height=5000, 
#                           col = "black",
#                           # lty = "dotted",
#                           lwd = 4,
#                           # fill = c("yellow", "green", "darkorchid1", "cornflowerblue"),
#                           # fill = c("#A6CEE3", "#1F78B4", "darkgreen", "lightgreen"),
#                           fill = c("#A6CEE3", "purple4", "darkgreen", "lightgreen"),
#                           # label.col = c("red", "red", "black", "red", "red", 
#                           label.col = c("black", "black", "black", "black", "black", 
#                                         "black", "black", "black", "black", "black",
#                                         "black", "black", "black", "black", "black"),
#                           alpha = 0.50,
#                           # cex = 1.5,
#                           fontfamily = "serif",
#                           fontface = "bold",
#                           # cat.col = c("orange", "darkgreen", "darkorchid4", "darkblue"),
#                           cat.cex = 1.5,
#                           cat.fontfamily = "serif"
#                           )
# dev.off(); grid.draw(venn.plot)
all_vals <- union(union(res_all$FRASER, res_all$FRASER2), res_all$pathogenicSplicing)
venn_dt <- data.table(values=all_vals, 
                      "FRASER\noutliers"=all_vals %in% res_all$FRASER, 
                      "FRASER2\noutliers"=all_vals %in% res_all$FRASER2, 
                      "pathogenic\nsplicing"=all_vals %in% res_all$pathogenicSplicing,
                      "aberrantly\nexpressed"=all_vals %in% res_all$aberrantExpression)
g_venn <- ggplot(venn_dt) +
    geom_venn(aes(A = `FRASER\noutliers`, B = `FRASER2\noutliers`, C = `pathogenic\nsplicing`, D=`aberrantly\nexpressed`), 
              show_percentage = FALSE,
              fill_color = c("#A6CEE3", "purple4", "darkgreen", "lightgreen"),
              stroke_size = 0.75,
              text_size=6) +
    theme_void()
# ggsave(g_venn, filename=snakemake@output$venn, width=10, height=6)

ggplots <- list()
ggplots[["g_venn"]] <- g_venn

### outlier rank per sample (pathogenic cases)
res_dt_old[, outlierSampleRank:=frank(padjust, ties.method="min"), by="sampleID"]
res_dt_jaccard[, outlierSampleRank:=frank(padjust, ties.method="min"), by="sampleID"]

sampleRank_dt <- merge(res_dt_old[sampleGene %in% res_all$pathogenic, .(sampleGene, outlierSampleRank)],
                       res_dt_jaccard[sampleGene %in% res_all$pathogenic, .(sampleGene, outlierSampleRank)], all=T, by="sampleGene")
setnames(sampleRank_dt, "outlierSampleRank.x", "FRASER")
setnames(sampleRank_dt, "outlierSampleRank.y", "jaccard")
g_sampleRank <- ggplot(sampleRank_dt, aes(FRASER, jaccard)) + 
    geom_point() + 
    geom_abline(intercept=0, slope=1, linetype="dotted") + 
    labs(x="outlier sample rank (FRASER)", y="outlier sample rank (jaccard)") + 
    theme_bw()
ggsave(g_sampleRank, filename=snakemake@output$patho_sampleRank)
ggplots[["g_sampleRank"]] <- g_sampleRank
# sampleRank_dt <- melt(sampleRank_dt, id.vars="sampleGene", variable.name="metric", value.name="rank")
# ggplot(sampleRank_dt, aes(metric, rank)) + geom_boxplot() + theme_bw()
pval_patho_dt <- merge(res_dt_old[sampleGene %in% res_all$pathogenic, .(sampleGene, pValue)],
                       res_dt_jaccard[sampleGene %in% res_all$pathogenic, .(sampleGene, pValue)], all=T, by="sampleGene")
setnames(pval_patho_dt, "pValue.x", "FRASER")
setnames(pval_patho_dt, "pValue.y", "jaccard")
g_pval_patho <- ggplot(pval_patho_dt, aes(-log10(FRASER), -log10(jaccard))) + 
    geom_point() + 
    geom_abline(intercept=0, slope=1, linetype="dotted") + 
    labs(x="-log10(pval of pathogenic event) (FRASER)", 
         y="-log10(pvalue of pathogenic event) (jaccard)") + 
    theme_bw()
ggsave(g_pval_patho, filename=snakemake@output$patho_pvals)
ggplots[["g_pval_patho"]] <- g_pval_patho
###


# subset 
res_junction_jaccard <- res_junction_jaccard[sampleID %in% sample_anno$RNA_ID,]
res_junction_old <- res_junction_old[sampleID %in% sample_anno$RNA_ID,]

# remove duplicate entries for same outlier when overlapping several genes
res_junction_jaccard <- res_junction_jaccard[!duplicated(res_junction_jaccard, by=c("sampleID","seqnames","start","end","strand")), ]

# get nr of outliers per sample (gene-level)
res_dt_jaccard[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
res_dt_jaccard[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
res_dt_old[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
res_dt_old[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]

plot_dt <- merge(res_dt_jaccard[,unique(numOutliersPerSample), by="sampleID"],
                 res_dt_old[,unique(numOutliersPerSample), by="sampleID"],
                 by="sampleID", all=TRUE)
setnames(plot_dt, "V1.x", "numOutlierGenes_jaccard")
setnames(plot_dt, "V1.y", "numOutlierGenes_FRASER")
plot_dt[is.na(numOutlierGenes_jaccard), numOutlierGenes_jaccard := 0]
plot_dt[is.na(numOutlierGenes_FRASER), numOutlierGenes_FRASER := 0]
# g_numOut_gene <- ggplot(plot_dt, aes(numOutlierGenes_FRASER, numOutlierGenes_jaccard)) +
#     geom_point() +
#     geom_abline(intercept=0, slope=1, linetype="dotted") +
#     geom_hline(yintercept=plot_dt[,median(numOutlierGenes_jaccard)], linetype="dashed", col="firebrick") +
#     geom_vline(xintercept=plot_dt[,median(numOutlierGenes_FRASER)], linetype="dashed", col="firebrick") +
#     labs(x="FRASER outliers (gene-level)", y="FRASER2 outliers (gene-level)") +
#     scale_x_log10() +
#     scale_y_log10() +
#     annotation_logticks(sides="bl") +
#     theme_bw() +
#     theme(text=element_text(size=14))
g_numOut_gene <- ggplot(plot_dt, aes(numOutlierGenes_FRASER+1, numOutlierGenes_jaccard + 1)) +
    geom_hex(bins = 50) +
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    geom_hline(yintercept=plot_dt[,median(numOutlierGenes_jaccard)], linetype="dashed", col="firebrick") +
    geom_vline(xintercept=plot_dt[,median(numOutlierGenes_FRASER)], linetype="dashed", col="firebrick") +
    labs(x="FRASER splicing outliers per sample + 1", y="FRASER2 splicing outliers per sample + 1") +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl") +
    scale_fill_continuous(type = "viridis") +
    #scale_fill_distiller(palette= "Spectral") +
    theme_classic() +
    theme(text=element_text(size=14))
ggsave(g_numOut_gene, filename=snakemake@output$numOutGene)
ggplots[["g_numOut_gene"]] <- g_numOut_gene

# get nr of outliers per sample (junction-level)
res_junction_jaccard[, numOutlierJuncsPerSample := .N, by = "sampleID"]
res_junction_jaccard[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
res_junction_jaccard[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
res_junction_old[, numOutlierJuncsPerSample := .N, by = "sampleID"]
res_junction_old[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
res_junction_old[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]

plot_dt_junc <- merge(res_junction_jaccard[,unique(numOutlierJuncsPerSample), by="sampleID"],
                 res_junction_old[,unique(numOutlierJuncsPerSample), by="sampleID"],
                 by="sampleID", all=TRUE)
setnames(plot_dt_junc, "V1.x", "numOutlierIntrons_jaccard")
setnames(plot_dt_junc, "V1.y", "numOutlierIntrons_FRASER")
plot_dt[is.na(numOutlierGenes_jaccard), numOutlierGenes_jaccard := 0]
plot_dt[is.na(numOutlierGenes_FRASER), numOutlierGenes_FRASER := 0]
g_numOut_junc <- ggplot(plot_dt_junc, aes(numOutlierIntrons_FRASER + 1, numOutlierIntrons_jaccard + 1)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    geom_hline(yintercept=plot_dt_junc[,median(numOutlierIntrons_jaccard)], linetype="dashed", col="firebrick") +
    geom_vline(xintercept=plot_dt_junc[,median(numOutlierIntrons_FRASER)], linetype="dashed", col="firebrick") +
    labs(x="FRASER splicing outliers per sample + 1", y="FRASER2 splicing outliers per sample + 1") +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl") +
    theme_bw()
ggsave(g_numOut_junc, filename=snakemake@output$numOutJunc)
ggplots[["g_numOut_junc"]] <- g_numOut_junc

#+ save ggplots as rds
saveRDS(ggplots, file=snakemake@output$ggplots)
