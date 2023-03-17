#'---
#' title: Compare FDR on subset of genes to genome-wide FDR
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/12d_FDR_comparison_MAF{maf}_{varSubset}_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 5
#'  resources:
#'   - mem_mb: 20000
#'  input:
#'   - fraser2_FDR_sub: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + 
#'                config["mito_annotation"] + "__FDR_sub_MAF{maf}_{varSubset}" +
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - fraser2_fds: '`sm config["mito_processed_results"] + 
#'                "/datasets/{implementation}/savedObjects/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "--" + 
#'                config["mito_annotation"] + 
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   - fraser2_res_gene: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - fraser2_res_junc: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_per_junction_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - fraser1_res_gene_old: '`sm config["mito_fraser1_results"] + "/results/" + config["mito_annotation"] + "/fraser/fib/results.tsv"`'
#'   - sample_anno: '`sm config["mito_sample_anno"]`'
#'   - full_sample_anno: '`sm config["mito_full_sample_anno"]`'
#'   - omim_genes: '`sm config["omim_genes"]`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/fdr_subset_comparison_MAF{maf}_{varSubset}_F1_old_filter.html"`'
#'   - venn: '`sm expand(config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/venn_qopt_deltaJaccard{delta}_FDR_subset_MAF{maf}_{varSubset}_F1_old_filter.png", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - ggplots: '`sm config["mito_processed_results"] + 
#'                  "/ggplot_rds/{implementation}/" + config["mito_annotation"] + "/" +
#'                  config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                  "/FDR_subset_MAF{maf}_{varSubset}_ggplots_F1_old_filter.Rds"`'
#'  type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load packages
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggvenn)
# library(VennDiagram)
library(UpSetR)
library(ggupset)

res_gene_fraser1 <- fread(snakemake@input$fraser1_res_gene_old)

#+ load FRASER2 objects
fds2_FDRsub <- loadFraserDataSet(file=snakemake@input$fraser2_FDR_sub)
fds2_FDRfull <- loadFraserDataSet(file=snakemake@input$fraser2_fds)

#+ read in list of omim genes
omim_genes <- fread(snakemake@input$omim_genes)

#+ get nr of outliers per sample on FDR subset
FDR_cutoff <- snakemake@config$padjCutoff
FDR_subset <- metadata(fds2_FDRsub)$FDR_rare_omim
nr_junc_outliers_per_sample <- FDR_subset[, .SD[FDR_sub <= FDR_cutoff, .N], by="sampleID"]
nr_junc_outliers_per_sample[, summary(V1)]
nr_gene_outliers_per_sample <- FDR_subset[, .SD[FDR_sub_gene <= FDR_cutoff, uniqueN(gene_symbol)], by="sampleID"]
nr_gene_outliers_per_sample[, summary(V1)]

#+ read in FRASER2 results (genome-wide FDR)
res_gene <- fread(snakemake@input$fraser2_res_gene)
res_junc <- fread(snakemake@input$fraser2_res_junc)
sample_anno <- fread(snakemake@input$sample_anno)
sample_anno[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] # symbol was updated
pathogenic_vars <- sample_anno[!is.na(FRASER_padj),]
pathogenic_vars[, sampleGene:=paste(RNA_ID, KNOWN_MUTATION, sep="__")]
pathogenic_vars <- pathogenic_vars[!grepl("deletion|cnv", VARIANT_EFFECT)]

res_gene_fraser1 <- res_gene_fraser1[sampleID %in% sample_anno$RNA_ID,]
res_gene_fraser1[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]

#+ subset and get sample-gene pairs
res_gene <- res_gene[sampleID %in% sample_anno$RNA_ID,]
res_gene[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]
res_junc <- res_junc[sampleID %in% sample_anno$RNA_ID,]
res_junc[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]

# remove duplicate entries for same outlier when overlapping several genes
res_gene <- res_gene[!duplicated(res_gene, by=c("sampleID","seqnames","start","end","strand")), ]

#+ compare nr of splicing outliers to genome-wide FDR
genome_wide_outlier_junc <- res_junc[, .SD[padjust <= FDR_cutoff, .N], by="sampleID"]
genome_wide_outlier_gene <- res_gene[, .SD[padjustGene <= FDR_cutoff, .N], by="sampleID"]

#+ compare nr of junction-level splicing outliers to genome-wide FDR
outliers_junc <- merge(genome_wide_outlier_junc, nr_junc_outliers_per_sample, by="sampleID", all=TRUE)
setnames(outliers_junc, "V1.x", "genome_wide")
setnames(outliers_junc, "V1.y", "rare_variants_in_omim_gene")
outliers_junc[is.na(genome_wide), genome_wide := 0]
outliers_junc[is.na(rare_variants_in_omim_gene), rare_variants_in_omim_gene := 0]
outliers_junc[, summary(genome_wide)]
outliers_junc[, summary(rare_variants_in_omim_gene)]
# outliers_junc
# gg_scatter_junc <- ggscatter(outliers_junc, x="genome_wide", y="rare_variants_in_omim_gene") +
#     geom_abline(slope=1, intercept=0, linetype="dotted") + 
#     labs(x="junction-level splicing outliers per sample (genome-wide FDR)",
#          y="junction-level splicing outliers per sample\n(FDR only on junctions in OMIM genes with a rare variant)") #+
#     # geom_vline(xintercept=outliers_junc[, median(genome_wide)], linetype="dotted", col="firebrick") +
#     # geom_hline(yintercept=outliers_junc[, median(rare_variants_in_omim_gene)], linetype="dotted", col="firebrick")
gg_scatter_junc <- ggplot(outliers_junc, aes(x=genome_wide + 1, y=rare_variants_in_omim_gene + 1) ) +
    geom_hex(bins = 50) +
    geom_abline(slope=1, intercept=0, linetype="dotted") + 
    geom_hline(yintercept=outliers_junc[,median(rare_variants_in_omim_gene)], linetype="dashed", col="firebrick") +
    geom_vline(xintercept=outliers_junc[,median(genome_wide)], linetype="dashed", col="firebrick") +
    labs(x="FRASER2 splicing outliers per sample + 1\n(transcriptome-wide)",
         y="FRASER2 splicing outliers per sample + 1\n(only on junctions in OMIM genes with a rare variant)") +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl") +
    scale_fill_continuous(type = "viridis") +
    #scale_fill_distiller(palette= "Spectral") +
    theme_classic()
gg_scatter_junc

#+ compare nr of gene-level splicing outliers to genome-wide FDR
outliers_gene <- merge(genome_wide_outlier_gene, nr_gene_outliers_per_sample, by="sampleID", all=TRUE)
setnames(outliers_gene, "V1.x", "genome_wide")
setnames(outliers_gene, "V1.y", "rare_variants_in_omim_gene")
outliers_gene[is.na(genome_wide), genome_wide := 0]
outliers_gene[is.na(rare_variants_in_omim_gene), rare_variants_in_omim_gene := 0]
outliers_gene[, summary(genome_wide)]
outliers_gene[, summary(rare_variants_in_omim_gene)]
# outliers_gene
# gg_scatter_gene <- ggscatter(outliers_gene, x="genome_wide", y="rare_variants_in_omim_gene") +
#     geom_abline(slope=1, intercept=0, linetype="dotted") +
#     labs(x="gene-level splicing outliers per sample (genome-wide FDR)",
#          y="gene-level splicing outliers per sample\n(FDR only on junctions in OMIM genes with a rare variant)") #+
#     # geom_vline(xintercept=outliers_junc[, median(genome_wide)], linetype="dotted", col="firebrick") +
#     # geom_hline(yintercept=outliers_junc[, median(rare_variants_in_omim_gene)], linetype="dotted", col="firebrick")
gg_scatter_gene <- ggplot(outliers_gene, aes(x=genome_wide + 1, y=rare_variants_in_omim_gene + 1) ) +
    geom_hex(bins = 50) +
    geom_abline(slope=1, intercept=0, linetype="dotted") +
    geom_hline(yintercept=outliers_gene[,median(rare_variants_in_omim_gene)], linetype="dashed", col="firebrick") +
    geom_vline(xintercept=outliers_gene[,median(genome_wide)], linetype="dashed", col="firebrick") +
    labs(x="FRASER2 splicing outliers per sample + 1\n(transcriptome-wide FDR)",
         y="FRASER2 splicing outliers per sample + 1\n(only junctions in OMIM genes with a rare variant)") +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl") +
    scale_fill_continuous(type = "viridis") +
    #scale_fill_distiller(palette= "Spectral") +
    theme_classic()
gg_scatter_gene

ggplots <- list()
ggplots[["FDR_comparison_gene"]] <- gg_scatter_gene
ggplots[["FDR_comparison_junc"]] <- gg_scatter_junc


#+ boxplots of outliers per sample
outliers_all <- merge(genome_wide_outlier_gene, nr_gene_outliers_per_sample, by="sampleID", all=TRUE)
setnames(outliers_all, "V1.x", "FRASER2")
setnames(outliers_all, "V1.y", "FRASER2\n(OMIM + RV)")
fraser1_outliers_gene <- res_gene_fraser1[, .SD[padjustGene <= FDR_cutoff, .N], by="sampleID"]
fraser1_outliers_gene[, summary(V1)]
outliers_all <- merge(outliers_all, fraser1_outliers_gene, by="sampleID", all=TRUE)
setnames(outliers_all, "V1", "FRASER")
outliers_all[is.na(FRASER), FRASER := 0]
outliers_all[is.na(FRASER2), FRASER2 := 0]
outliers_all[is.na(`FRASER2\n(OMIM + RV)`), `FRASER2\n(OMIM + RV)` := 0]
outliers_all <- melt(outliers_all, id.vars=c("sampleID"), variable.name="method", value.name="num_out")
outliers_all[, method:=factor(method, levels=c("FRASER", "FRASER2", "FRASER2\n(OMIM + RV)"))]
g_boxplot <- ggplot(outliers_all, aes(x=method, y=num_out + 1, col=method) ) +
    geom_violin() +
    geom_boxplot(width=0.25, alpha=0.2) + # , color="black"
    labs(x="method",
         y="Splicing outliers per sample + 1") +
    scale_y_log10() +
    annotation_logticks(sides="l") +
    labs(x="") +
    # scale_fill_brewer(palette="Paired") + 
    scale_color_manual(values=c("lightblue", "purple4", "purple1")) +
    # scale_fill_manual(values=c("lightblue", "purple4", "purple4")) + 
    theme_pubr() + 
    theme(legend.position="none")
g_boxplot

ggplots[["boxplots_numOut_gene"]] <- g_boxplot

#+ Venn diagram on FDR subset
# FDR_subset[FDR_sub <= FDR_cutoff & sampleID %in% pathogenic_vars$RNA_ID,]
# FDR_subset[FDR_sub_gene <= FDR_cutoff & sampleID %in% pathogenic_vars$RNA_ID &
#                gene_symbol %in% pathogenic_vars$KNOWN_MUTATION,][, .SD[pval==min(pval)], by="gene_symbol,sampleID"]
FDR_subset[, sampleGene:=paste(sampleID, gene_symbol, sep="__")]
FDR_subset[sampleGene %in% pathogenic_vars$sampleGene,][, .SD[pval==min(pval)], by="gene_symbol,sampleID"]
ggplots[["dt_patho_pvals"]] <- FDR_subset[sampleGene %in% pathogenic_vars$sampleGene, .SD[which.min(pval),], by='sampleGene']
# pathogenic genes not present in subset for FDR
not_founc <- pathogenic_vars[!sampleGene %in% FDR_subset$sampleGene, .(RNA_ID, KNOWN_MUTATION)]
not_founc[, is_omim_gene:=(KNOWN_MUTATION %in% omim_genes$gene_v29 | KNOWN_MUTATION %in% c("MICOS13", "PTCD3", "MRPS25"))]
ggplots[["dt_patho_not_found"]] <- not_founc
not_founc

# venn_ls <- list("FRASER2 outliers\nFDR_rareOMIM"=FDR_subset[FDR_sub_gene <= FDR_cutoff,][, .SD[pval==min(pval)], by="gene_symbol,sampleID"][, sampleGene],
#                 "FRASER2 outliers\nFDR_genomewide"=res_gene$sampleGene,
#                 "pathogenic\nsplicing"=pathogenic_vars[,sampleGene],
#                 "aberrantly\nexpressed"=pathogenic_vars[AE_signif == TRUE, sampleGene])
# all_vals <- union(union(venn_ls[["FRASER2 outliers\nFDR_rareOMIM"]], venn_ls[["FRASER2 outliers\nFDR_genomewide"]]), venn_ls[["pathogenic\nsplicing"]])
venn_ls <- list("FRASER"=res_gene_fraser1$sampleGene,
                "FRASER2"=res_gene$sampleGene,
                "FRASER2_rareOMIM"=FDR_subset[FDR_sub_gene <= FDR_cutoff,][, .SD[pval==min(pval)], by="gene_symbol,sampleID"][, sampleGene],
                "pathogenic_splicing"=pathogenic_vars[,sampleGene],
                "aberrantly_expressed"=pathogenic_vars[AE_signif == TRUE, sampleGene])
# res_all$pathogenicSplicing[!res_all$pathogenicSplicing %in% res_all$FRASER2]

all_vals <- union(union(venn_ls[["FRASER"]], union(venn_ls[["FRASER2_rareOMIM"]], venn_ls[["FRASER2"]])), venn_ls[["pathogenic_splicing"]])
venn_dt <- data.table(values=all_vals, 
                      FRASER= all_vals %in% venn_ls[["FRASER"]],
                      FRASER2= all_vals %in% venn_ls[["FRASER2"]],
                      "FRASER2 (OMIM + RV)"= all_vals %in% venn_ls[["FRASER2_rareOMIM"]],
                      "pathogenic splice defect"= all_vals %in% venn_ls[["pathogenic_splicing"]]) #,
                      # "aberrantly expressed"= all_vals %in% venn_ls[["aberrantly_expressed"]])
upset_dt <- melt(venn_dt, id.vars="values", variable.name="Set", value.name="present")
upset_dt <- upset_dt[present == TRUE,]
upset_dt[, comb:=list(list(Set)), by="values"]
upset_dt[, Set:=NULL]
upset_dt[, present:=NULL]
upset_dt <- upset_dt[!duplicated(values)]
# upset_dt[, comb2:=paste(Set, collapse=","), by="values"]

gg_upset <- ggplot(upset_dt, aes(x=comb)) +
    geom_bar() + 
    theme_pubr() +
    scale_x_upset() + 
    # scale_x_upset(ytrans="log1p") + 
    # scale_y_log10() + 
    # annotation_logticks(sides="l") +
    labs(x="") + 
    geom_text(aes(label = ..count..), stat = 'count', nudge_y = 250)
gg_upset
# upset_plot <-UpSetR::upset(fromList(venn_ls), scale.sets="log10",
#                            order.by="freq",
#                            intersections=list(
#                                list("FRASER", "FRASER", "FRASER2_rareOMIM"),
#                                list("FRASER", "FRASER2"),
#                                list("FRASER", "FRASER2_rareOMIM"),
#                                list("FRASER2", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER", "FRASER2", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER", "FRASER2"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER2", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER2"),
#                                list("pathogenic_splicing", "aberrantly_expressed",  "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "aberrantly_expressed", "FRASER"),
#                                list("pathogenic_splicing", "FRASER", "FRASER2", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "FRASER", "FRASER2"),
#                                list("pathogenic_splicing", "FRASER", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "FRASER2", "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "FRASER2"),
#                                list("pathogenic_splicing",  "FRASER2_rareOMIM"),
#                                list("pathogenic_splicing", "FRASER")
#                                ),
#                            empty.intersections=TRUE
#                            # nintersects=15
#                            )
# upset_plot

# venn.plot <- venn.diagram(venn_ls, 
#                           # filename=NULL,
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
#                           fill = c("purple1", "purple4", "darkgreen", "lightgreen"),
#                           # fill = c("cornflowerblue", "darkblue", "darkgreen", "lightgreen"),
#                           # label.col = c("red", "red", "black", "black", "black", 
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
# )
# dev.off(); grid.draw(venn.plot)
# g_venn <- ggvenn(venn_ls,
#                  show_percentage=FALSE,
#                  fill_color = c("purple1", "purple4", "darkgreen", "lightgreen"),
#                  stroke_size = 0.75,
#                  text_size=6)

venn_dt <- data.table(values=all_vals, 
                      "FRASER"=all_vals %in% venn_ls[["FRASER"]],
                      "FRASER2 (OMIM + RV)"=all_vals %in% venn_ls[["FRASER2_rareOMIM"]], 
                      "FRASER2"=all_vals %in% venn_ls[["FRASER2"]], 
                      "pathogenic\nsplice defect"=all_vals %in% venn_ls[["pathogenic_splicing"]]) #,
                      # "aberrantly\nexpressed"=all_vals %in% venn_ls[["aberrantly_expressed"]])
g_venn_w_F1 <- ggplot(venn_dt) +
    geom_venn(aes(A = `FRASER`, B = `FRASER2`, C = `pathogenic\nsplice defect`),
              show_percentage = FALSE,
              fill_color = c("lightblue", "purple4", "darkgreen"),
              stroke_size = 0.75,
              text_size=5) +
    theme_void()
g_venn_w_F1
g_venn_wo_F1 <- ggplot(venn_dt) +
    geom_venn(aes(A = `FRASER2`, B = `FRASER2 (OMIM + RV)`, D = `pathogenic\nsplice defect`),
              show_percentage = FALSE,
              fill_color = c("purple4", "purple1", "darkgreen"),
              stroke_size = 0.75,
              text_size=5) +
    theme_void()
g_venn_wo_F1
g_venn_all <- ggplot(venn_dt) +
    geom_venn(aes(A = `FRASER`, B = `FRASER2`, C = `FRASER2 (OMIM + RV)`, D = `pathogenic\nsplice defect`),
              show_percentage = FALSE,
              fill_color = c("lightblue", "purple4", "purple1", "darkgreen"),
              stroke_size = 0.75,
              text_size=5) +
    theme_void()
g_venn_all

ggsave(g_venn_w_F1, filename=snakemake@output$venn, width=8, height=6)
ggplots[["FDR_venn_all"]] <- g_venn_all
ggplots[["FDR_venn_w_F1"]] <- g_venn_w_F1
ggplots[["FDR_venn_wo_F1"]] <- g_venn_wo_F1
ggplots[["gg_upset"]] <- gg_upset
# ggplots[["upset"]] <- upset_plot

saveRDS(ggplots, file=snakemake@output$ggplots)

# # read full proksich sample anno
# full_sa <- fread(snakemake@input$full_sample_anno)
