#'---
#' title: Get all jaccard variant enrichment results
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/gtex_nrOutliers_delta{delta}.Rds"`'
#'   threads: 3
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - nr_outliers_fraser2:  '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_25_minN10/PCA__pc0.1/optQ/delta{{delta}}/FRASER2_nrOutliers.tsv", dataset=config["tissues_for_reproducibility"])`'
#'     - nr_outliers_fraser1:  '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/minK20_95_minN1/PCA/optQ/FRASER1_nrOutliers.tsv", dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/nrOutliers_FRASER_FRASER2_comparison.html"`'
#'     - combined_tsv: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/nr_outliers_FRASER_FRASER2.tsv"`'
#'     - png: '`sm config["htmlOutputPath"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/nrOutliers_comparison_FRASER_FRASER2.png"`'
#'     - gg_rds: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/PCA__pc0.1/delta{delta}/nrOutliers_comparison_ggplot.Rds"`'
#'   type: noindex
#' output:
#'   html_document
#'---

# #'     - png_newFilt: '`sm config["htmlOutputPath"] + "/GTEx_v8/nrOutliers_FRASER_FRASER2_comparison/nrOutliers_comparison_FRASER_newFilt_FRASER2.png"`' 
# #'     - png_allFilt: '`sm config["htmlOutputPath"] + "/GTEx_v8/nrOutliers_FRASER_FRASER2_comparison/nrOutliers_comparison_FRASER_allFilt_FRASER2.png"`' 
     
saveRDS(snakemake, snakemake@log$snakemake)

library(data.table)
library(ggplot2)
library(ggpubr)

#+ read in all tsvs and combine
dt_combined <- rbindlist(mapply(snakemake@input$nr_outliers_fraser2, snakemake@input$nr_outliers_fraser1, FUN=function(file_f2, file_f1){
    dt_f2 <- fread(file_f2)
    dt_f2 <- dt_f2[, .(sampleID, FRASER2_nrOutlierGenes, FRASER2_nrOutlierJunctions_jaccard, tissue)]
    dt_f1 <- fread(file_f1)
    dt_f1 <- dt_f1[, .(sampleID, FRASER_nrOutlierGenes, FRASER_nrOutlierJunctions_psi5, 
                       FRASER_nrOutlierJunctions_psi3, FRASER_nrOutlierJunctions_theta, tissue)]
    
    dt <- merge(dt_f2, dt_f1, by=c("sampleID", "tissue"), all=TRUE)
    return(dt)
}, SIMPLIFY=FALSE))

#+ write combined file
fwrite(dt_combined, file=snakemake@output$combined_tsv)

# get dt ready for plotting nr of outliers across tissues (FRASER vs FRASER2)
dt_plot <- melt(dt_combined, id.vars=c("sampleID", "tissue"), value.name="nr_outliers_per_sample", variable.name="method")
dt_plot <- dt_plot[method %in% c("FRASER_nrOutlierGenes", "FRASER2_nrOutlierGenes")]
dt_plot[method == "FRASER_nrOutlierGenes", method:="FRASER"]
dt_plot[method == "FRASER2_nrOutlierGenes", method:="FRASER2"]

label_tissue_clean <- function(t){
    t <- gsub("_-_", " ", t)
    t <- gsub("_", " ", t)
    return(t)
}

dt_plot[,tissue:=label_tissue_clean(tissue)]
tissue_order <- dt_plot[,median(nr_outliers_per_sample),by="tissue,method"][method == "FRASER"][order(-V1), tissue]
dt_plot[,tissue:=factor(tissue, levels=tissue_order)]
dt_plot[,method:=factor(method, levels=c("FRASER", "FRASER2"))]

#+ check significance of difference
dt_plot[, wilcox.test(nr_outliers_per_sample ~ method, data = .SD)$p.value, by="tissue"]
dt_plot[, t.test(nr_outliers_per_sample ~ method, data = .SD)$p.value, by="tissue"]

#+ nr of outliers across tissues (FRASER vs FRASER2), fig.height=9, fig.width=18
g <- ggplot(dt_plot, 
            aes(tissue, nr_outliers_per_sample+1, col=method)) +
    geom_boxplot(outlier.size=1) +
    labs(x="GTEx tissue", y="number of outlier genes per sample + 1") +
    scale_y_log10() + 
    # scale_fill_brewer(palette="Paired") +
    scale_color_manual(values=c("lightblue", "purple4")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.margin = unit(c(1,1,1,2.5), "cm"))
# g <- g + ggpubr::stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test")
g

#+ save ggplot
saveRDS(g, file=snakemake@output$gg_rds)

#+ save as png
ggsave(g, filename=snakemake@output$png, height=9, width=18)

# # plot nr of outliers across tissues (FRASER_newFilt vs FRASER2_newFilt)
# dt_plot <- melt(dt_combined, id.vars=c("sampleID", "tissue"), value.name="nr_outliers_per_sample", variable.name="method")
# dt_plot <- dt_plot[method %in% c("FRASER_newFilt_nrOutlierGenes", "FRASER2_newFilt_nrOutlierGenes")]
# dt_plot[method == "FRASER_newFilt_nrOutlierGenes", method:="FRASER"]
# dt_plot[method == "FRASER2_newFilt_nrOutlierGenes", method:="FRASER2"]
# tissue_order <- dt_plot[,median(nr_outliers_per_sample),by="tissue,method"][method == "FRASER"][order(-V1), tissue]
# dt_plot[,tissue:=factor(tissue, levels=tissue_order)]
# dt_plot[,method:=factor(method, levels=c("FRASER2", "FRASER"))]
# 
# #+ nr of outliers across tissues (FRASER_newFilt vs FRASER2), fig.height=6, fig.width=10
# g_newFilt <- ggplot(dt_plot, 
#             aes(tissue, nr_outliers_per_sample+1, fill=method)) +
#     geom_boxplot() +
#     labs(x="GTEx tissue", y="number of outlier genes per sample + 1") +
#     scale_y_log10() + 
#     scale_fill_brewer(palette="Dark2") +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     theme(plot.margin = unit(c(1,1,1,2.5), "cm"))
# g_newFilt
# 
# #+ save as png (newFilter)
# ggsave(g_newFilt, filename=snakemake@output$png_newFilt, height=7, width=12)
# 
# # all FRASER vs FRASER2 in one plot
# dt_plot <- melt(dt_combined, id.vars=c("sampleID", "tissue"), value.name="nr_outliers_per_sample", variable.name="method")
# dt_plot <- dt_plot[method %in% c("FRASER_nrOutlierGenes", "FRASER2_nrOutlierGenes", "FRASER_newFilt_nrOutlierGenes", "FRASER2_newFilt_nrOutlierGenes")]
# dt_plot[method == "FRASER_newFilt_nrOutlierGenes", method:="FRASER_new_filtering"]
# dt_plot[method == "FRASER2_newFilt_nrOutlierGenes", method:="FRASER2_new_filtering"]
# dt_plot[method == "FRASER_nrOutlierGenes", method:="FRASER_old_filtering"]
# dt_plot[method == "FRASER2_nrOutlierGenes", method:="FRASER2_old_filtering"]
# tissue_order <- dt_plot[,median(nr_outliers_per_sample),by="tissue,method"][method == "FRASER_old_filtering"][order(-V1), tissue]
# dt_plot[,tissue:=factor(tissue, levels=tissue_order)]
# dt_plot[,method:=factor(method, levels=c("FRASER2_old_filtering", "FRASER_old_filtering", "FRASER2_new_filtering", "FRASER_new_filtering"))]
# 
# #+ nr of outliers across tissues (all filtering options), fig.height=6, fig.width=10
# g_all_filt <- ggplot(dt_plot, 
#                     aes(tissue, nr_outliers_per_sample+1, fill=method)) +
#     geom_boxplot() +
#     labs(x="GTEx tissue", y="number of outlier genes per sample + 1") +
#     scale_y_log10() + 
#     scale_fill_brewer(palette="Dark2") +
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     theme(plot.margin = unit(c(1,1,1,2.5), "cm"), legend.position="top") 
# g_all_filt
# 
# #+ check significance of difference
# wilcox.test(nr_outliers_per_sample ~ method, data = dt_plot[grep('new', method)])
# wilcox.test(nr_outliers_per_sample ~ method, data = dt_plot[grep('old', method)])
# t.test(nr_outliers_per_sample ~ method, data = dt_plot[grep('new', method)])
# 
# #+ save as png (allFilter)
# ggsave(g_all_filt, filename=snakemake@output$png_allFilt, height=7, width=16)