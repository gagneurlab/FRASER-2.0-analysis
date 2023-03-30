#'---
#' title: Comparison Jaccard Index to psi3/5 
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/GTEx_v8/{dataset}/minK{k}_{q}_minN{n}/{implementation}/delta{delta}/compare_jaccard_psi.Rds"`'
#'   threads: 4
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'     - res_fraser1: '`sm  config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/results/gencode34/fraser/{dataset}_old_filter/results.tsv"`'
#'     - res_fraser2: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK{k}_{q}_minN{n}/{implementation}/{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv"`'
#'     - fds_fraser2: '`sm config["DATADIR"] + "/GTEx_v8/fds/minK{k}_{q}_minN{n}/{implementation}/savedObjects/{dataset}__optQ__newFilt/fds-object.RDS"`'
#'     - fds_fraser1: '`sm  "config["general_data_dir] + "/gtex_genetic_diagnosis/v8/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}_old_filter--gencode34/fds-object.RDS"`'
#'   output:
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_v8/fraser2_improvements/{dataset}/minK{k}_{q}_minN{n}/{implementation}/delta{delta}/optQ/jaccard_to_psi_comparison.html"`'
#'     - res_fraser1_anno: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/{dataset}/minK{k}_{q}_minN{n}/{implementation}/delta{delta}/optQ/fraser1_res_with_jaccard.tsv"`'
#'   type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ load packages
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

dataset <- snakemake@wildcards$dataset
dataset <- gsub("_", " ", gsub("_-_", " ", dataset))

#+ read in fraser results
res_f1 <- fread(snakemake@input$res_fraser1)
res_f1 <- res_f1[type != "theta", .(sampleID, seqnames, start, end, strand, type, psiValue, deltaPsi)]
res_f1

#+ read in fraser2 fds
fds_f2 <- loadFraserDataSet(file=snakemake@input$fds_fraser2)
dim(fds_f2)
jaccard <- assay(fds_f2, "jaccard")
deltaJaccard <- deltaPsiValue(fds_f2, type="jaccard")

#+ find junction index in fraser2 fds object
rr_f2 <- rowRanges(fds_f2, type="j")
res_gr_f1 <- makeGRangesFromDataFrame(res_f1, keep.extra.columns=TRUE)
hits <- findOverlaps(res_gr_f1, rr_f2, type="equal")
res_f1[from(hits), jidx_f2 := to(hits)]
res_f1

#+ get col index of sampleIDs
s_vec <- seq_len(ncol(fds_f2))
names(s_vec) <- colnames(fds_f2)
res_f1[, sidx_f2 := s_vec[sampleID]]

#+ extract jaccard value
res_f1[!is.na(jidx_f2), jaccardValue := jaccard[cbind(jidx_f2, sidx_f2)] ]
res_f1[!is.na(jidx_f2), deltaJaccard := deltaJaccard[cbind(jidx_f2, sidx_f2)] ]
res_f1

#+ add fraser2 outlier status
res_f2 <- fread(snakemake@input$res_fraser2)
res_f2[, FRASER2_outlier := TRUE]
res_f2 <- res_f2[,.(sampleID, seqnames, start, end, strand, FRASER2_outlier)]
res_f1 <- merge(res_f1, res_f2, by=c("sampleID", "seqnames", "start", "end", "strand"), all.x=TRUE)
res_f1[is.na(FRASER2_outlier), FRASER2_outlier := FALSE]
res_f1[, table(FRASER2_outlier)]
res_f1

#+ relationship of delta psi and delta jaccard by scatter plotting
ggplot(res_f1, aes(deltaPsi, deltaJaccard)) +
    facet_grid(FRASER2_outlier~type, labeller=label_both) +
    geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_hex(binwidth=0.025) +
    labs(x="delta psi", y="delta Intron Jaccard Index", title=dataset)

#+ scatter plot of psi vs jaccard with density, fig.width=8, fig.heigth=4
ggplot(res_f1[,], aes(psiValue, jaccardValue)) +
    facet_wrap(~type) +
    geom_point() +
    geom_hex(binwidth=0.025) +
    labs(x="psi_3/5", y="Intron Jaccard Index", title=dataset)

#+ scatter plot of psi vs jaccard with density in log scale, fig.width=8, fig.heigth=4
ggplot(res_f1[,], aes(psiValue, jaccardValue)) +
    facet_wrap(~type) +
    geom_point() +
    geom_hex(binwidth=0.01) +
    geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
    geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
    scale_fill_gradient(trans = "log10") +
    labs(x="psi_3/5", y="Intron Jaccard Index", title=dataset)

#+ scatter plot of psi vs jaccard without density, fig.width=8, fig.heigth=4
ggplot(res_f1[,], aes(psiValue, jaccardValue)) +
    facet_wrap(~type) +
    geom_point(size=0.75) + #, alpha=0.25
    labs(x="psi_3/5", y="Intron Jaccard Index", title=dataset)

#+ boxplot of jaccard values for fraser1 outliers with large psi3/5 value, fig.width=8, fig.heigth=4
ggplot(res_f1[psiValue >= 0.9], aes(jaccardValue)) +
    facet_wrap(~type) +
    geom_boxplot() + 
    labs(x="Intron Jaccard Index of FRASER1 outliers\nwith psi3/5 >= 0.9",
         title=dataset)
ggplot(res_f1[psiValue >= 0.75], aes(jaccardValue)) +
    facet_wrap(~type) +
    geom_boxplot() + 
    labs(x="Intron Jaccard Index of FRASER1 outliers\nwith psi3/5 >= 0.75", 
         title=dataset)

#+ scatterhist for one tissue, fig.width=6, fig.height=4
library("RColorBrewer")
colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100))))
p <- ggplot(res_f1[type == "psi5"], aes(psiValue, jaccardValue)) +
    geom_point(size=0.75) + #, alpha=0.25
    # geom_hex(binwidth=0.01) +
    geom_bin2d(binwidth=0.03) +
    geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
    geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
    scale_fill_gradient(trans = "log10") +
    labs(x=expression(bold(psi[5])), y="Intron Jaccard Index",
         title=dataset) + 
    colorscale
ggExtra::ggMarginal(p, type = "histogram") # not sure if it will work with the facets

#+ different options for plotting, fig.width=6, fig.height=4
p <- ggplot(res_f1[type == "psi5"], aes(psiValue, jaccardValue)) + 
    ggtitle(dataset) + labs(x="psi_5", y="Intron Jaccard Index") 
# option stat_density2d
p + stat_density2d(h = 0.1, bins = 100, aes( fill = ..level..), geom = "polygon") + 
    colorscale +
    geom_abline(slope=1, intercept=0, col="black") +
    geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
    geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")
# option 2 geom_hex
p + geom_hex() +  
    colorscale +
    geom_abline(slope=1, intercept=0, col="black") +
    geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
    geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")
# option 3 geom_bin2d
p + geom_bin2d() + 
    colorscale +
    geom_abline(slope=1, intercept=0, col="black") +
    geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
    geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")


#+ table of high psi value vs low jaccard value
res_f1[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.9` = psiValue >= 0.9, `jaccard <= 0.1` = jaccardValue <= 0.1)]
res_f1[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.9` = psiValue >= 0.9, `jaccard <= 0.1` = jaccardValue <= 0.1)]
tb_5 <- table(res_f1[, .(`psi5 >= 0.9`, `jaccard <= 0.1`)] )
tb_3 <- table(res_f1[, .(`psi3 >= 0.9`, `jaccard <= 0.1`)] )
tb_5[2,2] / sum(tb_5)
tb_3[2,2] / sum(tb_3)
res_f1[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.85` = psiValue >= 0.85, `jaccard <= 0.15` = jaccardValue <= 0.15)]
res_f1[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.85` = psiValue >= 0.85, `jaccard <= 0.15` = jaccardValue <= 0.15)]
tb_5 <- table(res_f1[, .(`psi5 >= 0.85`, `jaccard <= 0.15`)] )
tb_3 <- table(res_f1[, .(`psi3 >= 0.85`, `jaccard <= 0.15`)] )
tb_5[2,2] / sum(tb_5)
tb_3[2,2] / sum(tb_3)
res_f1[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.75` = psiValue >= 0.75, `jaccard <= 0.25` = jaccardValue <= 0.25)]
res_f1[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.75` = psiValue >= 0.75, `jaccard <= 0.25` = jaccardValue <= 0.25)]
tb_5 <- table(res_f1[, .(`psi5 >= 0.75`, `jaccard <= 0.25`)] )
tb_3 <- table(res_f1[, .(`psi3 >= 0.75`, `jaccard <= 0.25`)] )
tb_5[2,2] / sum(tb_5)
tb_3[2,2] / sum(tb_3)

#+ barplots of table of high psi value vs low jaccard value, fig.width=6, fig.heigth=4
res_f1[!is.na(jidx_f2), `:=`(psiBig = psiValue >= 0.9, jaccardSmall = jaccardValue <= 0.1)]
ggplot(as.data.frame(table(res_f1[, .(psiBig, jaccardSmall, type)])), aes(x=psiBig, y = Freq, fill=jaccardSmall)) + 
    facet_wrap(~type) +
    geom_bar(stat="identity") +
    labs(x="psi >= 0.9", title=dataset) +
    guides(fill=guide_legend(title="jaccard <= 0.1"))

#+ write fraser1 results annotated with jaccard information to file
fwrite(res_f1, file=snakemake@output$res_fraser1_anno, sep="\t")
