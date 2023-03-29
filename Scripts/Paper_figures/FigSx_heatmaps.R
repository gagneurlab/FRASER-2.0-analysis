#'---
#' title: Paper figure Sx (GTEx heatmaps)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_heatmaps.Rds"`'
#'   threads: 8
#'   resources:
#'     - mem_mb: 64000
#'   input:
#'     - fds_files: '`sm expand(config["DATADIR"] + "/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/{dataset}__optQ__newFilt/predictedMeans_jaccard.h5", dataset=["Muscle_-_Skeletal", "Skin_-_Not_Sun_Exposed_Suprapubic", "Whole_Blood"])`'
#'     - gtex_annotation: '`sm config["gtex_sample_anno"]`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_heatmaps.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_heatmaps.pdf"`'
#'    - outSvg_legend: '`sm config["PAPER_FIGDIR"] + "/heatmaps/legend.svg"`'
#'    - heatmap_rds: '`sm config["DATADIR"] + "/GTEx_v8/heatmaps/minK20_25_minN10/PCA_pc0.1/heatmaps.Rds"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gtable)
library(grid)
library(RColorBrewer)
library(forcats)
register(MulticoreParam(snakemake@threads))

source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

fdsFiles    <- snakemake@input$fds_files[1:3]
datasets    <- basename(dirname(fdsFiles))
workingDir  <- dirname(dirname(dirname(fdsFiles[1])))
ptype       <- "jaccard"
BPPARAM     <- MulticoreParam(snakemake@threads)


#+ echo=FALSE
fds_ls <- lapply(datasets, loadFraserDataSet, dir=workingDir)
tissueNamesClean <- gsub("_", " ", gsub("_-_", " ", gsub("__optQ__newFilt", "", datasets)))
tissueNamesClean[grepl("Skin Not Sun Exposed", tissueNamesClean)] <- "Suprapubic Skin"
names(fds_ls) <- tissueNamesClean

#+ read in gtex annotation
anno <- fread(snakemake@input$gtex_annotation)
anno[, rin_number := cut(rin_number, breaks=seq(floor(min(rin_number)), ceiling(max(rin_number)), by=1))]
anno[, rin_number := gsub("[()]", "", rin_number)]
anno[, rin_number := gsub("\\]", "", rin_number)]
anno[, rin_number := gsub(",", "-", rin_number)]
anno[, age_value := cut(age_value, breaks=seq(floor(min(age_value)), ceiling(max(age_value)), by=10))]
anno[, age_value := gsub("[()]", "", age_value)]
anno[, age_value := gsub("\\]", "", age_value)]
anno[, age_value := gsub(",", "-", age_value)]
anno[hardy_scale == "", hardy_scale := NA]
anno[hardy_scale == "Fast death of natural causes", hardy_scale := "Fast death (natural)"]
anno[hardy_scale == "Violent and fast death", hardy_scale := "Fast death (violent)"]
# anno[hardy_scale == "Intermediate death", hardy_scale := "Intermediate\ndeath"]

#'
#'# create annotation coloring scheme
#'
annotation_cols_2plot <- list(
    hardy_scale   = c("DTHHRDY", "Set2"),
    age_value       = c("AGE",     "white royalblue"),
    SEX    = c("GENDER",  "Paired"),
    rin_number     = c("RIN",     "white seagreen")#,
    # SMNABTCHT = c("BATCH",   "Accent"),
    # SMCENTER  = c("CENTER",  "Set2")
    )

annotation_full_col <- rbindlist(lapply(fds_ls, function(fds){
    sample_info <- as.data.table(colData(fds))
    sample_info <- merge(sample_info, anno, by.x="sampleID", by.y="RNA_ID", all.x=TRUE)
    col_anno <- data.frame(sample_info[, c("sampleID", names(annotation_cols_2plot)), with=FALSE])
    for(i in names(annotation_cols_2plot)){
        if(any(is.na(col_anno[, i]))){
            col_anno[,i] <- as.factor(paste0("", col_anno[,i]))
            col_anno[,i] <- fct_relevel(col_anno[,i], "NA")
        }
        if(is.integer(col_anno[,i]) && length(levels(as.factor(col_anno[,i]))) <= 10){
            col_anno[,i] <- as.factor(paste0("", col_anno[,i]))
        }
    }
    rownames(col_anno) <- col_anno$sampleID
    as.data.frame(col_anno[,names(annotation_cols_2plot)]) 
    }))
colorSets <- sapply(annotation_cols_2plot, "[[", 2)

#'
#' Correct naming
#'
colnames(annotation_full_col) <- sapply(annotation_cols_2plot, "[[", 1)
names(colorSets) <- colnames(annotation_full_col)


#+ manually create annotation colors
needed_colors <- apply(annotation_full_col, 2, function(x){ length(levels(as.factor(x)))})
annotation_colors <- lapply(names(needed_colors), function(name){
    nrValues <- needed_colors[name]
    if(!grepl(" ", colorSets[name])){
        res <- brewer.pal(max(3, nrValues), colorSets[name])
        if(nrValues < 3){
            res <- res[1:nrValues]
        }
        if(any(is.na(annotation_full_col[,get(name)]))){
            res <- c(res, "white")
        }
    } else {
        res <- colorRampPalette(strsplit(colorSets[name], " ")[[1]])(nrValues)
    }
    if(!is.null(levels(annotation_full_col[,get(name)]))){
        names(res) <- levels(annotation_full_col[,get(name)])
    } else{
        names(res) <- levels(factor(paste0("", annotation_full_col[,get(name)])))
    }
    res
})
names(annotation_colors) <- names(colorSets)


#'
#' # Create heatmaps
#'
#+ creating heatmaps
topN <- 30000
topJ <- 10000
heatmap <- bplapply(names(fds_ls), BPPARAM=BPPARAM, function(tissue) {
    
    message(date(), ": ", tissue)
    fds <- fds_ls[[tissue]]
    colData <- as.data.table(colData(fds))
    colData <- merge(colData, anno, by.x="sampleID", by.y="RNA_ID", all.x=TRUE)
    colData <- data.frame(colData[, names(annotation_cols_2plot), with=FALSE])
    rownames(colData) <- colData$sampleID
    colnames(colData) <- sapply(annotation_cols_2plot, "[[", 1)
    
    for(i in colnames(colData)){
        colData[,i] <- factor(colData[,i], levels=levels(
            factor(annotation_full_col[,get(i)])))
    }
    colData(fds) <- cbind(colData(fds)[,c("sampleID", "bamFile")], colData)
    
    x <- plotCountCorHeatmap(
        fds,
        type = ptype,
        # logit = TRUE,
        topN = topN,
        topJ = topJ,
        # plotType = "junctionSample",
        sampleClustering = NA,
        normalized = FALSE,
        annotation_col = colnames(colData),
        annotation_row = NA,
        annotation_colors = annotation_colors,
        sampleCluster = NA,
        # plotMeanPsi = FALSE,
        # plotCov = FALSE,
        annotation_legend = TRUE,
        main=paste0(tissue, " (before)"),
        fontsize=font_size
    )
    
    x
})
heatmap[[1]]

heatmap_after <- bplapply(names(fds_ls), BPPARAM=BPPARAM, function(tissue) {
    
    message(date(), ": ", tissue)
    fds <- fds_ls[[tissue]]
    colData <- as.data.table(colData(fds))
    colData
    colData <- merge(colData, anno, by.x="sampleID", by.y="RNA_ID", all.x=TRUE)
    colData <- data.frame(colData[, names(annotation_cols_2plot), with=FALSE])
    rownames(colData) <- colData$sampleID
    colnames(colData) <- sapply(annotation_cols_2plot, "[[", 1)
    
    for(i in colnames(colData)){
        colData[,i] <- factor(colData[,i], levels=levels(
            factor(annotation_full_col[,get(i)], )))
    }
    colData(fds) <- cbind(colData(fds)[,c("sampleID", "bamFile")], colData)
    
    x <- plotCountCorHeatmap(
        fds,
        type = ptype,
        # logit = TRUE,
        topN = topN,
        topJ = topJ,
        # plotType = "junctionSample",
        sampleClustering = NA,
        normalized = TRUE,
        annotation_col = colnames(colData),
        annotation_row = NA,
        annotation_colors = annotation_colors,
        sampleCluster = NA,
        # plotMeanPsi = FALSE,
        # plotCov = FALSE,
        annotation_legend = TRUE,
        main=paste0(tissue, " (after)"),
        fontsize=font_size
    )
    
    x
})

saveRDS(file=snakemake@output$heatmap_rds, 
        object=list(before=heatmap, after=heatmap_after)
)

get_legend_grob <- function(gt, spacing=20){
    gt <- gtable_remove_grobs(gt, c(
        "main", "col_tree", "row_tree", "matrix", "col_annotation", "col_annotation_names"))
    gt <- gtable_add_rows(gt, unit(1, "bigpts"), pos=0)
    gt$heights[1] <- unit(spacing, "bigpts")
    gt
}

get_heatmap_grob <- function(gt, spacing=25){
    gt <- gtable_remove_grobs(gt, c(
        "legend", "annotation_legend"))
    gt <- gtable_remove_grobs(gt, c("row_tree"))
    gt$widths[1] <- unit(0.9, "grobwidth", data=gt$grobs[[2]])
    gt <- gtable_add_cols(gt, unit(1, "bigpts"))
    gt$widths[4] <- unit(spacing, "bigpts")
    gt
}


g_legend <- ggarrange(get_legend_grob(heatmap[[1]][4]$gtable, spacing=125))
g_all <- ggarrange(nrow=3, heights=c(1,1,1),
                   ggarrange(nrow=1, ncol=5, widths=c(0.1, 1, 0.15, 1, 0.1), 
                             labels=c(LETTERS[1], "", LETTERS[2], ""),
                             font.label=list(size=12, color = "black", face = "bold", family = font),
                             NULL,
                             get_heatmap_grob(heatmap[[1]][4]$gtable),
                             NULL,
                             get_heatmap_grob(heatmap_after[[1]][4]$gtable)),
                   ggarrange(nrow=1, ncol=5, widths=c(0.1, 1, 0.15, 1, 0.1), 
                             labels=c(LETTERS[3], "", LETTERS[4], ""),
                             font.label=list(size=12, color = "black", face = "bold", family = font),
                             NULL,
                             get_heatmap_grob(heatmap[[2]][4]$gtable),
                             NULL,
                             get_heatmap_grob(heatmap_after[[2]][4]$gtable)),
                   ggarrange(nrow=1, ncol=5, widths=c(0.1, 1, 0.15, 1, 0.1), 
                             labels=c(LETTERS[5], "", LETTERS[6], ""),
                             font.label=list(size=12, color = "black", face = "bold", family = font),
                             NULL,
                             get_heatmap_grob(heatmap[[3]][4]$gtable),
                             NULL,
                             get_heatmap_grob(heatmap_after[[3]][4]$gtable))
)
g_all

#+ save figure as png and pdf
ggsave(plot=g_all, filename=snakemake@output$outPng, width=page_width, height=1.5*page_width, unit=width_unit)
ggsave(plot=g_all, filename=snakemake@output$outPdf, width=page_width, height=1.5*page_width, unit=width_unit)
ggsave(plot=g_legend, filename=snakemake@output$outSvg_legend, width=3, height=7)
