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
#'    - outSvg_a: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_a.svg"`'
#'    - outSvg_b: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_b.svg"`'
#'    - outSvg_c: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_c.svg"`'
#'    - outSvg_d: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_d.svg"`'
#'    - outSvg_e: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_e.svg"`'
#'    - outSvg_f: '`sm config["PAPER_FIGDIR"] + "/heatmaps/heatmap_f.svg"`'
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

# library(ggplot2)
# library(ggpubr)
# library(cowplot)
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
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
    # setnames(sample_info, "hardy_scale", "DTHHRRDY")
    # setnames(sample_info, "age_value", "AGE")
    # setnames(sample_info, "SEX", "GENDER")
    # setnames(sample_info, "rin_number", "RIN")
    sample_info <- merge(sample_info, anno, by.x="sampleID", by.y="RNA_ID", all.x=TRUE)
    col_anno <- data.frame(sample_info[, c("sampleID", names(annotation_cols_2plot)), with=FALSE])
    for(i in names(annotation_cols_2plot)){
        if(any(is.na(col_anno[, i]))){
            col_anno[,i] <- as.factor(paste0("", col_anno[,i]))
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
# annotation_full_col <- SMNABTCHT4plot(annotation_full_col)
# annotation_full_col <- SMRIN4plot(annotation_full_col)
# annotation_full_col <- GENDER4plot(annotation_full_col)
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
    names(res) <- levels(factor(paste0("", annotation_full_col[,get(name)])))
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
    # fds <- fds[seqnames(fds) == "21", 1:30]
    
    # colData <- colData(fds)[,names(annotation_cols_2plot)]
    # colData <- SMNABTCHT4plot(colData)
    # colData <- SMRIN4plot(colData)
    # colData <- GENDER4plot(colData)
    colData <- as.data.table(colData(fds))
    # setnames(sample_info, "hardy_scale", "DTHHRRDY")
    # setnames(sample_info, "age_value", "AGE")
    # setnames(sample_info, "SEX", "GENDER")
    # setnames(sample_info, "rin_number", "RIN")
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
        annotation_legend = TRUE
    )
    
    x
})
heatmap[[1]]

heatmap_after <- bplapply(names(fds_ls), BPPARAM=BPPARAM, function(tissue) {
    
    message(date(), ": ", tissue)
    fds <- fds_ls[[tissue]]
    # fds <- fds[seqnames(fds) == "21", 1:30]
    
    # colData <- colData(fds)[,names(annotation_cols_2plot)]
    # colData <- SMNABTCHT4plot(colData)
    # colData <- SMRIN4plot(colData)
    # colData <- GENDER4plot(colData)
    colData <- as.data.table(colData(fds))
    colData
    # setnames(sample_info, "hardy_scale", "DTHHRRDY")
    # setnames(sample_info, "age_value", "AGE")
    # setnames(sample_info, "SEX", "GENDER")
    # setnames(sample_info, "rin_number", "RIN")
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
        normalized = TRUE,
        annotation_col = colnames(colData),
        annotation_row = NA,
        annotation_colors = annotation_colors,
        sampleCluster = NA,
        # plotMeanPsi = FALSE,
        # plotCov = FALSE,
        annotation_legend = TRUE
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

get_heatmap_grob <- function(gt, spacing=10){
    gt <- gtable_remove_grobs(gt, c(
        "main", "legend", "annotation_legend"))
    gt <- gtable_remove_grobs(gt, c("row_tree"))
    # gt$widths[2] <- unit(1, "grobwidth", data=gt$grobs[[2]])
    gt <- gtable_add_cols(gt, unit(1, "bigpts"))
    gt$widths[4] <- unit(spacing, "bigpts")
    gt
}


displayNames <- tissueNamesClean
displayNames[grepl("Skin Not Sun Exposed", displayNames)] <- "Suprapubic Skin"
g_legend <- ggarrange(get_legend_grob(heatmap[[1]][4]$gtable, spacing=125))
# g_all <- ggarrange(nrow=3, heights=c(6, 6, 6),
#                    ggarrange(nrow=1, ncol=3, widths=c(1,10,10), labels=c("", LETTERS[1:2], ""),
#                              grid.text("", rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[1]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[1]][4]$gtable)),
#                    ggarrange(nrow=1, ncol=3, widths=c(1,10,10), labels=c("", LETTERS[3:4], ""),
#                              grid.text("", rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[2]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[2]][4]$gtable)),
#                    ggarrange(nrow=1, ncol=3, widths=c(1,10,10), labels=c("", LETTERS[5:6], ""),
#                              grid.text("", rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[3]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[3]][4]$gtable))
# )
g_all <- ggarrange(nrow=3, heights=c(6, 6),
                   ggarrange(nrow=1, ncol=2, widths=c(10,10), labels=c(LETTERS[1:2], ""),
                             get_heatmap_grob(heatmap[[1]][4]$gtable),
                             get_heatmap_grob(heatmap_after[[1]][4]$gtable)),
                   ggarrange(nrow=1, ncol=2, widths=c(10,10), labels=c(LETTERS[3:4], ""),
                             get_heatmap_grob(heatmap[[2]][4]$gtable),
                             get_heatmap_grob(heatmap_after[[2]][4]$gtable)),
                   ggarrange(nrow=1, ncol=2, widths=c(10,10), labels=c(LETTERS[5:6], ""),
                             get_heatmap_grob(heatmap[[3]][4]$gtable),
                             get_heatmap_grob(heatmap_after[[3]][4]$gtable))
)
# g_all <- ggarrange(nrow=4, heights=c(6, 6, 6, 4),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[1:2]),
#                              # grid.text(paste(displayNames[1], "samples           "), rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[1]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[1]][4]$gtable)),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[3:4]),
#                              # grid.text(paste(displayNames[2], "samples           "), rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[2]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[2]][4]$gtable)),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[5:6]),
#                              # grid.text(paste(displayNames[3], "samples           "), rot=90, gp=gpar(fontsize=font_size)),
#                              get_heatmap_grob(heatmap[[3]][4]$gtable),
#                              get_heatmap_grob(heatmap_after[[3]][4]$gtable)),
#                    ggarrange(nrow=1,
#                              get_legend_grob(heatmap[[1]][4]$gtable, spacing=0))
# )
# g_all <- ggarrange(nrow=3, heights=c(6, 6, 6),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[1:2]),
#                              heatmap[[1]][4]$gtable,
#                              heatmap_after[[1]][4]$gtable),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[3:4]),
#                              heatmap[[2]][4]$gtable,
#                              heatmap_after[[2]][4]$gtable),
#                    ggarrange(nrow=1, ncol=2, labels=c(LETTERS[5:6]),
#                              heatmap[[3]][4]$gtable,
#                              heatmap_after[[3]][4]$gtable)
# )
g_all



# #+ input
# fds_file <- snakemake@input$fds_file
# 
# #+ load fds
# fds <- loadFraserDataSet(file=fds_file)
# message("Loaded fds.")
# 
# #+ read in gtex annotation
# anno <- fread(snakemake@input$gtex_annotation)
# 
# sample_info <- as.data.table(colData(fds))
# sample_info
# sample_info <- merge(sample_info, anno, by.x="sampleID", by.y="RNA_ID", all.x=TRUE)
# 
# col_anno <- data.frame(sample_info[, .(sampleID, rin_number, SEX, ancestry)])
# for(i in colnames(col_anno)){
#     if(any(is.na(col_anno[, i]))){
#         col_anno[,i] <- as.factor(paste0("", col_anno[,i]))
#     }
#     if(is.integer(col_anno[,i]) && length(levels(as.factor(col_anno[,i]))) <= 10){
#         col_anno[,i] <- as.factor(paste0("", col_anno[,i]))
#     }
# }
# rownames(col_anno) <- col_anno$sampleID
# 
# #+ heatmap before
# p_before <- plotCountCorHeatmap(fds, normalized=FALSE, topN=25000)
# p_before
# 
# #+ heatmap after
# p_after <- plotCountCorHeatmap(fds, normalized=TRUE, topN=25000,
#                                annotation_col=col_anno,
#                                sample)
# p_after

# #+ assemble to figure
# g_all <- ggarrange(p_before$gtable, 
#                    p_after$gtable,
#                    nrow=1, ncol=2,
#                    labels=LETTERS[1:2],
#                    align="hv")
# g_all

#+ save figure as png and pdf
# ggsave(plot=g_all, filename=snakemake@output$outPng, width=1.0*page_width, height=1.1*page_width, unit=width_unit)
# ggsave(plot=g_all, filename=snakemake@output$outPdf, width=1.0*page_width, height=1.1*page_width, unit=width_unit)
ggsave(plot=g_all, filename=snakemake@output$outPng, width=10, height=11)
ggsave(plot=g_all, filename=snakemake@output$outPdf, width=10, height=11)
ggsave(plot=ggarrange(heatmap[[1]][4]$gtable, labels=LETTERS[1]), filename=snakemake@output$outSvg_a, width=7, height=5)
ggsave(plot=ggarrange(heatmap[[2]][4]$gtable, labels=LETTERS[3]), filename=snakemake@output$outSvg_c, width=7, height=5)
ggsave(plot=ggarrange(heatmap[[3]][4]$gtable, labels=LETTERS[5]), filename=snakemake@output$outSvg_e, width=7, height=5)
ggsave(plot=ggarrange(heatmap_after[[1]][4]$gtable, labels=LETTERS[2]), filename=snakemake@output$outSvg_b, width=7, height=5)
ggsave(plot=ggarrange(heatmap_after[[2]][4]$gtable, labels=LETTERS[4]), filename=snakemake@output$outSvg_d, width=7, height=5)
ggsave(plot=ggarrange(heatmap_after[[3]][4]$gtable, labels=LETTERS[6]), filename=snakemake@output$outSvg_f, width=7, height=5)
ggsave(plot=g_legend, filename=snakemake@output$outSvg_legend, width=3, height=7)
