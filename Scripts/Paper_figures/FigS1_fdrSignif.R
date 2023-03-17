#'---
#' title: Paper figure S1 (FDR signif PR curve)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS1_fdrSignif.Rds"`'
#'   threads: 20
#'   resources:
#'     - mem_mb: 75000
#'   input:
#'     - variant_enrich_comparison_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER2_vs_competitors_fdrSignif_rv_recall_data_rare{snptype}.Rds", snptype=["Splicing", "MMSplice", "SpliceAI", "AbSplice"], dataset=config["tissues_for_reproducibility"])`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS1_fdrSignif_tissueSet1.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS1_fdrSignif_tissueSet1.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(data.table)
# .libPaths("~/R/4.1/FRASER2_BB_loss")
# library(FRASER)
library(ggplot2)
library(ggpubr)
# library(gridExtra)
library(cowplot)
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
maxLineNchar <- 20

#+ read in plots and data for the different panels
tissues <- snakemake@config$tissues_for_reproducibility
# tissues <- tissues[tissues != "Skin_-_Not_Sun_Exposed_Suprapubic"]
# tissues <- sort(sample(tissues, 7))
input_files  <- snakemake@input$variant_enrich_comparison_all
input_files <- input_files[which(Reduce("+", lapply(tissues, function(x) grepl(x, input_files))) == 1)]
var_enrich_comp_data <- rbindlist(bplapply(input_files, 
                                         function(file){
                                             message(date(), " file: ", file)
                                             data <- readRDS(file)
                                             recall_dt <- data$recallData
                                             snptype <- gsub(".Rds", "", strsplit(basename(file), "_")[[1]][8])
                                             if(snptype == "rareSplicing"){
                                                 recall_dt[, snptype:="rare splice site vicinity (VEP)"]
                                             }
                                             if(snptype == "rareMMSplice"){
                                                 recall_dt[, snptype:="rare MMSplice"]
                                             }
                                             if(snptype == "rareSpliceAI"){
                                                 recall_dt[, snptype:="rare SpliceAI"]
                                             }
                                             if(snptype == "rareAbSplice"){
                                                 recall_dt[, snptype:="rare AbSplice"]
                                             }
                                             dataset <- basename(dirname(dirname(dirname(file))))
                                             dataset <- gsub("_", " ", gsub("_-_", " ", dataset))
                                             if(nchar(dataset) > maxLineNchar){
                                                 fsplit <- strsplit(dataset, " ", fixed=TRUE)[[1]]
                                                 dataset <- fsplit[1]
                                                 nc <- nchar(dataset)
                                                 for(s in fsplit[-1]){
                                                     nc <- nc + nchar(s)
                                                     if(nc > maxLineNchar){
                                                         dataset <- paste(dataset, s, sep="\n")
                                                         nc <- nchar(s)
                                                     } else{
                                                         dataset <- paste(dataset, s, sep=" ")
                                                     }
                                                     
                                                 }
                                             }
                                             recall_dt[, dataset := dataset]
                                             return(recall_dt)
                                         }) )
method_order <- c("LeafcutterMD", "SPOT", "FRASER", "IntronJaccardIndex", "FRASER2")
var_enrich_comp_data[, Method := factor(Method, method_order)]
var_enrich_comp_data[, snptype:=factor(snptype, levels=c("rare AbSplice", "rare MMSplice",
                                       "rare SpliceAI", "rare splice site vicinity (VEP)"))]
cutoff_dt_full <- var_enrich_comp_data[!is.na(Cutoff) & Method != "IntronJaccardIndex",]
rrdt_full <- var_enrich_comp_data[is.na(Cutoff) & Method != "IntronJaccardIndex",]

getMetricLabels <- function(methods){
    vapply(methods, FUN = function(x) switch(x, 
                                             FRASER = "FRASER", 
                                             FRASER2 = "FRASER2",
                                             IntronJaccardIndex = "FRASER + Intron Jaccard Index",
                                             SPOT = "SPOT",
                                             LeafcutterMD = "LeafcutterMD"), 
           FUN.VALUE = "")
}

metric_labels <- getMetricLabels(methods=method_order[method_order != "IntronJaccardIndex"])

ggplots <- list()

all_tissues <- rrdt_full[, unique(dataset)]
tissues_per_plot <- 8
for(i in seq(0, length(all_tissues)/tissues_per_plot)){
    tissues_sub <- all_tissues[(i*tissues_per_plot+1):(i*tissues_per_plot+tissues_per_plot)]
    tissues_sub <- tissues_sub[!is.na(tissues_sub)]
    
    rrdt <- rrdt_full[dataset %in% tissues_sub]
    cutoff_dt <- cutoff_dt_full[dataset %in% tissues_sub]

    g_pr_fdrSignif_sup <- ggplot(rrdt[Method != "totalPossibleRank"], aes(x=recall, y=precision, col=Method)) +
        facet_grid(dataset ~ factor(snptype)) +
        geom_line() +
        geom_point(data=cutoff_dt, aes(x=recall, y=precision, color=Method, shape=Cutoff), size=3) +
        labs(x="Recall of rare splice altering variants", y="Precision") +
        grids(color="white") +
        ylim(0, 1) +
        scale_color_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4"),
                           labels=metric_labels) +
        scale_shape_discrete(labels=function(x)parse(text=x)) +
        guides(linetype = "none") + 
        guides(shape=guide_legend(title=ifelse(all(cutoff_dt[,Type == "FDR"]), "FDR cutoff", "Cutoff"), order = 2),
               color=guide_legend(title="Method", order = 1, nrow=2))  +
            theme_pubr() +
            theme(legend.position="bottom",
                  legend.box="vertical",
                  legend.title=element_text(size=font_size),
                  legend.text=element_text(size=font_size),
                  axis.title=element_text(face="bold"),
                  text=element_text(size=font_size)) +
            cowplot::background_grid(major="xy", minor="xy")
    ggplots[[paste0("tissue_subset_", i)]] <- g_pr_fdrSignif_sup
}

for(i in seq_along(ggplots)){
    g_pr_fdrSignif_sup <- ggplots[[i]]
    
    #+ compile supplemental figure
    gg_sup <- ggarrange(g_pr_fdrSignif_sup, labels=NULL)
    
    outPng <- gsub("tissueSet1", paste0("tissueSet", i), snakemake@output$outPng)
    outPdf <- gsub("tissueSet1", paste0("tissueSet", i), snakemake@output$outPdf)
    
    #+ save figure as png and pdf
    ggsave(plot=gg_sup, filename=outPng, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
    ggsave(plot=gg_sup, filename=outPdf, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
    
}

