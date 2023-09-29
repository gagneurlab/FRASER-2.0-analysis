#'---
#' title: Paper figure S1
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_recall_jaccard_vs_psi.Rds"`'
#'   threads: 30
#'   resources:
#'     - mem_mb: 150000
#'   input:
#'     - variant_enrich_comparison_all: '`sm expand(config["DATADIR"] + "/GTEx_v8/{dataset}/plot_rds/FRASER2_enrichment/FRASER_types_vs_jaccard_rv_recall_data_rare{snptype}.Rds", snptype=["SpliceSite", "MMSplice", "SpliceAI", "AbSplice"], dataset=config["tissues_for_detailed_analysis"])`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_recall_jaccard_vs_psi.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_recall_jaccard_vs_psi.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(BiocParallel)
source("src/R/ggplot_theme_for_manuscript.R")
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
maxLineNchar <- 20

#+ read in plots and data for the different panels
tissues <- snakemake@config$tissues_for_detailed_analysis
tissues <- tissues[tissues != "Skin_-_Not_Sun_Exposed_Suprapubic"]
# tissues <- sample(tissues, 7)
tissues <- c("Brain_-_Amygdala", "Brain_-_Cortex", "Heart_-_Left_Ventricle", "Liver", "Lung", "Muscle_-_Skeletal", "Whole_Blood") # fix to use same tissues as in initial submission
tissues <- sort(tissues)
input_files  <- snakemake@input$variant_enrich_comparison_all
input_files <- input_files[which(Reduce("+", lapply(tissues, function(x) grepl(x, input_files))) == 1)]
var_enrich_comp_data <- rbindlist(lapply(input_files, 
                                         function(file){
                                             message(date(), " file: ", file)
                                             data <- readRDS(file)
                                             recall_dt <- data$recallData
                                             snptype <- gsub(".Rds", "", strsplit(basename(file), "_")[[1]][8])
                                             if(snptype == "rareSplicing"){
                                                 recall_dt[, snptype:="rare splice site vicinity\n(VEP)"]
                                             }
                                             if(snptype == "rareSpliceSite"){
                                                 recall_dt[, snptype:="rare direct splice site\nvariant (VEP)"]
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
method_order <- c("psi5", "psi3", "theta", "FRASER", "IntronJaccardIndex")
var_enrich_comp_data[, Method := factor(Method, method_order)]
var_enrich_comp_data[, snptype := factor(snptype, 
                                         levels=c("rare direct splice site\nvariant (VEP)", #"rare splice site vicinity\n(VEP)", 
                                                    "rare MMSplice", 
                                                    "rare SpliceAI", 
                                                    "rare AbSplice"))]
cutoff_dt <- var_enrich_comp_data[!is.na(Cutoff),]
rrdt <- var_enrich_comp_data[is.na(Cutoff),]

getMetricLabels <- function(methods){
    vapply(methods, FUN = function(x) switch(x, 
                                             psi5=c(bquote(psi[5])), 
                                             psi3=c(bquote(psi[3])), 
                                             theta= c(bquote(theta)), 
                                             FRASER = c(bquote(FRASER~"("~psi[5]~","~psi[3]~","~theta~")")), 
                                             IntronJaccardIndex = c(bquote(Intron~Jaccard~Index)) ), 
           FUN.VALUE = c(bquote(psi[3])))
}

metric_labels <- getMetricLabels(methods=method_order)

# source("src/R/enrichment_helper.R")
maxRank <- 100000
dt <- rrdt
dt <- dt[rank<maxRank]
maxPoints <- 1e4
prob4Samp <- min(1, maxPoints / (nrow(dt)/nrow(unique(dt[,.(Method, Type)]))))
dt <- dt[
    rank < 1e3 |
        rank > max(rank) - 1e3 |
        sample(c(TRUE, FALSE), .N, replace = TRUE, prob = c(prob4Samp, 1 - prob4Samp))]


scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

gg <- ggplot(dt[,.(rank=rank, recall=recall, Method, Type, snptype, dataset)],
             aes(rank, recall, col=Method, linetype=Type)) +
    geom_line() + grids()

g_var_enrich_sup <- gg +
    facet_grid(dataset ~ factor(snptype), scales="free_y") +
    labs(x="Top N outliers", y="Recall of rare splice-disrupting candidate variants") +
    grids(color="white") +
    xlim(0, maxRank) +
    ylim(0, rrdt[rank < maxRank, max(recall)]) +
    geom_point(data=cutoff_dt, aes(x=rank, y=recall, color=Method, shape=Cutoff), size=3) +
    geom_abline(intercept=0, slope=rrdt[,1/max(rank)], col="firebrick", linetype="dashed") + 
    scale_color_manual(values=c(RColorBrewer::brewer.pal(5, "Blues")[-1], "darkorchid4"), 
                       labels = metric_labels) +
    scale_x_continuous(breaks=seq(0, maxRank, by=25000),
                       labels=c(0, 25000, "", 75000, ""),
                       limits=c(0, maxRank)) +
    scale_shape_discrete(labels=function(x)parse(text=x)) +
    guides(linetype = "none") + 
    guides(shape=guide_legend(title=ifelse(all(cutoff_dt[,Type == "P-value"]), "Nominal\np-value cutoff", "Cutoff"), order = 2),
           color=guide_legend(title="Splice metric", order = 1)) +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.position="bottom",
          legend.box="vertical",
          legend.title=element_text(size=font_size),
          legend.text=element_text(size=font_size)) +
    cowplot::background_grid(major="xy", minor="xy") 
g_var_enrich_sup

#+ compile supplemental figure
gg_sup <- ggarrange(g_var_enrich_sup, labels=NULL)

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=1.5*page_width, unit=width_unit, dpi=300)
