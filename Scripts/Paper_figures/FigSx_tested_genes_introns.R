#'---
#' title: Paper figure S8 (nr of genes to be tested on OMIM+RV)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_testedGenesIntrons.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'   - fraser2_FDR_sub: '`sm config["mito_processed_results"] + "/datasets/PCA__pc0.1/savedObjects/" + 
#'                "fib-minExpr20-quantile0.25-quantCoverage10--gencode34__FDR_sub_MAF0.001_no_utr" +
#'                "/padjBetaBinomial_rho0.1_jaccard.h5"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_testedGenesIntrons.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_testedGenesIntrons.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in fds with gene subsets and extract genes in subset per sample
fdsFile <- snakemake@input$fraser2_FDR_sub
fds <- loadFraserDataSet(file=fdsFile) 

#+ number of tested genes and introns
fdr_dt <- metadata(fds)$FDR_rare_omim
introns_tested_omim_rv <- fdr_dt[, uniqueN(jidx), by="sampleID"]
setnames(introns_tested_omim_rv, "V1", "Introns")
genes_tested_omim_rv <- fdr_dt[, uniqueN(gene_symbol), by="sampleID"]
setnames(genes_tested_omim_rv, "V1", "Genes")
dt <- merge(genes_tested_omim_rv, introns_tested_omim_rv, by="sampleID")
message("median tested genes per sample = ", dt[, median(Genes)])
message("median tested genes per sample = ", dt[, median(Introns)])

#+ plot nr of genes tested per sample
g_genes <- ggplot(dt, aes(Genes)) + 
    geom_histogram(binwidth=10) + 
    geom_vline(xintercept=dt[, median(Genes)], 
               linetype="dotted", col="firebrick") +
    labs(x="Tested genes per sample", y="Count") +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    cowplot::background_grid(major="xy", minor="xy")
# g_genes

g_introns <- ggplot(dt, aes(Introns)) + 
    geom_histogram(binwidth=250) + 
    geom_vline(xintercept=dt[, median(Introns)], 
               linetype="dotted", col="firebrick") +
    labs(x="Tested introns per sample", y="Count") +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    cowplot::background_grid(major="xy", minor="xy")
# g_introns

g_sup <- ggarrange(
    g_genes,
    g_introns,
    labels=LETTERS[1:2],
    font.label=list(size=12, color = "black", face = "bold", family = font),
    nrow=1, ncol=2,
    align="hv")
g_sup

#+ save figure as png and pdf
ggsave(plot=g_sup, filename=snakemake@output$outPng, width=page_width, height=0.33*page_width, unit=width_unit, dpi=300) 
ggsave(plot=g_sup, filename=snakemake@output$outPdf, width=page_width, height=0.33*page_width, unit=width_unit) 
    
