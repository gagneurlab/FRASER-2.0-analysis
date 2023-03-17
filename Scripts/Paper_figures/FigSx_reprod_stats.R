#'---
#' title: Tested nr of genes /individuals per method (reproducibility)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/FigSx_reprod_stats.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 16000
#'   input:
#'     - common_genes_rds: '`sm expand(config["DATADIR"] + "/GTEx_v8/fraser2_improvements/tissue_reproducibility_FDR/{method}/reproducibility_genes_rareSpliceAI.Rds", method=["LeafcutterMD", "SPOT", "FRASER", "FRASER2"])`'
#'     - gtex_anno: '`sm config["gtex_sample_anno"]`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_reprod_stats.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_reprod_stats.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in files
input_files <- snakemake@input$common_genes_rds
info_dt <- rbindlist(lapply(input_files, function(file){
    genes <- readRDS(file)
    method_name <- basename(dirname(file))
    return(data.table(method=method_name, testedGenes=length(genes)))
    }))
info_dt[, method:=factor(method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER2"))]
info_dt

#+ plot tested nr of genes
gg_genes <- ggplot(info_dt, aes(method, testedGenes, fill=method)) +
    geom_bar(stat="identity") + 
    labs(x="", y=paste0("Genes present in ≥ ", snakemake@config$required_tissues_per_gene, " tissues")) +
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
    theme_pubr() +
    theme(legend.position="none", 
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) +
    cowplot::background_grid(major="y", minor="y") +
    geom_text(aes(label = testedGenes), nudge_y = 750, 
              size=font_size/.pt)
# gg_genes

#+ plot tested nr of individuals
gtex_v8_anno <- fread(snakemake@input$gtex_anno)
all_tissues <- gtex_v8_anno[,uniqueN(INDIVIDUAL_ID),by="TISSUE"][V1 > 100, sort(TISSUE)]
gtex_v8_anno <- gtex_v8_anno[TISSUE %in% all_tissues,]
required_tissues_per_individual <- snakemake@config$required_tissues_per_individual 

plot_dt <- gtex_v8_anno[,.(tissuesPerIndividual=uniqueN(RNA_ID) ), by="INDIVIDUAL_ID"]
gg_indiv <- ggplot(plot_dt, aes(tissuesPerIndividual)) +
    geom_histogram(binwidth=1) +
    geom_vline(xintercept=required_tissues_per_individual, col="firebrick") +
    labs(x="Available tissues per individual", y="Count") +
    theme_pubr() +
    theme(legend.position="none", 
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) +
    cowplot::background_grid(major="y", minor="y")
# gg_indiv

#+ arrange into figure
gg_sup <- ggarrange(gg_genes,
                   gg_indiv,
                   labels=LETTERS[1:2],
                   nrow=1, ncol=2,
                   align="hv")
# gg_sup

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300) 
