#'---
#' title: Paper figure S2 (pseudocount example/motivation)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS2.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - fraser2_fds_pc1: '`sm config["DATADIR"] + "/GTEx_v8/fds/minK20_95_minN1/PCA/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'
#'     - fraser2_fds_pc01: '`sm config["DATADIR"] + "/GTEx_v8/fds/minK20_95_minN1/PCA__pc0.1/savedObjects/Skin_-_Not_Sun_Exposed_Suprapubic__optQ__newFilt/pvaluesBetaBinomial_junction_jaccard.h5"`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_pc.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_pc.pdf"`'
#'   type: script
#'---

# minK20_95_minN1

# #'    - S_outPng: '`sm config["PAPER_FIGDIR"] + "/FigS2.png"`'
# #'    - S_outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS2.pdf"`'

# #'     - res_fraser2_pc1: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK20_5_minN10/PCA/Muscle_-_Skeletal/optQ__newFilt/results_junction.tsv"`'
# #'     - res_fraser2_pc01: '`sm config["DATADIR"] + "/GTEx_v8/FRASER2_results/minK20_5_minN10/PCA__pc0.1/Muscle_-_Skeletal/optQ__newFilt/results_junction.tsv"`'

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(ggplot2)
library(ggpubr)
library(cowplot)
# library(gridExtra)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size - 1
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ read in fds objects
fds_pc1 <- loadFraserDataSet(file=snakemake@input$fraser2_fds_pc1)
# res_pc1 <- fread(snakemake@input$res_fraser2_pc1)
# res_pc01 <- fread(snakemake@input$res_fraser2_pc01)
# res_pc1[, tmp:=paste(sampleID, seqnames, start, end, strand, sep="_")]
# res_pc01[, tmp:=paste(sampleID, seqnames, start, end, strand, sep="_")]

# pc1_out <- which(padjVals(fds_pc1, type="jaccard", level="site") < snakemake@config$f, arr.ind=TRUE)
# rq90 <- rowQuantiles(N(fds_pc1, type="jaccard"), prob=0.9)
# idx_fp <- pc1_out[which.min(rq90[pc1_out[,1]]), 1]
# pc_example_fp <- plotExpression(fds_pc1, type="jaccard", idx=idx_fp, padjCutoff=0.1)
# pc_example_fp_predObs <- plotExpectedVsObservedPsi(fds_pc1, type="jaccard", idx=idx_fp, padjCutoff=0.1)

# jidx <- 1271 # 1576 # which filtering?
jidx <- 1529 # for qauntile=25% filtering
# jidx <- 1550 # for quantile=5% filtering
pc_example <- plotExpression(fds_pc1, type="jaccard", idx=jidx, padjCutoff=0.1)
# pc_example <- plotExpression(fds_pc1, result=res_pc1[! tmp %in% res_pc01$tmp][order(meanTotalCounts)][1,], padjCutoff=0.1)
pc_example <- pc_example + ggtitle("") +
    theme_pubr() +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# pc_example
predObs_example_pc1 <- plotExpectedVsObservedPsi(fds_pc1, type="jaccard", idx=jidx, padjCutoff=0.1)
predObs_example_pc1 <- predObs_example_pc1 +
    labs(title="Pseudocount = 1") +
    theme_pubr() +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          plot.title = element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# predObs_example_pc1 <- plotExpectedVsObservedPsi(fds_pc1, result=res_pc1[! tmp %in% res_pc01$tmp][order(meanTotalCounts)][1,], padjCutoff=0.1)
# predObs_example_pc1

fds_pc01 <- loadFraserDataSet(file=snakemake@input$fraser2_fds_pc01)
predObs_example_pc01 <- plotExpectedVsObservedPsi(fds_pc01, type="jaccard", idx=jidx, padjCutoff=0.1)
predObs_example_pc01 <- predObs_example_pc01 +
    labs(title="Pseudocount = 0.1") +
    theme_pubr() +
    theme(legend.position="none",
          axis.title=element_text(face="bold"),
          plot.title = element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# predObs_example_pc01

#+ combine panels into figure, width=22, height=15
gg_sup <- ggarrange(pc_example,
                       predObs_example_pc1,
                       predObs_example_pc01,
                       labels=letters[1:3],
                       nrow=1, ncol=3)
# gg_sup

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=0.5*page_width, unit=width_unit, dpi=300) 
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=0.5*page_width, unit=width_unit, dpi=300) 
