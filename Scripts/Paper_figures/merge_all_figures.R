#'---
#' title: Merge figures
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/merge_all.Rds"`'
#'   threads: 1
#'   resources:
#'    - mem_mb: 16000
#'   input:
#'    - Figure1:   '`sm config["PAPER_FIGDIR"] + "/Fig1_full_v4.pdf"`'
#'    - Figure2:   '`sm config["PAPER_FIGDIR"] + "/Fig2.pdf"`'
#'    - Figure3:   '`sm config["PAPER_FIGDIR"] + "/Fig3.pdf"`'
#'    - Figure4:   '`sm config["PAPER_FIGDIR"] + "/Fig4.pdf"`'
#'    - FigureS01:   '`sm config["PAPER_FIGDIR"] + "/FigSx_recall_jaccard_vs_psi.pdf"`'
#'    - FigureS02:   '`sm config["PAPER_FIGDIR"] + "/FigSx_paramOpt_recallAt20.pdf"`'
#'    - FigureS03:   '`sm config["PAPER_FIGDIR"] + "/FigSx_filterDeltaOpt_recallAt20.pdf"`'
#'    - FigureS04:   '`sm config["PAPER_FIGDIR"] + "/FigSx_rhoOpt_recallAt20.pdf"`'
#'    - FigureS05:   '`sm config["PAPER_FIGDIR"] + "/heatmaps/FigSx_heatmap_full.pdf"`'
#'    - FigureS06:   '`sm config["PAPER_FIGDIR"] + "/FigSx_comp_psi_jaccard.pdf"`'
#'    - FigureS07:   '`sm config["PAPER_FIGDIR"] + "/FigSx_qqplots.pdf"`'
#'    - FigureS08:   '`sm config["PAPER_FIGDIR"] + "/FigSx_precision_recall_FDRsignif_allTissues.pdf"`'
#'    - FigureS09:   '`sm config["PAPER_FIGDIR"] + "/FigSx_gtex_venn_recall.pdf"`'
#'    - FigureS10:   '`sm config["PAPER_FIGDIR"] + "/FigSx_reproducibility.pdf"`'
#'    - FigureS11:   '`sm config["PAPER_FIGDIR"] + "/FigSx_testedGenesIntrons.pdf"`'
#'    - FigureS12:   '`sm config["PAPER_FIGDIR"] + "/FigSx_power_analysis.pdf"`'
#'    - numbers_for_manuscript: '`sm config["htmlOutputPath"] + "/paper/paper_numbers.html"`'
#'   output:
#'   - figure:  '`sm config["PAPER_FIGDIR"] + "/Figure_all_main.pdf"`'
#'   - figureS: '`sm config["PAPER_FIGDIR"] + "/Figure_all_sup.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ get figure files
mainFigures <- unlist(unique(snakemake@input[grepl("Figure[0-9]",  names(snakemake@input))]))
supFigures  <- unlist(unique(snakemake@input[grepl("FigureS[0-9]", names(snakemake@input))]))

message(mainFigures)

message(date(), ': Start with main figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figure, ' ',
    paste(mainFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=',
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figure), ' ', snakemake@output$figure))



message(supFigures)


message(date(), ': Start with supplement figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figureS, ' ',
    paste(supFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=',
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figureS), ' ', snakemake@output$figureS))
