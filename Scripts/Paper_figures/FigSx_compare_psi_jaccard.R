#'---
#' title: Figure Comparison Jaccard Index to psi3/5 
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_comp_psi_jaccard.Rds"`'
#'   threads: 4
#'   resources:
#'     - mem_mb: 32000
#'   input:
#'     - res_fraser1_anno: '`sm expand(config["DATADIR"] + "/GTEx_v8/fraser2_improvements/{dataset}/minK20_25_minN10/PCA__pc0.1/delta0/optQ/fraser1_res_with_jaccard.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'   output:
#'     - res_fraser1_anno_comb: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/delta0/fraser1_res_with_jaccard_comb.tsv"`'
#'     - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_comp_psi_jaccard.png"`'
#'     - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_comp_psi_jaccard.pdf"`'
#'     - outPng_delta: '`sm config["PAPER_FIGDIR"] + "/FigSx_comp_delta_psi_jaccard.png"`'
#'     - outPdf_delta: '`sm config["PAPER_FIGDIR"] + "/FigSx_comp_delta_psi_jaccard.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(BiocParallel)
library(forcats)
library(cowplot)
source("src/R/ggplot_theme_for_manuscript.R")
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", 
                                              scales::scientific_format()(x))) }

#+ read in fraser1 results with jaccard information
res_comb <- rbindlist(bplapply(snakemake@input$res_fraser1_anno, FUN=function(res_file){
    res <- fread(res_file)
    t <- basename(dirname(dirname(dirname(dirname(res_file)))))
    res[, tissue := t]
    return(res)
}) )

#+ save combined fraser1 results
fwrite(res_comb, file=snakemake@output$res_fraser1_anno_comb, sep="\t")

#+ psi5 vs jaccard scatterhist with 2d density across all tissue, fig.width=6, fig.height=4
p_ls <- list() 
for(psiType in c("psi5", "psi3")){
    p_f2outlier <- ggplot(res_comb[type == psiType & FRASER2_outlier == TRUE], 
                          aes(psiValue, jaccardValue)) +
        geom_point(size=0.5) + #, alpha=0.25
        geom_hex(binwidth=0.03) +
        # geom_bin2d(binwidth=0.03) +
        # geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        # geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100))),
            breaks = c(1, 500, 1000, 1500, 2000),
            labels = c(1, 500, 1000, 1500, 2000),
            limits = c(1, 2000)
            ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(psi[5])), expression(bold(psi[3]))), 
             y="Intron Jaccard Index",
             fill="FRASER\nsplicing\noutliers") + 
        theme_manuscript(fig_font_size=font_size, fig_font=font) + 
        theme(legend.position="bottom", 
              legend.direction="horizontal" # "vertical"
        ) + 
        cowplot::background_grid(major="xy", minor="xy") #+
    p_legend1 <- get_legend(p_f2outlier)
    p_f2outlier <- ggExtra::ggMarginal(p_f2outlier + theme(legend.position="none"), 
                                  type = "histogram") 
    p_ls[[paste0("f2outlier_", psiType)]] <- p_f2outlier
    p_ls[[paste0("legend_f2outlier_", psiType)]] <- p_legend1
    
    p_noF2outlier <- ggplot(res_comb[type == psiType & FRASER2_outlier == FALSE], 
                          aes(psiValue, jaccardValue)) +
        geom_point(size=0.5) + #, alpha=0.25
        geom_hex(binwidth=0.03) +
        # geom_bin2d(binwidth=0.03) +
        # geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        # geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100))),
            breaks = c(1, 5e4, 10e4, 15e4, 2e5),
            labels = function(x) ifelse(!x %in% c(1, 1e5, 2e5), "", ifelse(x == 1, 1, scientific_10(x))),
            limits= c(1, 2.3e5)
        ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(psi[5])), expression(bold(psi[3]))), 
             y="Intron Jaccard Index",
             fill="FRASER\nsplicing\noutliers") + 
        theme_pubr() + 
        theme(legend.position="bottom", 
              legend.direction="horizontal", # "vertical",
              # plot.margin=unit(c(1,1,5,1), unit="cm"),
              axis.title=element_text(face="bold"),
              text=element_text(size=font_size)
        ) + 
        cowplot::background_grid(major="xy", minor="xy")
    p_legend2 <- get_legend(p_noF2outlier)
    p_noF2outlier <- ggExtra::ggMarginal(p_noF2outlier + theme(legend.position="none"), 
                                       type = "histogram") 
    p_ls[[paste0("noF2outlier_", psiType)]] <- p_noF2outlier
    p_ls[[paste0("legend_noF2outlier_", psiType)]] <- p_legend2
    
    # scatter plot delta psi vs delta jaccard
    p_f2outlier_deltaScatter <- ggplot(res_comb[type == psiType & FRASER2_outlier == TRUE], 
                          aes(deltaPsi, deltaJaccard)) +
        geom_point(size=0.5) + #, alpha=0.25
        geom_hex(binwidth=0.03) +
        # geom_bin2d(binwidth=0.03) +
        geom_abline(slope=1, xintercept=0, linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100)))
        ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(Delta~psi[5])), expression(bold(Delta~psi[3]))), 
             y=expression(bold(Delta(Intron~Jaccard~Index))),
             fill="FRASER\nsplicing\noutliers") + 
        theme_pubr() + 
        theme(legend.position="bottom", 
              legend.direction="horizontal", # "vertical",
              # plot.margin=unit(c(1,1,5,1), unit="cm"),
              axis.title=element_text(face="bold"),
              text=element_text(size=font_size)
        ) + 
        cowplot::background_grid(major="xy", minor="xy")
    p_legend3 <- get_legend(p_f2outlier_deltaScatter)
    p_f2outlier_deltaScatter <- ggExtra::ggMarginal(p_f2outlier_deltaScatter + theme(legend.position="none"), 
                                       type = "histogram") 
    p_ls[[paste0("f2outlier_deltaScatter_", psiType)]] <- p_f2outlier_deltaScatter
    p_ls[[paste0("legend_f2outlier_deltaScatter_", psiType)]] <- p_legend3
    
    p_noF2outlier_deltaScatter <- ggplot(res_comb[type == psiType & FRASER2_outlier == FALSE], 
                                       aes(deltaPsi, deltaJaccard)) +
        geom_point(size=0.5) + #, alpha=0.25
        geom_hex(binwidth=0.03) +
        # geom_bin2d(binwidth=0.03) +
        geom_abline(slope=1, xintercept=0, linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100)))
        ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(Delta~psi[5])), expression(bold(Delta~psi[3]))), 
             y=expression(bold(Delta(Intron~Jaccard~Index))),
             fill="FRASER\nsplicing\noutliers") + 
        theme_pubr() + 
        theme(legend.position="bottom", 
              legend.direction="horizontal", # "vertical",
              # plot.margin=unit(c(1,1,5,1), unit="cm"),
              axis.title=element_text(face="bold"),
              text=element_text(size=font_size)
        ) + 
        cowplot::background_grid(major="xy", minor="xy")
    p_legend4 <- get_legend(p_noF2outlier_deltaScatter)
    p_noF2outlier_deltaScatter <- ggExtra::ggMarginal(p_noF2outlier_deltaScatter + theme(legend.position="none"), 
                                                    type = "histogram") 
    p_ls[[paste0("noF2outlier_deltaScatter_", psiType)]] <- p_noF2outlier_deltaScatter
    p_ls[[paste0("legend_noF2outlier_deltaScatter_", psiType)]] <- p_legend4
}

#+ compile into figure
gg_sup <- ggdraw() +
    draw_plot(grid::textGrob("FRASER 2.0 outliers"), x=0.1, y=0.95, width=0.3, height=0.05) +
    draw_plot(grid::textGrob("not FRASER 2.0 outliers"), x=0.6, y=0.95, width=0.3, height=0.05) +
    draw_plot(p_ls[["f2outlier_psi5"]], x = 0., y = 0.525, width=0.5, height=0.425) +
    draw_plot(p_ls[["noF2outlier_psi5"]], x = .5, y = 0.525, width=0.5, height=0.425) +
    draw_plot(p_ls[["f2outlier_psi3"]], x = 0, y = 0.1, width=0.5, height=0.425) +
    draw_plot(p_ls[["noF2outlier_psi3"]], x = .5, y = 0.1, width=0.5, height=0.425) +
    draw_plot(p_ls[["legend_f2outlier_psi5"]], x= 0, y=0, width=0.5, height=0.1) +
    draw_plot(p_ls[["legend_noF2outlier_psi5"]], x= 0.5, y=0, width=0.5, height=0.1) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = 12,
                    x = c(0, 0.5, 0., 0.5), y = c(1, 1, 0.5, 0.5))
gg_sup

# delta scatter plot
gg_sup2 <- ggdraw() +
    draw_plot(grid::textGrob("FRASER 2.0 outliers"), x=0.1, y=0.95, width=0.3, height=0.05) +
    draw_plot(grid::textGrob("not FRASER 2.0 outliers"), x=0.6, y=0.95, width=0.3, height=0.05) +
    draw_plot(p_ls[["f2outlier_deltaScatter_psi5"]], x = 0., y = 0.525, width=0.5, height=0.425) +
    draw_plot(p_ls[["noF2outlier_deltaScatter_psi5"]], x = .5, y = 0.525, width=0.5, height=0.425) +
    draw_plot(p_ls[["f2outlier_deltaScatter_psi3"]], x = 0, y = 0.1, width=0.5, height=0.425) +
    draw_plot(p_ls[["noF2outlier_deltaScatter_psi3"]], x = .5, y = 0.1, width=0.5, height=0.425) +
    draw_plot(p_ls[["legend_f2outlier_psi5"]], x= 0, y=0, width=0.5, height=0.1) +
    draw_plot(p_ls[["legend_noF2outlier_deltaScatter_psi5"]], x= 0.5, y=0, width=0.5, height=0.1) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = 12,
                    x = c(0, 0.5, 0., 0.5), y = c(1, 1, 0.5, 0.5))
gg_sup2

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=page_width, unit=width_unit)
ggsave(plot=gg_sup2, filename=snakemake@output$outPng_delta, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup2, filename=snakemake@output$outPdf_delta, width=page_width, height=0.8*page_width, unit=width_unit)


