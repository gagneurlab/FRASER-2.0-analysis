#'---
#' title: Figure Comparison Jaccard Index to psi3/5 
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS1_comp_psi_jaccard.Rds"`'
#'   threads: 4
#'   resources:
#'     - mem_mb: 32000
#'   input:
#'     - res_fraser1_anno: '`sm expand(config["DATADIR"] + "/GTEx_v8/fraser2_improvements/{dataset}/minK20_25_minN10/PCA__pc0.1/delta0/optQ/fraser1_res_with_jaccard.tsv", dataset=config["tissues_for_reproducibility"], allow_missing=True)`'
#'   output:
#'     - res_fraser1_anno_comb: '`sm config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/delta0/fraser1_res_with_jaccard_comb.tsv"`'
#'     - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS1_comp_psi_jaccard.png"`'
#'     - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS1_comp_psi_jaccard.pdf"`'
#'     - outPng_delta: '`sm config["PAPER_FIGDIR"] + "/FigS1_comp_delta_psi_jaccard.png"`'
#'     - outPdf_delta: '`sm config["PAPER_FIGDIR"] + "/FigS1_comp_delta_psi_jaccard.pdf"`'
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
register(MulticoreParam(snakemake@threads))

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

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
        geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100)))
            ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(psi[5])), expression(bold(psi[3]))), 
             y="Intron Jaccard Index",
             fill="FRASER\nsplicing\noutliers") + 
        theme_pubr() + 
        theme(legend.position="bottom", 
              legend.direction="vertical",
              # plot.margin=unit(c(1,1,5,1), unit="cm"),
              axis.title=element_text(face="bold"),
              text=element_text(size=font_size)
        ) + 
        cowplot::background_grid(major="xy", minor="xy") #+
        # annotation_custom(grid::textGrob("FRASER2 outlier", rot = 90), 
        #                   xmin = -1, xmax = -1, ymin = 0)
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
        geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
        scale_fill_gradientn(
            # trans = "log10",
            colors = rev(brewer.pal(9, "YlGnBu")),
            values = c(0, exp(seq(-5, 0, length.out = 100)))
        ) +
        labs(x=ifelse(psiType == "psi5", expression(bold(psi[5])), expression(bold(psi[3]))), 
             y="Intron Jaccard Index",
             fill="FRASER\nsplicing\noutliers") + 
        theme_pubr() + 
        theme(legend.position="bottom", 
              legend.direction="vertical",
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
        # geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        # geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
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
              legend.direction="vertical",
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
        # geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
        # geom_vline(xintercept=0.75, col="firebrick", linetype="dashed") +
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
              legend.direction="vertical",
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

#+ create psi and jaccard decile bins
res_comb[, psiBin := cut(psiValue, breaks=c(-0.1, seq(0.1, 1, by=0.1)))]
res_comb[, psiBin := fct_recode(psiBin, "[0,0.1]" = "(-0.1,0.1]")]
# res_comb[, table(psiBin)]
# res_comb[is.na(psiBin)]
res_comb[, jaccardBin := cut(jaccardValue, breaks=c(-0.1, seq(0.1, 1, by=0.1)))]
res_comb[, jaccardBin := fct_recode(jaccardBin, "[0,0.1]" = "(-0.1,0.1]")]
# res_comb[, table(jaccardBin)]
# res_comb[!is.na(jidx_f2), table(is.na(jaccardBin))]

#+ barplots of table of psi deciles vs low jaccard value, fig.width=8, fig.height=4
plot_dt <- as.data.table(table(res_comb[, .(psiBin, jaccardBin == "[0,0.1]", type)]))
plot_dt[, V2 := factor(V2, levels=c("TRUE", "FALSE"))]
plot_dt[, psiBin := fct_relevel(psiBin, "[0,0.1]")]
plot_dt[, levels(psiBin)]
plot_dt[, prop:=N/sum(N), by="psiBin,type"]
plot_dt
# p_bar <- ggplot(as.data.frame(table(res_comb[, .(psiBin, jaccardBin == "[0,0.1]", type)])), aes(x=psiBin, y = Freq, fill=V2)) +
p_bar <- ggplot(plot_dt, aes(x=psiBin, y = prop, fill=V2)) +
    facet_wrap(~type) +
    geom_bar(stat="identity") +
    labs(x=expression(bold(psi["5/3"] ~ "deciles")),
         # title="All GTEx tissues",
         y="" ) +
    guides(fill=guide_legend(title="jaccard <= 0.1")) +
    scale_fill_brewer(palette="Paired") +
    scale_y_continuous(labels = scales::percent) +
    theme_pubr() + 
    theme(legend.position="bottom", 
          legend.spacing.x=unit(0.5, "cm"),
          axis.title=element_text(face="bold"),
          axis.text.x=element_text(angle=45, hjust=0.5, vjust=0.5),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# p_bar

#+ histogram of psi5/3/jaccard values of fraser outliers
melt_dt <- melt(res_comb[, .(sampleID, seqnames,start,end,strand,type,psiValue, jaccardValue,FRASER2_outlier)], id.vars=c("sampleID", "seqnames","start","end","strand","type", "FRASER2_outlier"),
                value.name="metricValue", variable.name="metric")
melt_dt[type == "psi5" & metric == "psiValue", metric := "psi5"]
melt_dt[type == "psi3" & metric == "psiValue", metric := "psi3"]
melt_dt[metric == "jaccardValue", metric := "Intron Jaccard Index"]
p_hist <- ggplot(melt_dt, aes(metricValue, fill=metric)) + 
    facet_grid(type~FRASER2_outlier, labeller=labeller(FRASER2_outlier=label_both)) +
    geom_histogram(alpha=0.7, aes(y = ..density..), position="identity") +
    # geom_histogram(alpha=0.5) +
    labs(x="Splice metric value", fill = "Splice metric") +
    theme_pubr() + 
    theme(legend.position="bottom", 
          axis.title=element_text(face="bold"),
          text=element_text(size=font_size)) + 
    cowplot::background_grid(major="xy", minor="xy")
# p_hist

#+ compile into figure
# gg_sup <- ggarrange(print(p_psi5), 
#                     print(p_psi3),
#                     labels=LETTERS[1:2],
#                     nrow=1, ncol=2,
#                     common.legend=TRUE)
# gg_sup <- patchwork::wrap_plots(p_psi5, p_psi3, nrow = 1, guides="collect", )
library("cowplot")
# gg_sup <- ggdraw() +
#     draw_plot(p_psi5, x = 0, y = 0, width=0.45) +
#     draw_plot(p_psi3, x = .45, y = 0, width=0.45) +
#     draw_plot(p_legend, x= .9, y=0, width=0.1) +
#     draw_plot_label(label = c("A", "B"), size = 15,
#                     x = c(0, 0.5), y = c(1, 1))
# gg_sup <- ggdraw() +
#     draw_plot(p_psi5, x = 0, y = .5, width=0.475, height=0.5) +
#     draw_plot(p_psi3, x = .475, y = .5, width=0.475, height=0.5) +
#     draw_plot(p_legend, x= .9, y= .5, width=0.1, height=0.5) +
#     draw_plot(p_bar, x = 0, y = 0, width=1, height=0.5) +
#     draw_plot_label(label = c("A", "B", "C"), size = 15,
#                     x = c(0, 0.5, 0), y = c(1, 1, 0.5))
# gg_sup <- ggdraw() +
#     draw_plot(p_ls[[3]], x = 0, y = 0.5, width=0.45, height=0.5) +
#     draw_plot(p_ls[[7]], x = .45, y = 0.5, width=0.45, height=0.5) +
#     draw_plot(p_ls[[1]], x = 0, y = 0, width=0.45, height=0.5) +
#     draw_plot(p_ls[[5]], x = .45, y = 0, width=0.45, height=0.5) +
#     draw_plot(p_ls[[4]], x= .9, y=0.5, width=0.1, height=0.5) +
#     draw_plot(p_ls[[2]], x= .9, y=0, width=0.1, height=0.5) +
#     draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
#                     x = c(0, 0.45, 0, 0.45), y = c(1, 1, 0.5, 0.5))
# gg_sup

gg_sup <- ggdraw() +
    draw_plot(p_ls[["f2outlier_psi5"]], x = 0.05, y = 0.5, width=0.425, height=0.5) +
    draw_plot(p_ls[["f2outlier_psi3"]], x = .475, y = 0.5, width=0.425, height=0.5) +
    draw_plot(p_ls[["noF2outlier_psi5"]], x = 0.05, y = 0, width=0.425, height=0.5) +
    draw_plot(p_ls[["noF2outlier_psi5"]], x = .475, y = 0, width=0.425, height=0.5) +
    draw_plot(p_ls[["legend_noF2outlier_psi5"]], x= .9, y=0.5, width=0.1, height=0.5) +
    draw_plot(p_ls[["legend_noF2outlier_psi3"]], x= .9, y=0, width=0.1, height=0.5) +
    draw_plot(grid::textGrob("FRASER 2.0 outliers", rot=90), x=0, y=0.6, width=0.05, height=0.3) +
    draw_plot(grid::textGrob("not FRASER 2.0 outliers", rot=90), x=0, y=0.1, width=0.05, height=0.3) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                    x = c(0.05, 0.45, 0.05, 0.45), y = c(1, 1, 0.5, 0.5))
gg_sup

# delta scatter plot
gg_sup2 <- ggdraw() +
    draw_plot(p_ls[["f2outlier_deltaScatter_psi5"]], x = 0.05, y = 0.5, width=0.425, height=0.5) +
    draw_plot(p_ls[["f2outlier_deltaScatter_psi3"]], x = .475, y = 0.5, width=0.425, height=0.5) +
    draw_plot(p_ls[["noF2outlier_deltaScatter_psi5"]], x = 0.05, y = 0, width=0.425, height=0.5) +
    draw_plot(p_ls[["noF2outlier_deltaScatter_psi3"]], x = .475, y = 0, width=0.425, height=0.5) +
    draw_plot(p_ls[["legend_noF2outlier_deltaScatter_psi5"]], x= .9, y=0.5, width=0.1, height=0.5) +
    draw_plot(p_ls[["legend_noF2outlier_deltaScatter_psi3"]], x= .9, y=0, width=0.1, height=0.5) +
    draw_plot(grid::textGrob("FRASER 2.0 outliers", rot=90), x=0, y=0.6, width=0.05, height=0.3) +
    draw_plot(grid::textGrob("not FRASER 2.0 outliers", rot=90), x=0, y=0.1, width=0.05, height=0.3) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                    x = c(0.05, 0.45, 0.05, 0.45), y = c(1, 1, 0.5, 0.5))
gg_sup2

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300)
# ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300)
# ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=0.4*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup2, filename=snakemake@output$outPng_delta, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup2, filename=snakemake@output$outPdf_delta, width=page_width, height=0.8*page_width, unit=width_unit, dpi=300)

# #+ different options for plotting, fig.width=6, fig.height=4
# p_psi5 <- ggplot(res_f1[type == "psi5"], aes(psiValue, jaccardValue)) + 
#     ggtitle(dataset) + labs(x="psi_5", y="Intron Jaccard Index") 
# # option stat_density2d
# p_psi5 + stat_density2d(h = 0.1, bins = 100, aes( fill = ..level..), geom = "polygon") + 
#     colorscale +
#     geom_abline(slope=1, intercept=0, col="black") +
#     geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
#     geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")
# # option 2 geom_hex
# p_psi5 + geom_hex() +  
#     colorscale +
#     geom_abline(slope=1, intercept=0, col="black") +
#     geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
#     geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")
# # option 3 geom_bin2d
# p_psi5 + geom_bin2d() + 
#     colorscale +
#     geom_abline(slope=1, intercept=0, col="black") +
#     geom_hline(yintercept=0.25, col="firebrick", linetype="dashed") +
#     geom_vline(xintercept=0.75, col="firebrick", linetype="dashed")

# #+ boxplot of jaccard values for fraser1 outliers with large psi3/5 value, fig.width=8, fig.heigth=4
# ggplot(res_comb[psiValue >= 0.9], aes(jaccardValue)) +
#     facet_wrap(~type) +
#     geom_boxplot() + 
#     labs(x="Intron Jaccard Index of FRASER1 outliers\nwith psi3/5 >= 0.9",
#          title="All GTEx tissues")
# ggplot(res_comb[psiValue >= 0.75], aes(jaccardValue)) +
#     facet_wrap(~type) +
#     geom_boxplot() + 
#     labs(x="Intron Jaccard Index of FRASER1 outliers\nwith psi3/5 >= 0.75",
#          title="All GTEx tissues")


# #+ table of high psi value vs low jaccard value
# res_comb[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.9` = psiValue >= 0.9, `jaccard <= 0.1` = jaccardValue <= 0.1)]
# res_comb[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.9` = psiValue >= 0.9, `jaccard <= 0.1` = jaccardValue <= 0.1)]
# tb_5 <- table(res_comb[, .(`psi5 >= 0.9`, `jaccard <= 0.1`)] )
# tb_3 <- table(res_comb[, .(`psi3 >= 0.9`, `jaccard <= 0.1`)] )
# tb_5[2,2] / sum(tb_5)
# tb_3[2,2] / sum(tb_3)
# res_comb[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.85` = psiValue >= 0.85, `jaccard <= 0.15` = jaccardValue <= 0.15)]
# res_comb[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.85` = psiValue >= 0.85, `jaccard <= 0.15` = jaccardValue <= 0.15)]
# tb_5 <- table(res_comb[, .(`psi5 >= 0.85`, `jaccard <= 0.15`)] )
# tb_3 <- table(res_comb[, .(`psi3 >= 0.85`, `jaccard <= 0.15`)] )
# tb_5[2,2] / sum(tb_5)
# tb_3[2,2] / sum(tb_3)
# res_comb[type=="psi5" & !is.na(jidx_f2), `:=`(`psi5 >= 0.75` = psiValue >= 0.75, `jaccard <= 0.25` = jaccardValue <= 0.25)]
# res_comb[type=="psi3" & !is.na(jidx_f2), `:=`(`psi3 >= 0.75` = psiValue >= 0.75, `jaccard <= 0.25` = jaccardValue <= 0.25)]
# tb_5 <- table(res_comb[, .(`psi5 >= 0.75`, `jaccard <= 0.25`)] )
# tb_3 <- table(res_comb[, .(`psi3 >= 0.75`, `jaccard <= 0.25`)] )
# tb_5[2,2] / sum(tb_5)
# tb_3[2,2] / sum(tb_3)
# 
# #+ barplots of table of high psi value vs low jaccard value, fig.width=6, fig.heigth=4
# res_comb[!is.na(jidx_f2), `:=`(psiBig = psiValue >= 0.9, jaccardSmall = jaccardValue <= 0.1)]
# ggplot(as.data.frame(table(res_comb[, .(psiBig, jaccardSmall, type)])), aes(x=psiBig, y = Freq, fill=jaccardSmall)) +
#     facet_wrap(~type) +
#     geom_bar(stat="identity") +
#     labs(x="psi >= 0.9",
#          title="All GTEx tissues") +
#     guides(fill=guide_legend(title="jaccard <= 0.1"))
# res_comb[!is.na(jidx_f2), `:=`(psiBig = psiValue >= 0.75, jaccardSmall = jaccardValue <= 0.25)]
# ggplot(as.data.frame(table(res_comb[, .(psiBig, jaccardSmall, type)])), aes(x=psiBig, y = Freq, fill=jaccardSmall)) +
#     facet_wrap(~type) +
#     geom_bar(stat="identity") +
#     labs(x="psi >= 0.75",
#          title="All GTEx tissues") +
#     guides(fill=guide_legend(title="jaccard <= 0.25"))

