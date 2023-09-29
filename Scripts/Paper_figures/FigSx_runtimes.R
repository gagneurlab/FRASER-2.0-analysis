#'---
#' title: Paper figure Sx (runtimes of FRASER 2.0 steps)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figSx_runtimes.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 24000
#'   input:
#'    - timings_fraser2: '`sm expand(config["DATADIR"] + "/GTEx_v8/runtimes/fraser_v2/PCA__pc0.1/timings/{dataset}_runtime.Rds", dataset=["Liver", "Skin_-_Not_Sun_Exposed_Suprapubic", "Heart_-_Left_Ventricle"])`'
#'    - timings_fraser2_oldFilt: '`sm expand(config["DATADIR"] + "/GTEx_v8/runtimes/fraser_v2_oldFilt/PCA__pc0.1/timings/{dataset}_runtime.Rds", dataset=["Liver", "Skin_-_Not_Sun_Exposed_Suprapubic", "Heart_-_Left_Ventricle"])`'
#'    - timings_fraser1: '`sm expand(config["DATADIR"] + "/GTEx_v8/runtimes/fraser_v1/PCA/timings/{dataset}_runtime.Rds", dataset=["Liver", "Skin_-_Not_Sun_Exposed_Suprapubic", "Heart_-_Left_Ventricle"])`'
#'   output:
#'    - outPng:  '`sm config["PAPER_FIGDIR"] + "/figSx_runtimes.png"`'
#'    - outPdf:  '`sm config["PAPER_FIGDIR"] + "/figSx_runtimes.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# slurm jobs corresponding to the results:
# squeue --me
# 7606264 noninterr snakejob scheller  R       0:38      1 ouga08
# 7606265 noninterr snakejob scheller  R       0:38      1 ouga08
# 7606266 noninterr snakejob scheller  R    1:06:27      1 ouga08
# 7606267 noninterr snakejob scheller  R       9:14      1 ouga08

# 7606267 noninterr snakejob scheller PD       0:00      1 (Priority)
# 7606268 noninterr snakejob scheller PD       0:00      1 (Priority)
# 7606269 noninterr snakejob scheller PD       0:00      1 (Priority)
# 7606270 noninterr snakejob scheller PD       0:00      1 (Priority)
# 7606271 noninterr snakejob scheller PD       0:00      1 (Priority)
# 7606272 noninterr snakejob scheller PD       0:00      1 (Priority)


# library(FRASER)
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

#+ create data table from time points
runtime_dt <- rbindlist(lapply(c(snakemake@input$timings_fraser2, snakemake@input$timings_fraser2_oldFilt, snakemake@input$timings_fraser1), function(timings_file){
    timings <- readRDS(timings_file)
    t_name <- timings[["tissue"]]# gsub("_runtime.Rds", "", basename(timings_file))
    t_name <- gsub("_", " ", gsub("_-_", " ", t_name))
    dt <- data.table(tissue=t_name, N=timings[["sample_size"]], introns=timings[["nr_introns"]], genes=timings[["nr_genes"]], method=timings[["method"]],
                     step=c("Filtering", "Hyperparameter\noptimization", "Autoencoder fit", "P value\ncalculation", "Results\nintron-level", "Results\ngene-level", "Results", "Full analysis"),
                     time_sec=c(timings[["end_filtering"]] - timings[["begin_filtering"]], 
                                timings[["end_hyperpar"]] - timings[["begin_hyperpar"]],
                                timings[["end_fit"]] - timings[["begin_fit"]],
                                timings[["end_pvals"]] - timings[["begin_pvals"]],
                                timings[["end_res_intronlevel"]] - timings[["begin_res_intronlevel"]],
                                timings[["end_res_genelevel"]] - timings[["begin_res_genelevel"]],
                                timings[["end_results"]] - timings[["begin_results"]],
                                timings[["end"]] - timings[["begin"]] 
                     ))
}))
runtime_dt[, time_min := time_sec / 60]
runtime_dt[, time_h := time_min / 60]

runtime_dt <- runtime_dt[!step %in% c("Results\nintron-level", "Results\ngene-level"),]
runtime_dt[, step:=factor(step, levels=c("Full analysis", "Filtering", "Hyperparameter\noptimization", "Autoencoder fit", "P value\ncalculation", "Results"))]
runtime_dt[, tissue_long:=factor(paste0(tissue, " (N=", N, ")"), levels=runtime_dt[!duplicated(runtime_dt, by=c("tissue", "N")),paste0(tissue, " (N=", N, ")")])]

# only show default fraser 2.0 
runtime_dt <- runtime_dt[method %in% c("FRASER", "FRASER 2.0 (recommmended filtering settings for FRASER 2.0)")]
runtime_dt[method == "FRASER 2.0 (recommmended filtering settings for FRASER 2.0)", method := "FRASER 2.0"]

#+ visualize runtimes
# g_sup <- ggplot(runtime_dt, aes(step, time_h, fill=method)) + 
#     facet_wrap(~tissue_long) +
#     geom_col(position="dodge") +
#     scale_fill_manual(values=c("FRASER"="dodgerblue3", "FRASER 2.0"="purple4")) +
#     labs(x="Computational step", y="Runtime (h)") +
#     theme_manuscript(fig_font_size=font_size, fig_font=font) + 
#     theme(legend.title=element_blank(),
#           axis.text.x=element_text(angle=45, hjust=1)) + 
#     cowplot::background_grid(major="y", minor="y") 
# g_sup

g_sup <- ggplot(runtime_dt, aes(N, time_h, color=method)) + 
    facet_wrap(~step, scales="free_y") +
    # facet_wrap(~step) +
    geom_point() +
    geom_line() +
    scale_color_manual(values=c("FRASER"="dodgerblue3", 
                                "FRASER 2.0"="purple4")) +
    scale_x_continuous(limits=c(150, 650)) +
    labs(x="Sample size", y="Runtime (h)") +
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    theme(legend.title=element_blank()) + 
    cowplot::background_grid(major="xy", minor="xy") 
g_sup

#+ save figure as png and pdf
ggsave(plot=g_sup, filename=snakemake@output$outPng, width=page_width, height=0.66*page_width, unit=width_unit, dpi=300) 
ggsave(plot=g_sup, filename=snakemake@output$outPdf, width=page_width, height=0.66*page_width, unit=width_unit, dpi=300) 
