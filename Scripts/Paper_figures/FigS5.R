#'---
#' title: Paper figure S5 (GTEx QQ plots)
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/figS5.Rds"`'
#'   threads: 1
#'   resources:
#'     - mem_mb: 12000
#'   input:
#'     - global_qq: '`sm expand(config["DATADIR"] + "/GTEx_v8/fraser2_improvements/minK20_25_minN10/optQ/PCA__pc0.1/{dataset}/global_qqPlots.Rds", dataset=config["tissues_for_detailed_analysis"])`'
#'   output:
#'    - outPng: '`sm config["PAPER_FIGDIR"] + "/FigS5.png"`'
#'    - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigS5.pdf"`'
#'   type: script
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE 
library(ggplot2)
library(ggpubr)
library(cowplot)

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit
point_size <- 0.5
maxLineNchar <- 20

#+ read in plots and data for the different panels
input_files <- snakemake@input$global_qq
if(length(input_files) > 15){
    input_files <- input_files[!grepl("Skin_-_Not_Sun_Exposed_Suprapubic", input_files)]
}
gg_ls <- lapply(sort(sample(input_files, 15)), FUN=function(qq_plot_rds){
    global_qq_plot <- readRDS(qq_plot_rds)
    tissue <- basename(dirname(qq_plot_rds))
    tissue <- gsub("_", " ", gsub("_-_", " ", tissue))
    if(nchar(tissue) > maxLineNchar){
        fsplit <- strsplit(tissue, " ", fixed=TRUE)[[1]]
        tissue <- fsplit[1]
        nc <- nchar(tissue)
        for(s in fsplit[-1]){
            nc <- nc + nchar(s)
            if(nc > maxLineNchar){
                tissue <- paste(tissue, s, sep="\n")
                nc <- nchar(s)
            } else{
                tissue <- paste(tissue, s, sep=" ")
            }
            
        }
    }
    global_qq_plot <- global_qq_plot + ggtitle(tissue) + 
        labs(
            x = expression(bold(-log[10]("expected"~"P"))),
            y = expression(bold(-log[10]("observed"~"P")))
        ) + 
        guides(col = guide_legend(nrow = 1, title="")) +
        theme_pubr() +
        theme(axis.title=element_text(face="bold"), 
              legend.position="top",
              title=element_text(face="bold"),
              text=element_text(size=font_size),
              axis.text=element_text(size=10),
              legend.text=element_text(size=10)) +
        cowplot::background_grid(major="xy", minor="xy")
    global_qq_plot$layers[[1]]$aes_params$size <- point_size
        
    return(global_qq_plot)
})

#+ combine into one figure
gg_sup <- ggarrange(plotlist=gg_ls, 
                    nrow=5, ncol=3,
                    common.legend=TRUE,
                    legend="bottom",
                    labels=LETTERS[seq_along(gg_ls)]
                    )
gg_sup

#+ save figure as png and pdf
ggsave(plot=gg_sup, filename=snakemake@output$outPng, width=page_width, height=1.6*page_width, unit=width_unit, dpi=300)
ggsave(plot=gg_sup, filename=snakemake@output$outPdf, width=page_width, height=1.6*page_width, unit=width_unit, dpi=300)

