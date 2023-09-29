#'---
#' title: Paper figure S17 reproduibilty
#' author: Christian Mertes
#' wb:
#'  log:
#'   - snakemake: '`sm config["log_dir"] + "/snakemake/paper_figures/FigSx_reproducibility.Rds"`' 
#'  threads: 10
#'  resources:
#'   - mem_mb: 64000
#'  input:
#'   - data:   '`sm config["DATADIR"] + "/GTEx_v8/reproducibility_fraser1paper/rareSpliceSite__0.0_reproducibility.RDS"`'
#'   - table:  '`sm config["DATADIR"] + "/GTEx_v8/reproducibility_fraser1paper/rareSpliceSite__0.0_reproducibility.tsv.gz"`'
#'  output:
#'   - outPng: '`sm config["PAPER_FIGDIR"] + "/FigSx_reproducibility.png"`'
#'   - outPdf: '`sm config["PAPER_FIGDIR"] + "/FigSx_reproducibility.pdf"`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)

#+ echo=FALSE
# source("./src/r/config.R")
.libPaths("~/R/4.1/FRASER2/")
library(FRASER)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggpubr)
source("src/R/ggplot_theme_for_manuscript.R")

#+ read in figure font size and width params from config
font_size <- snakemake@config$font_size
font <- snakemake@config$font
page_width <- snakemake@config$page_width
width_unit <- snakemake@config$width_unit

#+ input
data_file  <- snakemake@input$data
table_file <- snakemake@input$table
AE_NAME    <- "FRASER"
outPng     <- snakemake@output$outPng
outPdf     <- snakemake@output$outPdf


data_file
table_file
AE_NAME
outPng


#' 
#' Read in data
#' 
data <- readRDS(data_file)
res <- fread(table_file)

names(data$datatables) <- c("dt2p", "dt2p9p", "dt2p7p", "dt2p5p")
attach(data$datatables)


#' 
#' rename methods
renameTheMethod <- function(dt){
    dt[, Method := gsub("_p", "", Method)]
    dt[, Method := factor(Method, levels=c("LeafcutterMD", "SPOT", "FRASER", "FRASER2"))]
}
dt2p <- renameTheMethod(dt2p)
dt2p5p <- renameTheMethod(dt2p5p)
dt2p7p <- renameTheMethod(dt2p7p)
dt2p9p <- renameTheMethod(dt2p9p)


#'
#'
#'
gt1 <- ggplot(res[,.(Method, totalTested)], aes(totalTested, fill=Method)) + 
    geom_histogram(position="dodge") + 
    scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
    ylab("Number of tested events") + 
    xlab(bquote("Number of tissues event is tested")) + 
    theme_manuscript(fig_font_size=font_size, fig_font=font) + 
    grids(axis="y") + 
    scale_y_log10()
gt1

# plot total number of events passsing cutoffs
plotTotalNumberOfEvents <- function(dt, value){
    plot_col <- paste0("hits", value)
    ggplot(dt[get(plot_col) != 0], aes(x=.data[[plot_col]], fill=Method)) + 
        geom_bar(position="dodge") + 
        theme_bw() + 
        facet_wrap(~spliceVariant) + 
        scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
        ylab("Number of events") + 
        xlab(bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-.(value) ~ ")")) + 
        theme_manuscript(fig_font_size=font_size, fig_font=font) + 
        grids(axis="y") + 
        scale_y_log10()
}

plotTotalNumberOfEventsByCategory <- function(dt, value){
    plot_col <- paste0("hits", value)
    dt_melt <- melt(dt, id.vars=c("subjectID", "geneID", "Method", "totalTested", "IMPACT", "spliceVariant"))
    dt_melt <- dt_melt[value > 0]
    # combine tissues into categories 1, 2-5, > 5
    dt_melt[value == 1, category:="1"]
    dt_melt[value >= 2 & value <= 5, category:="2-5"]
    dt_melt[value >= 6 , category:="> 5"]
    dt_melt <- dt_melt[, .N, by="Method,spliceVariant,variable,category"]
    dt_melt[, category:=factor(category, levels=c("1", "2-5", "> 5"))]
    
    ggplot(dt_melt[variable == plot_col], aes(x=category, y=N, fill=Method)) + 
        geom_bar(position="dodge", stat="identity") + 
        theme_bw() + 
        facet_wrap(~spliceVariant) + 
        scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
        ylab("Number of events") + 
        xlab(bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-.(value) ~ ")")) + 
        theme_manuscript(fig_font_size=font_size, fig_font=font) + 
        grids(axis="y") + 
        scale_y_log10()
}

g1 <- plotTotalNumberOfEvents(dt2p, "5")
g1 <- plotTotalNumberOfEventsByCategory(dt2p, "5")
g1

g2 <- plotTotalNumberOfEvents(dt2p, "7")
g2 <- plotTotalNumberOfEventsByCategory(dt2p, "7")
g2

g3 <- plotTotalNumberOfEvents(dt2p, "9")
g3 <- plotTotalNumberOfEventsByCategory(dt2p, "9")
g3



#' 
#' Plot percentages
#' 
plotPercentage <- function(dt, value){
    ggplot(dt, aes(y=freq*100, x=names, fill=Method)) + 
        geom_bar(stat="identity", position="dodge") + 
        labs(x=bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-.(value) ~ ")"),
             y="Percentage\nwithin Method") + 
        scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
        theme_manuscript(fig_font_size=font_size, fig_font=font) + 
        grids(axis="y")
}

plotPercentageByCategory <- function(dt, value){
    # combine tissues into categories 1, 2-5, > 5
    dt[as.numeric(names) == 1, category:="1"]
    dt[as.numeric(names) >= 2 & as.numeric(names) <= 5, category:="2-5"]
    dt[as.numeric(names) >= 6 , category:="> 5"]
    dt[, freq_category := sum(freq), by="Method,category"]
    dt_freq <- dt[, .SD[1,.(freq_category=freq_category)], by="Method,category"]
    dt_freq[, category:=factor(category, levels=c("1", "2-5", "> 5"))]
    
    # plot
    ggplot(dt_freq, aes(y=freq_category*100, x=category, fill=Method)) + 
        geom_bar(stat="identity", position="dodge") + 
        labs(x=bquote("Number of tissues outlier is present (" ~ italic(P) < 10^-.(value) ~ ")"),
             y="Percentage\nwithin Method") + 
        scale_fill_manual(values=c("orange", "darkolivegreen", "dodgerblue3", "purple4")) +
        theme_manuscript(fig_font_size=font_size, fig_font=font) + 
        grids(axis="y")
}

gg9p <- plotPercentage(dt2p9p, "9")
gg9p <- plotPercentageByCategory(dt2p9p, "9")
gg9p

gg7p <- plotPercentage(dt2p7p, "7")
gg7p <- plotPercentageByCategory(dt2p7p, "7")
gg7p

gg5p <- plotPercentage(dt2p5p, "5")
gg5p <- plotPercentageByCategory(dt2p5p, "5")
gg5p



#'
#' Arrange the plots
#'
g <- ggarrange(labels=LETTERS[1:4], ncol=1, common.legend=TRUE, legend="bottom",
               font.label=list(size=12, color = "black", face = "bold", family = font),
               g1,
               # g3,
               gg5p,
               gg7p,
               gg9p)
g


#+ save figure
factor <- 0.65
outPng
ggsave(plot=g, filename=outPng, width=page_width, height=1.3*page_width, unit=width_unit, dpi=300)
ggsave(plot=g, filename=outPdf, width=page_width, height=1.3*page_width, unit=width_unit)
