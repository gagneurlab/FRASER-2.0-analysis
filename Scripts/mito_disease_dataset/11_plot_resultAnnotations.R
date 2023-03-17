#'---
#' title: Get result annotation of pathogenic event
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/11_resultAnnotation_{implementation}_minExpr{minK}-quantile{quant}-quantCoverage{minN}.Rds"`'
#'  threads: 1
#'  resources:
#'   - mem_mb: 10000
#'  input:
#'   - resultTableGene_jaccard: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - resultTableJunction_jaccard: '`sm expand(config["mito_processed_results"] + 
#'                          "/results/{implementation}/" + config["mito_annotation"] + "/" +
#'                          config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + 
#'                          "/deltaJaccard{delta}/results_per_junction_blacklist.tsv", delta=config["deltaCutoff"], allow_missing=True)`'
#'   - sample_anno: '/s/project/genetic_diagnosis/rna_paper/sa_solved.tsv'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/outlier_annotation.html"`'
#'   - outlier_annotation: '`sm config["figdir"] + 
#'                "/FRASER_vs_FRASER2/{implementation}/" + 
#'                config["mito_dataset_name"] + "-minExpr{minK}-quantile{quant}-quantCoverage{minN}" + "/" +
#'                config["mito_annotation"] + 
#'                "/outlier_annotation.png"`'
#'  type: noindex
#' output:
#'   html_document
#'---

saveRDS(snakemake, snakemake@log$snakemake)
# snakemake <- readRDS("logs/mito/09_compareVenn_PCA_minExpr20-quantile0.95-quantCoverage1.Rds")

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
# library(VennDiagram)
library(ggplot2)
library(RColorBrewer)

# read in gene-level results tables
res_dt_gene <- fread(snakemake@input$resultTableGene_jaccard)
res_dt_junc <- fread(snakemake@input$resultTableJunction_jaccard)

# read in sample anno of subset and known pathogenic variants
sample_anno <- fread(snakemake@input$sample_anno)
sample_anno[KNOWN_MUTATION == "C19ORF70", KNOWN_MUTATION := "MICOS13"] # symbol was updated
pathogenic_vars <- sample_anno[!is.na(FRASER_padj),]
pathogenic_vars[, sampleGene:=paste(RNA_ID, KNOWN_MUTATION, sep="__")]

# subset and get sample-gene pairs
res_dt_gene <- res_dt_gene[sampleID %in% sample_anno$RNA_ID,]
res_dt_gene[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]
res_dt_gene <- res_dt_gene[!duplicated(res_dt_gene, by=c("sampleID","seqnames","start","end","strand")), ]
res_dt_junc <- res_dt_junc[sampleID %in% sample_anno$RNA_ID,]
res_dt_junc[, sampleGene:=paste(sampleID, hgncSymbol, sep="__")]


# get data for plotting
tdt <- as.data.table(res_dt_junc[, table(spliceType, blacklist)])
tdt[, propTotal := sum(N)/tdt[,sum(N)], by="spliceType"]
tdt[, propBlacklist := .SD[blacklist==TRUE, sum(N)]/tdt[blacklist==TRUE,sum(N)], by="spliceType"]
# tdt[, .(overall=unique(propTotal), blacklist=unique(propBlacklist), N), by="spliceType"]
t2p <- melt(tdt[, .(overall=unique(propTotal), blacklist=unique(propBlacklist)), by="spliceType"], id.vars=c("spliceType"), value.name="proportion")
tdt[blacklist==FALSE, variable:="overall"]
tdt[blacklist==TRUE, variable:="blacklist"]
t2p <- merge(t2p, tdt[, .(spliceType, variable, N)], by=c("spliceType", "variable"))
t2p[, spliceType := factor(spliceType, levels=t2p[variable == "overall",][t2p[variable == "overall", order(-proportion)]][, spliceType])]

# plot
g <- ggplot(t2p, aes(spliceType, proportion, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") + 
    geom_text(aes(label = N, x = spliceType, y = proportion), 
              position = position_dodge(width = 0.8), vjust = -0.6) +
    theme_bw() + 
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          plot.margin = margin(l = 40)) 
g

# pie chart
# Overall piechart
ggplot(t2p[variable == "overall",], aes(x="", y=N, fill=spliceType)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Paired") + 
    ggtitle("Overall proportion in outliers") + 
    theme_void() + # remove background, grid, numeric labels
    theme(plot.title = element_text(hjust = 0.5))

# Blacklist piechart
ggplot(t2p[variable == "blacklist",], aes(x="", y=N, fill=spliceType)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Paired") + 
    ggtitle("Proportion in outliers in blacklist regions") + 
    theme_void() + # remove background, grid, numeric labels
    theme(plot.title = element_text(hjust = 0.5))
    

# get computed splice annotation for pathogenic events
patho_types <- res_dt_gene[sampleGene %in% pathogenic_vars$sampleGene,]
ggplot( as.data.table(patho_types[, table(spliceType)]), 
       aes(x="", y=N, fill=spliceType)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Paired") + 
    ggtitle("Proportion in pathogenic outliers\n(gene level)") + 
    theme_void() + # remove background, grid, numeric labels
    theme(plot.title = element_text(hjust = 0.5))
patho_types_junc <- 
    res_dt_junc[sapply(sampleID, function(sid) any(grepl(sid, pathogenic_vars$sampleGene))) &
            sapply(hgncSymbol, function(gene) any(sapply(strsplit(gene, ";", fixed=T)[[1]], function(g){grepl(g, pathogenic_vars$KNOWN_MUTATION)}))),]
ggplot( as.data.table(patho_types_junc[, table(spliceType)]), 
        aes(x="", y=N, fill=spliceType)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette="Paired") + 
    ggtitle("Proportion in pathogenic outliers\n(junction level)") + 
    theme_void() + # remove background, grid, numeric labels
    theme(plot.title = element_text(hjust = 0.5))

# save plot
out_file <- snakemake@output$outlier_annotation
ggsave(g, filename=out_file)
