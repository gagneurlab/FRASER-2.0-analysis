#'---
#' title: Compare with FRASER2 with old FRASER results on UDN dataset
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/snakemake/udn/compare_nrOutliers_F2_delta{delta}_{implementation}_minK{k}_{quant}_minN{n}.Rds"`'
#'  threads: 4
#'  resources:
#'   - mem_mb: 16000
#'  input:
#'   - res_tables_fraser2: '`sm expand(config["DATADIR"] + 
#'                          "/udn/FRASER2_results/minK{k}_{quant}_minN{n}/{implementation}/" +
#'                          "{dataset}/optQ__newFilt/delta{delta}/results_gene.tsv", 
#'                          dataset=config["udn_datasets"], allow_missing=True)`'
#'   - res_tables_fraser2_omim: '`sm expand(config["DATADIR"] + 
#'                          "/udn/FRASER2_results/minK{k}_{quant}_minN{n}/{implementation}/" +
#'                          "{dataset}/optQ__newFilt/delta{delta}/results_gene_FDRomim.tsv", 
#'                          dataset=config["udn_datasets"], allow_missing=True)`'
#'   - res_tables_fraser1: '`sm expand(config["general_data_dir"] + "/genetic_diagnosis/udn/processed_results/aberrant_splicing/results/v29/fraser/{dataset}/results.tsv", 
#'                          dataset=config["udn_datasets"])`'
#'   - fds_fraser2: '`sm expand(config["DATADIR"] + 
#'                          "/udn/fds/minK{k}_{quant}_minN{n}/{implementation}/savedObjects/" +
#'                          "{dataset}__optQ__newFilt/fds-object.RDS", 
#'                          dataset=config["udn_datasets"], allow_missing=True)`'
#'  output:
#'   - ggplots: '`sm config["DATADIR"] + "/udn/fraser2_improvements/minK{k}_{quant}_minN{n}/{implementation}/delta{delta}/nrOutliers_comparison_ggplot.Rds"`'
#'  type: script
#'---


saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(FRASER)
library(data.table)
library(ggplot2)
library(ggpubr)

# read in gene-level results tables
res_tables_fraser2 <- rbindlist(lapply(snakemake@input$res_tables_fraser2, FUN=function(filename){
    dt <- fread(filename)
    # remove duplicate entries for same outlier when overlapping several genes
    dt <- dt[!duplicated(dt, by=c("sampleID","seqnames","start","end","strand")), ]
    # get nr of outliers per sample (gene-level)
    dt[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
    dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    # add name of dataset
    ds_name <- basename(dirname(dirname(dirname(filename))))
    dt[, dataset_name:=ds_name]
    dt[, method:="FRASER2"]
    return(dt[,unique(numOutliersPerSample), by="sampleID,method,dataset_name"])
    }))
res_tables_fraser2_omim <- rbindlist(lapply(snakemake@input$res_tables_fraser2_omim, FUN=function(filename){
    dt <- fread(filename)
    # remove duplicate entries for same outlier when overlapping several genes
    dt <- dt[!duplicated(dt, by=c("sampleID","seqnames","start","end","strand")), ]
    # get nr of outliers per sample (gene-level)
    dt[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
    dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    # add name of dataset
    ds_name <- basename(dirname(dirname(dirname(filename))))
    dt[, dataset_name:=ds_name]
    dt[, method:="FRASER2 (OMIM)"]
    return(dt[,unique(numOutliersPerSample), by="sampleID,method,dataset_name"])
}))
res_tables_fraser1 <- rbindlist(lapply(snakemake@input$res_tables_fraser1, FUN=function(filename){
    dt <- fread(filename)
    # get nr of outliers per sample (gene-level)
    dt[, numOutliersPerSample := uniqueN(hgncSymbol), by = "sampleID"]
    dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    # add name of dataset
    ds_name <- basename(dirname(filename))
    dt[, dataset_name:=ds_name]
    dt[, method:="FRASER"]
    return(dt[,unique(numOutliersPerSample), by="sampleID,method,dataset_name"])
}))

#+ read fraser2 fds objects to get total nr of samples
fds_ls <- lapply(snakemake@input$fds_fraser2, function(filename){
    fds <- loadFraserDataSet(file=filename)
    return(fds)
})
names(fds_ls) <- sapply(snakemake@input$fds_fraser2, function(filename){
    fds_name <- basename(dirname(filename))
    ds_name <- strsplit(fds_name, "__", fixed=TRUE)[[1]][1]
    return(ds_name)
})
n_samples <- lapply(fds_ls, ncol)
all_samples_dt <- rbindlist(lapply(names(fds_ls), function(ds_name){
    fds <- fds_ls[[ds_name]]
    return(data.table(sampleID=samples(fds), dataset_name=ds_name))
}))

#+ create plot
plot_dt <- dcast(rbind(res_tables_fraser2, res_tables_fraser2_omim, res_tables_fraser1), sampleID + dataset_name ~ method, value.var="V1")
plot_dt <- merge(plot_dt, 
                 all_samples_dt,
                 by=c("sampleID", "dataset_name"), all.y=TRUE)
plot_dt[is.na(FRASER2), FRASER2 := 0]
plot_dt[is.na(`FRASER2 (OMIM)`), `FRASER2 (OMIM)` := 0]
plot_dt[is.na(FRASER), FRASER := 0]
plot_dt[, dataset_name := sapply(dataset_name, FUN=function(lab){
                                    switch(lab, 
                                        blood_polya = "Blood poly(A)",
                                        blood_total_rna = "Blood totalRNA",
                                        lab)
                                })]
plot_dt[, dataset_label := paste0(dataset_name, " (N=", uniqueN(sampleID), ")"), by="dataset_name"]
g_numOut_gene <- ggplot(plot_dt, aes(FRASER+1, FRASER2+1)) +
    facet_wrap(~dataset_label) +
    geom_hex(bins = 50) +
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    geom_hline(data=plot_dt[,median(FRASER2), by="dataset_label"], aes(yintercept=V1), linetype="dashed", col="firebrick") +
    geom_vline(data=plot_dt[,median(FRASER), by="dataset_label"], aes(xintercept=V1), linetype="dashed", col="firebrick") +
    labs(x="FRASER splicing outliers per sample + 1", y="FRASER2 splicing outliers per sample + 1") +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks(sides="bl") +
    scale_fill_continuous(type = "viridis") +
    #scale_fill_distiller(palette= "Spectral") +
    theme_pubr() +
    theme(text=element_text(size=14), 
          axis.title=element_text(size=14, face="bold"))


ggplots <- list()
ggplots[["g_numOut_gene"]] <- g_numOut_gene

#+ boxplots of outliers per sample
boxplot_dt <- copy(plot_dt)
boxplot_dt <- melt(boxplot_dt[, .(sampleID, dataset_name, dataset_label, FRASER, FRASER2, `FRASER2 (OMIM)`)], id.vars=c("sampleID", "dataset_name", "dataset_label"), variable.name="method", value.name="num_out")
boxplot_dt[, method:=factor(method, levels=c("FRASER", "FRASER2", "FRASER2 (OMIM)"))]
boxplot_dt[, dataset_label:=factor(dataset_label, levels=sort(unique(dataset_label)))]
g_boxplot <- ggplot(boxplot_dt, aes(x=method, y=num_out + 1, col=method) ) +
    facet_wrap(~dataset_label) + 
    geom_violin() +
    geom_boxplot(width=0.25, alpha=0.2) + # , color="black"
    labs(x="method",
         y="Splicing outliers per sample + 1") +
    scale_y_log10() +
    # annotation_logticks(sides="l") +
    labs(x="") +
    # scale_fill_brewer(palette="Paired") + 
    scale_color_manual(values=c("dodgerblue3", "purple4", "purple1")) + 
    theme_pubr() + 
    theme(legend.position="none")
a <- annotation_logticks(sides='l')
a$data <- data.frame(x=NA, dataset_label=c(sort(boxplot_dt[, unique(dataset_label)])[1]))
g_boxplot <- g_boxplot + a
g_boxplot

ggplots[["boxplots_numOut_gene"]] <- g_boxplot

#+ save plot as rds
saveRDS(ggplots, file=snakemake@output$ggplots)
