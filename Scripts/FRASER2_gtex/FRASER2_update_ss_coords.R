#'---
#' title: Update splice site coordinates to add strand info for gene annotation
#' author: Ines Scheller
#' wb:
#'   threads: 5
#'   resources:
#'     - mem_mb: 50000
#'   input:
#'     - fds_raw: '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/jaccard.h5"`'
#'   output:
#'     - ss_update_done: '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/ss_update.done"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)

#+ load raw fds (with jaccard)
fds <- loadFraserDataSet(file=snakemake@input$fds_raw)
ss_counts <- K(fds, "theta")
gr_prev <- rowRanges(fds, type="j")
gr_ss_prev <- rowRanges(fds, type="ss")

#+ create updated ss coordinates (with strand info)
# get intron genomic ranges and annotate with startID/endID
intronRanges_all <- FRASER:::annotateSpliceSite(granges(gr_prev))
# get ss ranges for introns passing K filter
minExpressionInOneSample <- 20
splitCounts <- K(fds, type="psi3")
maxCount <- rowMaxs(splitCounts)
passed <- maxCount >= minExpressionInOneSample
intronRanges <- intronRanges_all[passed,]
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(intronRanges)

#+ update non split counts to match new ss coords
map <- findOverlaps(spliceSiteCoords, gr_ss_prev, type="equal")
ss_counts_new <- ss_counts[to(map),]

#+ update startID and endID to new ones
startIDs_new <- intronRanges_all$startID
endIDs_new <- intronRanges_all$endID
mcols(fds, type="j")$startID <- startIDs_new
mcols(fds, type="j")$endID <- endIDs_new

#+ update fds with new ss coords and ss counts
# need to create new nonSplicedReads summarizedExperiment object due to dimension change
aName <- "rawCountsSS"
h5 <- FRASER:::saveAsHDF5(fds, aName, as.matrix(ss_counts_new))
colnames(h5) <- samples(fds)
aList <- list(a=h5)
names(aList) <- aName
ss_se <- SummarizedExperiment(
    colData=colData(fds), assays=aList, rowRanges=spliceSiteCoords)
nonSplicedReads(fds) <- ss_se

#+ save updated fds object
fds <- saveFraserDataSet(fds)
file.create(snakemake@output$ss_update_done)
