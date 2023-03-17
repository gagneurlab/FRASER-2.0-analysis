#'---
#' title: Filter junctions
#' author: Ines Scheller
#' wb:
#'   py:
#'   - |
#'    def get_input_fds(wildcards):
#'     return config["datasets"][wildcards.dataset]["raw_fds_file"]
#'   threads: 10
#'   resources:
#'     - mem_mb: 125000
#'   input:
#'     - fds_file: '`sm get_input_fds`'
#'   output:
#'     - fds_out: '`sm config["DATADIR"] + "/{dataset_group}/fds/savedObjects/raw-{dataset}__jaccard/jaccard.h5"`'
#'   type: script
#'---

#+ load FRASER
.libPaths("~/R/4.1/FRASER2")
library(FRASER)

#+ input
fds_file <- snakemake@input$fds_file
nthreads <- snakemake@threads
register(MulticoreParam(nthreads))

#+ output
out_file <- snakemake@output$fds_out[1]
out_dir <- dirname(dirname(dirname(out_file)))
out_name <- basename(dirname(out_file))

#+ load fds
fds <- loadFraserDataSet(file=fds_file)
message("Loaded fds.")

#+ set names and dir to new ones
message("Setting new  name of fds to ", out_name)
name(fds) <- as.character(out_name)
message("Setting workingDir of fds to ", out_dir)
workingDir(fds) <- out_dir

#+ remove assays that are not needed or will be recomputed
anames <- assayNames(fds)
assays_to_remove <- c(anames[grepl("padj", anames)], anames[grepl("pvalues", anames)], anames[grepl("zScores", anames)], anames[grepl("predictedMeans", anames)])
for(aname in assays_to_remove){
    assay(fds, aname, withDimnames = FALSE) <- NULL
}
message("Removed unneccesarry assays.")

##### quick fix for missing startIDs and endIDs:
spliceSiteCoords <- rowRanges(fds, type="ss")
# # check how many spliceSiteIDs are mismapped
# table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
#            "startID_found"=mcols(fds, type="j")$startID %in% spliceSiteCoords$spliceSiteID))
# table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
#            "endID_found"=mcols(fds, type="j")$endID %in% spliceSiteCoords$spliceSiteID))

# correct mapping
correct_mapping <- lapply(c("startID", "endID"), function(id_type){
    id_missing <- which(!mcols(fds, type="j")[, id_type] %in% spliceSiteCoords$spliceSiteID)
    missing_gr <- rowRanges(fds, type="j")[id_missing,]
    if(id_type == "startID"){
        end(missing_gr) <- start(missing_gr)
        start(missing_gr) <- start(missing_gr)-1
    } else{
        start(missing_gr) <- end(missing_gr)
        end(missing_gr) <- end(missing_gr) + 1        
    }
    missing_found <- findOverlaps(missing_gr, spliceSiteCoords, type="equal")
    return(list(mapping=missing_found, missing_ids=id_missing))
})
names(correct_mapping) <- c("startID", "endID")
mcols(fds, type="j")$startID[correct_mapping$startID$missing_ids[from(correct_mapping$startID$mapping)]] <- 
    spliceSiteCoords$spliceSiteID[to(correct_mapping$startID$mapping)]
mcols(fds, type="j")$endID[correct_mapping$endID$missing_ids[from(correct_mapping$endID$mapping)]] <- 
    spliceSiteCoords$spliceSiteID[to(correct_mapping$endID$mapping)]

# check if successful
message(table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
    "startID_found"=mcols(fds, type="j")$startID %in% spliceSiteCoords$spliceSiteID)))
message(table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
    "endID_found"=mcols(fds, type="j")$endID %in% spliceSiteCoords$spliceSiteID)))

# just jaccard index calculation (psi3/5 and theta already done)
# calculate intron jaccard index
fds <- FRASER:::calculateIntronNonsplitSum(fds, overwriteCts=TRUE)
fds <- FRASER:::calculateJaccardIntronIndex(fds, overwriteCts=TRUE)
fds <- FRASER:::calculateDeltaPsiValue(fds, "jaccard", "delta_jaccard")

#+ write output
fds <- saveFraserDataSet(fds)
