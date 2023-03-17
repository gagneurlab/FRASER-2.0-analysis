#'---
#' title: Calculate Jaccard Intron Index values
#' author: Ines Scheller
#' wb:
#'  log:
#'   - snakemake: '`sm config["log_dir"] + "/mito/02_Jaccardcalc.Rds"`'
#'  params:
#'   - workingDir: '`sm config["mito_processed_data"] + "/datasets/"`'
#'  resources:
#'   - mem_mb: 125000
#'  threads: 10
#'  input:
#'   - theta:     '`sm config["mito_processed_data"] +
#'                    "/datasets/savedObjects/raw-" + 
#'                    config["mito_dataset_name"] + "/theta.h5"`'
#'  output:
#'  - intron_jaccard: '`sm config["mito_processed_data"] + 
#'                    "/datasets/savedObjects/raw-" + 
#'                    config["mito_dataset_name"] + "/jaccard.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, snakemake@log$snakemake)

# load FRASER version with jaccard index
.libPaths("~/R/4.1/FRASER2")
library(DelayedArray)
library(FRASER)

dataset    <- snakemake@config$dataset_name
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-", dataset))

#+ quick fix for missing startIDs and endIDs:
spliceSiteCoords <- rowRanges(fds, type="ss")
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

# #+ check if successful
# table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
#            "startID_found"=mcols(fds, type="j")$startID %in% spliceSiteCoords$spliceSiteID))
# table(list("passedExpression"=mcols(fds, type="j")$passedExpression, 
#            "endID_found"=mcols(fds, type="j")$endID %in% spliceSiteCoords$spliceSiteID))


# Calculating jaccard intron index values (part of psi value calc function)
# fds <- calculatePSIValues(fds)

# just jaccard index calculation (psi3/5 and theta already done)
# calculate intron jaccard index
fds <- FRASER:::calculateIntronNonsplitSum(fds, overwriteCts=TRUE)
fds <- FRASER:::calculateJaccardIntronIndex(fds, overwriteCts=TRUE)
fds <- FRASER:::calculateDeltaPsiValue(fds, "jaccard", "delta_jaccard")

# FRASER object after PSI value calculation
fds <- saveFraserDataSet(fds)
