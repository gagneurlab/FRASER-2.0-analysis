### Helper functions for variant enrichment

getEnrichmentPvalMatrix <- function(fds, type, deltaPsiCutoff, minCoverageCutoff, 
                                    rhoCutoff, filterActivation=FALSE,
                                    filterBlacklist=FALSE, 
                                    blacklist_gr=NULL, BPPARAM,
                                    geneColumn="hgnc_symbol", invertFilters=FALSE){
    # retrieve junction level pvalues
    pvals <- as.matrix(pVals(fds, type=type, level="junction"))
    
    # get filter mask (TRUE means that junction does not pass at least one 
    # filter and should be filtered out)
    message(date(), ":   computing filter mask ...")
    filter_mask <- getFilterMask(fds=fds, type=type, 
                                    deltaPsiCutoff=deltaPsiCutoff, 
                                    minCoverageCutoff=minCoverageCutoff, 
                                    rhoCutoff=rhoCutoff, 
                                    filterBlacklist=filterBlacklist, 
                                    blacklist_gr=blacklist_gr,
                                    filterActivation=filterActivation)
    if(isFALSE(invertFilters)){
        pvals[filter_mask] <- NA
    } else{
        pvals[!filter_mask] <- NA
    }
    
    # aggregate to pvalues to splice site level
    message(date(), ":   site-level pval computation ...")
    index <- FRASER:::getSiteIndex(fds, type=type)
    site_pvals <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                          pvals=pvals, index, BPPARAM=BPPARAM,
                          method="holm")
    site_pvals <- do.call(cbind, site_pvals)
    
    # aggregate to pvalues to splice gene level
    message(date(), ":   gene-level pval computation ...")
    gene_pvals <- aggregatePvalToGeneLevel(fds, type=type, BPPARAM=BPPARAM,
                                            site_pvals=site_pvals, index=index,
                                            geneColumn=geneColumn)
    return(gene_pvals)
}

getFilterMask <- function(fds, type, deltaPsiCutoff, minCoverageCutoff, 
                            rhoCutoff, filterBlacklist, blacklist_gr,
                            filterActivation){
    dpsi <- abs(deltaPsiValue(fds, type=type))
    k    <- K(fds, type=type)
    n    <- N(fds, type=type)
    rho  <- matrix(rho(fds, type=type), 
                    nrow=nrow(dpsi), ncol=ncol(dpsi))
    rowMeansK <- matrix(rowMeans(K(fds, type=type)),
                    nrow=nrow(fds), ncol=ncol(fds))
    
    # filter out: 
    #   delta psi < cutoff
    #   coverage (N) < cutoff
    #   rho > cutoff
    mask <- as.matrix(dpsi < deltaPsiCutoff) |
            as.matrix(n < minCoverageCutoff) |
            as.matrix(rho > rhoCutoff) 
    
    # filter out activation outliers if requested
    if(isTRUE(filterActivation)){
        mask <- mask | (as.matrix(k >= 20) & (rowMeansK <= 5))
    }
    
    # get blacklist mask
    if(isTRUE(filterBlacklist)){
        blacklist_mask <- getBlacklistMask(fds, type, blacklist_gr)
        mask <- mask | blacklist_mask
    }
    
    # no filters applied -> cutoff settings:
    #   dPsi      >= 0 (abs(dPsi) always >= 0)
    #   coverage  >= 0 (coverage always >= 0)
    #   rho       <= 1 (rho always <= 1)
    
    return(mask)
}

adjust_FWER_PValues <- function(i, pvals, index, method="holm"){
    dt <- data.table(p=pvals[,i], idx=index, rho=rho)
    suppressWarnings(dt2 <- dt[,.(pa=min(p.adjust(p, method=method), 
                                         na.rm=TRUE)),by=idx])
    dt2[is.infinite(pa), pa:=NA]
    setkey(dt2, "idx")[J(index)][,pa]
}

aggregatePvalToGeneLevel <- function(fds, type, site_pvals, index, geneColumn, 
                                        BPPARAM){
    samples <- samples(fds)
    if(is.null(colnames(site_pvals))){
        colnames(site_pvals) <- samples
    }
    dt <- data.table(
        idx=index,
        geneID=getGeneIDs(fds, type=type, unique=FALSE, 
                          geneColumn=geneColumn),
        as.data.table(site_pvals))
    dt <- dt[!is.na(geneID)]
    geneIDs <- getGeneIDs(fds, type=type, unique=TRUE, 
                          geneColumn=geneColumn)
    
    # separate geneIDs into individual rows
    dt[, dt_idx:=seq_len(.N)]
    dt_tmp <- dt[, splitGenes(geneID), by="dt_idx"]
    dt <- dt[dt_tmp$dt_idx,]
    dt[,`:=`(geneID=dt_tmp$V1, dt_idx=NULL)]
    setkey(dt, geneID)
    
    # aggregate pvalues to gene level per sample
    pvalsPerGene <- genePvalsByGeneID(dt, samples=samples, geneIDs=geneIDs, 
                                      method="holm", BPPARAM=BPPARAM)
    return(pvalsPerGene)
}

getGeneIDs <- function(fds, type, unique=TRUE, geneColumn="hgnc_symbol"){
    geneIDs <- mcols(fds, type=type)[[geneColumn]]
    if(isTRUE(unique)){
        geneIDs <- unique(unlist(lapply(geneIDs, FUN=function(g){
            unlist(strsplit(g, ";"))}) ))
        geneIDs <- geneIDs[!is.na(geneIDs)]
    }
    geneIDs
}

genePvalsByGeneID <- function(dt, samples, geneIDs, method, BPPARAM){
    pvalsPerGene <- bplapply(geneIDs, BPPARAM=BPPARAM,
        FUN=function(g) {
            dttmp <- dt[geneID == g][!duplicated(idx)]
            suppressWarnings(
                pval_g <- apply(as.matrix(dttmp[,-c("idx", "geneID")]), 2,
                    function(x) min(p.adjust(x, method=method), na.rm = TRUE) )
                )
            pval_g[is.infinite(pval_g)] <- NA
            pval_g
    })
    pvalsPerGene <- do.call(rbind, pvalsPerGene)
    rownames(pvalsPerGene) <- geneIDs
    return(pvalsPerGene)
}

splitGenes <- function(x, sep=";"){
    return(unlist(strsplit(as.character(x), sep, fixed=TRUE)))
}

# add the blacklist information
getBlacklistMask <- function(fds, type, blacklist_gr, minoverlap=5){

    
    if(type == "theta"){
        junctions_gr <- rowRanges(fds, type="ss")    
    } else{
        junctions_gr <- rowRanges(fds, type="j")
    } 
    # junctions_gr <- makeGRangesFromDataFrame(junctions_dt)
    
    # get gr with start/end positions of each intron
    gr_start_ss <- junctions_gr
    end(gr_start_ss) <- start(gr_start_ss) + minoverlap - 1
    start(gr_start_ss) <- start(gr_start_ss) - minoverlap
    gr_end_ss <- junctions_gr
    start(gr_end_ss) <- end(gr_end_ss) - minoverlap + 1
    end(gr_end_ss) <- end(gr_end_ss) + minoverlap
    
    # set to the same seqlevelsstyle
    seqlevelsStyle(blacklist_gr) <- seqlevelsStyle(junctions_gr)
    
    ## create overlap with blacklist and annotate extra column
    black_hits_start_ss <- unique(from(findOverlaps(gr_start_ss, blacklist_gr)))
    black_hits_end_ss <- unique(from(findOverlaps(gr_end_ss, blacklist_gr)))
    
    # create mask with junctions lying in blacklist regions
    blacklist_mask <- matrix(FALSE, nrow=nrow(K(fds, type)), ncol=ncol(fds))
    blacklist_mask[black_hits_start_ss,] <- TRUE
    blacklist_mask[black_hits_end_ss,] <- TRUE
    
    return(blacklist_mask)
}
