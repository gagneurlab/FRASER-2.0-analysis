

calculateEnrichment <- function(dt, cols, cutOff, isZscore=TRUE, na.rm=FALSE, multiCall=0){
    dttmp <- copy(dt[,.(subjectID, geneID, tissue, simple_conseq, score=get(cols))])
    if(isTRUE(isZscore)){
        if(isTRUE(na.rm)){
            dttmp[is.na(score), c(cols):=0]
        }
        dttmp[,cutoff:=abs(score) > cutOff]
    } else {
        if(isTRUE(na.rm)){
            dttmp[is.na(score), c(cols):=1]
        }
        dttmp[,cutoff:=score < cutOff]
    }
    if(multiCall > 0){
        cols_multi <- paste0(cols, "_numMultiCall")
        if(!cols_multi %in% colnames(dt)){
            warning(cols_multi, " not there in table!")
            return(NULL)
        }
        dttmp[,mtoc:=dt[,get(cols_multi) >= multiCall]]
        dttmp[cutoff == TRUE & is.na(mtoc),   cutoff:=FALSE]
        dttmp[cutoff == TRUE & mtoc == FALSE, cutoff:=FALSE]
        dttmp[cutoff == TRUE]
    }
    
    ans <- dttmp[, .(nRareEvent=sum(!is.na(simple_conseq)), 
                     total=.N, nNA=sum(is.na(score))), by=cutoff]
    if(!any(ans$cutoff == FALSE, na.rm = TRUE)){
        ans <- rbind(ans, data.table(cutoff=FALSE, nRareEvent=0, total=0, nNA=0))
    }
    if(!any(ans$cutoff == TRUE, na.rm = TRUE)){
        ans <- rbind(ans, data.table(cutoff=TRUE, nRareEvent=0, total=0, nNA=0))
    }
    ans <- ans[,.(nRareEvent, total, fraction=(nRareEvent)/(total), nNA), by=cutoff]
    enrich <- ans[cutoff==TRUE,fraction]/ans[cutoff==FALSE,fraction]
    
    # get bound as done in Li et al
    out.var      <- ans[cutoff == TRUE,  nRareEvent]
    out.total    <- ans[cutoff == TRUE,  total]
    nonout.var   <- ans[cutoff == FALSE, nRareEvent]
    nonout.total <- ans[cutoff == FALSE, total]
    log.se = sqrt(1/out.var - 1/out.total + 1/nonout.var - 1/nonout.total)
    max.ci = enrich * exp(1.96*log.se)
    min.ci = enrich * exp(-1.96*log.se)
    
    return(list(dt=ans, enrichment=enrich, max.ci=max.ci, min.ci=min.ci, 
                multiCall=multiCall))
}

getEnrichmentForMethod <- function(dt, method, cutoffs=c(0.005, 2)){
    message(method)
    isZscore <- grepl('_z$', method)
    curCutoff <- cutoffs[isZscore + 1]
    curType <- paste0(c('P-value (<', 'Z score (>'), cutoffs, ')')[isZscore + 1]
    ans <- calculateEnrichment(dt, method, curCutoff, isZscore)
    return(data.table(
        Method=gsub('_[pz]$', '', method), Type=curType, Cutoff=curCutoff,
        enrichment=ans$enrichment, max.ci=ans$max.ci, min.ci=ans$min.ci,
        out.var      = ans$dt[cutoff == TRUE,  nRareEvent],
        out.total    = ans$dt[cutoff == TRUE,  total],
        nonout.var   = ans$dt[cutoff == FALSE, nRareEvent],
        nonout.total = ans$dt[cutoff == FALSE, total]))
}

getEnrichmentForTissues <- function(tissue, rdsFiles, cutoffs){
    message(tissue)
    rds <- readRDS(rdsFiles[[tissue]])
    methods <- grep('_[pz]$', colnames(rds), value=TRUE)
    enrichments <- rbindlist(lapply(methods, getEnrichmentForMethod,
                                    dt=rds, cutoffs=cutoffs))
    enrichments[,Tissue:=tissue]
    return(enrichments)
}

plotEnrichment <- function(dt, numVarOffset=0.8){
    
    ggplot(dt, aes(x=Method, y=enrichment, col=Type)) +
        geom_point(position=position_dodge(.75)) +
        geom_errorbar(aes(ymin=min.ci, ymax=max.ci, col=Type), width=.2,
                      position=position_dodge(.75)) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        geom_hline(yintercept=1, col='orange') +
        geom_point(aes(Method, numVarOffset, col=Type, shape=numVars),
                   position=position_dodge(.75)) +
        grids() +
        facet_wrap('Tissue')
}

calculateRecallRank <- function(dt, cols, isPvalue=TRUE){
    dt <- dt[order(abs(get(cols)), decreasing=!isPvalue)]
    dt[,c(paste0(cols, '_recall')):=cumsum(!is.na(simple_conseq))]
    dt[,c(paste0(cols, '_rank')):=seq_len(.N)]
    
    return(dt)
}

simplifyConsequences <- function(vardt, splicingFirst=FALSE){
    SORanking <- c(
        # HIGH
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        # MODERAT
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        # LOW
        "splice_region_variant",
        "synonymous_variant")
    if(isTRUE(splicingFirst)){
        idx <- startsWith(SORanking, "splice_")
        SORanking <- c(SORanking[idx], SORanking[!idx])
    }
    
    vardt[,simple_conseq:=NA_character_]
    sapply(SORanking, function(x){
        vardt[is.na(simple_conseq) & grepl(x, Consequence),
              c("rank", "simple_conseq"):=list(which(SORanking==x), x)]
    })
    
    vardt[,simple_conseq:=factor(simple_conseq, 
                                 levels=unique(c(SORanking, simple_conseq)))]
    return(vardt)
}

plotRecallRankForEnrichment <- function(dt, maxRank, maxPoints, logy=FALSE, logx=FALSE){
    if(!missing(maxRank)){
        dt <- dt[rank<maxRank]
    }
    if(!missing(maxPoints)){
        prob4Samp <- min(1, maxPoints / (nrow(dt)/nrow(unique(dt[,.(Method, Type)]))))
        dt <- dt[
            rank < 1e3 |
                rank > max(rank) - 1e3 |
                sample(c(TRUE, FALSE), .N, replace = TRUE, prob = c(prob4Samp, 1 - prob4Samp))]
    }
    
    gg <- ggplot(dt[,.(rank=rank, recall=recall, Method, Type)],
                 aes(rank, recall, col=Method, linetype=Type)) +
        geom_line() + grids()
    
    if(isTRUE(logx)){
        gg <- gg + scale_x_log10()
    }
    if(isTRUE(logy)){
        gg <- gg + scale_y_log10()
    }
    gg
}

maf2number <- function(maf, aggregateFun=max, naMAF=FALSE){
    mafls <- strsplit(maf, '&')
    ans <- unlist(bplapply(mafls, agrFun=aggregateFun,
                           function(x, agrFun){ agrFun(as.numeric(gsub('[A-Z-]+:', '', x)))}))
    if(isScalarNumeric(naMAF)){
        ans[is.na(ans)] <- naMAF
    }
    return(ans)
}

readVcfParallel <- function(vcfFile, ranges, threads, nPerChunk=25,
                            vcfParam=ScanVcfParam()){
    rangesOfInt <- reduce(ranges)
    nchunks <- max(10, ceiling(length(rangesOfInt)/nPerChunk))
    BPPARAM <- MulticoreParam(threads, nchunks, progressbar=TRUE)
    chunks <- chunk(seq_along(rangesOfInt), n.chunks=nchunks)
    
    message(date(), ':',
            ' Run with number of chunks: ', length(chunks),
            ' and size of chunks: ', length(chunks[[1]]),
            ' with threads in parallel: ', threads)
    
    if(isFALSE(bpisup(BPPARAM))){
        bpstart(BPPARAM)
    }
    # read vcf file in chunks
    vcfgtls <- bplapply(chunks, f=vcfFile, ranges=rangesOfInt, param=vcfParam,
                        BPPARAM=BPPARAM,
                        function(x, ranges, f, param){
                            vcfWhich(vcfParam) <- ranges[x]
                            readVcf(TabixFile(f), param=vcfParam) })
    message(date(), ': done extracting genotypes')
    vcfgt <- do.call(rbind, vcfgtls)
    
    bpstop(BPPARAM)
    
    # remove unwanted variants lying in the same region
    vcfgt <- vcfgt[unique(names(ranges))]
    return(vcfgt)
}
