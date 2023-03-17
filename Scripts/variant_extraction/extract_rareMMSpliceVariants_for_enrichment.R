#'---
#' title: Extract rare MMSplice variants for enrichment
#' author: Christian Mertes, Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/variant_extraction/extract_rareMMSpliceVars.Rds"`'   
#'   threads: 1
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - variantTable_dir:  "/s/project/absplice/data/results/gtex_v8_with_general_workflow/splicing_predictions/pred_variant_level/mmsplice_baseline"
#'   output:
#'     - variantTable: '`sm expand(config["DATADIR"] + "/{dataset_vcf_group}/variant_extraction/rareMMSplice_filtered_VariantsTable.tsv.gz", dataset_vcf_group=config["vcf_files"].keys())`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
# source("./src/r/config.R")
source("src/R/enrichment_helper.R")
library(data.table)

snptype <- "rareMMSplice"
mmspliceDir <- snakemake@input$variantTable_dir
mmsplice_cutoff <- 2 # as.numeric(snakemake@wildcards$cutoff)
outFile <- snakemake@output$variantTable
# tissues <- snakemake@config$EnrichmentTissues

# input mmsplice variants are already filtered for rare variants only (maf <= 0.001 and max 2 samples with the variant)
MAF_LIMIT <- 0.001

mmspliceDir
outFile
snptype


###########################################
#'
#' # VCF parsing
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- rbindlist(lapply(list.files(mmspliceDir, recursive=TRUE), function(file){
    return(fread(file.path(mmspliceDir, file)))
}))

# apply cutoff on predicted delta_logit_psi
variants <- variants[abs(delta_logit_psi) >= mmsplice_cutoff,]

# cleanup table
variants <- variants[,.(variantID=variant, subjectID=sample, MAF=0.001,
                        delta_logit_psi, Gene=gene_id, SYMBOL=gene_name,
                        simple_conseq="MMSplice", IMPACT="HIGH")]

# hist(variants[,.(MAF=max(MAF)), by=c('subjectID', 'variantID')][,log10(MAF)],
#      breaks=30, main='MAF distribution over all variants of interest',
#      xlab='log10(MAF)')
# abline(v=log10(MAF_LIMIT), col='red')


#' #' MAF filter
#' table(variants[,MAF < MAF_LIMIT])
#' variants <- variants[MAF < MAF_LIMIT]


#' Number of variants per sample
hist(variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N],
     main='Number of variants per sample', xlab='Number of variants')
variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,
   .(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' 
#' * Total number of variants/subject combinations considered in analysis
nrow(variants[,.N,by="variantID,subjectID"])

#'
#' * Delta Psi distribution in logit space
#' 
hist(variants[,delta_logit_psi], breaks=100)

#' Keep only the gene level and most severe variant annotation
#' (remove transcript level)
variants <- variants[order(subjectID, Gene, -abs(delta_logit_psi))]
dupVars <- duplicated(variants[,.(subjectID, Gene)])
table(dupVars)
variantsByGene <- variants[!dupVars]
sort(table(variantsByGene[,.N,by=c('variantID', 'simple_conseq')][,
    simple_conseq], useNA='always'))


#'
#' # Final var overview
#'

#' * Final number of variants per sample
hist(variantsByGene[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N])
variantsByGene[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,
    .(mean=round(mean(N), 2), sd=round(sd(N), 2))]

#' Number of variants sample combinations
nrow(variantsByGene[,.N,by=c('variantID', 'subjectID')])

#' * Number of affected samples:
length(unique(variantsByGene$subjectID))

#' * Number of affected genes:
length(unique(variantsByGene$SYMBOL))
length(unique(variantsByGene$Gene))

#'
#' # Save variant table
#' 
fwrite(variantsByGene[,.(variantID, subjectID, MAF, Gene, 
    SYMBOL, simple_conseq, IMPACT, Consequence="MMSplice")], outFile)
