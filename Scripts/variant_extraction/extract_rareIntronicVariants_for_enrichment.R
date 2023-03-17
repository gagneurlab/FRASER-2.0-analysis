#'---
#' title: Extract rare intronic variants for enrichment
#' author: Christian Mertes
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/variant_extraction/extract_rareIntronVars.Rds"`'   
#'   threads: 1
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - variantTable: '`sm expand(config["DATADIR"] + "/{dataset_vcf_group}/variant_extraction/rareIntronicVariantsTable.tsv", dataset_vcf_group=config["vcf_files"].keys())`' 
#'   output:
#'     - variantTable: '`sm expand(config["DATADIR"] + "/{dataset_vcf_group}/variant_extraction/rareIntronic_filtered_VariantsTable.tsv.gz", dataset_vcf_group=config["vcf_files"].keys())`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

saveRDS(snakemake, snakemake@log$snakemake)

# source config
library(data.table)
source("src/R/enrichment_helper.R")

snptype <- "rareIntronic"
vcfFile <- snakemake@input$variantTable
outFile <- snakemake@output$variantTable

MAF_LIMIT <- as.numeric(snakemake@config$MAF_LIMIT)

vcfFile
outFile
snptype

###########################################
#'
#' # VCF parsing
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(vcfFile)

#' Remove NA variants
table(variants$GT)
variants <- variants[GT != "./."]

#' Get only splice region variants
sort(table(gsub("&.*", "", variants$Consequence)))
variants <- variants[grepl("splice_", Consequence)]

#' Process the MAF column
#' * MAF for standard annotated variants
#' * 1/(2*#samples) for not annotated variants
#' * max(freq(allels)) for multi allele locations
# variants[,MAF:=pmax(AF, maf2number(MAX_AF), na.rm=TRUE)]
variants[,MAF:=pmax(AF, MAX_AF, na.rm=TRUE)]

hist(variants[,.(MAF=max(MAF)), by=c('subjectID', 'variantID')][,log10(MAF)],
     breaks=30, main='MAF distribution over all variants of interest',
     xlab='log10(MAF)')
abline(v=log10(MAF_LIMIT), col='red')

#' MAF filter
table(variants[,MAF < MAF_LIMIT])
variants <- variants[MAF < MAF_LIMIT]


#' Number of variants per sample
hist(variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N],
     main='Number of variants per sample', xlab='Number of variants')
variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,
                                                               .(mean=round(mean(N), 2), sd=round(sd(N), 2))]


#' 
#' * Total number of variants/subject combinations considered in analysis
nrow(variants[,.N,by="variantID,subjectID"])

#'
#' * Total number of gene/subject combinations
nrow(variants[,.N,by="SYMBOL,subjectID"])

#'
#' * Simplify Consequence annotation
#+ simplify consequences
variants <- simplifyConsequences(variants, splicingFirst=FALSE)

#' Keep only the gene level and most severe variant annotation
#' (remove transcript level)
variants <- variants[order(subjectID, Gene, simple_conseq, rank)]
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
                         SYMBOL, simple_conseq, IMPACT, Consequence)], outFile)



