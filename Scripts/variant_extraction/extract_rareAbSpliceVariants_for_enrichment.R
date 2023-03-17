#'---
#' title: Extract rare AbSplice variants for enrichment
#' author: Ines Scheller
#' wb:
#'   log:
#'    - snakemake: '`sm config["log_dir"] + "/variant_extraction/extract_rareAbSpliceVars.Rds"`'   
#'   threads: 1
#'   resources:
#'     - mem_mb: 15000
#'   input:
#'     - variantTable:  "/s/project/absplice/data/results/gtex_v8_with_general_workflow/benchmark_revision/benchmark_fraser_1_outliers.csv"
#'     - gencode: "/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode34/gene_name_mapping_gencode34.tsv"
#'   output:
#'     - variantTable: '`sm expand(config["DATADIR"] + "/{dataset_vcf_group}/variant_extraction/rareAbSplice_filtered_VariantsTable.tsv.gz", dataset_vcf_group=config["vcf_files"].keys(), allow_missing=True)`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# source config
# source("./src/r/config.R")
source("src/R/enrichment_helper.R")
library(data.table)

snptype <- "rareAbSplice"
abspliceFile <- snakemake@input$variantTable
absplice_cutoff <- 0.05 # high: 0.2, medium: 0.05, low: 0.01
outFile <- snakemake@output$variantTable
# tissues <- snakemake@config$EnrichmentTissues

# input mmsplice variants are already filtered for rare variants only (maf <= 0.001 and max 2 samples with the variant)
MAF_LIMIT <- 0.001

###########################################
#'
#' # VCF parsing
#'
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(abspliceFile)

# ignore testis tissue
variants <- variants[tissue != "Testis",]

# apply cutoff on predicted absplice score (max score across tissues)
variants[, AbSplice_DNA_max := max(AbSplice_DNA_absplice), by="gene_id,sample"]
variants[!duplicated(variants, by=c("gene_id", "sample"))]
variants <- variants[AbSplice_DNA_max >= absplice_cutoff,]

#
gencode <- fread(snakemake@input$gencode)
#gencode <- fread("/s/project/gtex_genetic_diagnosis/v8/processed_data/aberrant_expression/gencode34/gene_name_mapping_gencode34.tsv")
gencode <- gencode[, .(gene_id, gene_name)]
gencode[, gene_id := gsub("\\..*$", "", gene_id)]
variants <- merge(variants, gencode, by="gene_id")

# cleanup table
variants <- variants[,.(variantID=NA, subjectID=sample, MAF=0.001,
                        AbSplice_DNA_max, Gene=gene_id, SYMBOL=gene_name,
                        simple_conseq="AbSplice", IMPACT="HIGH")]

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
hist(variants[,AbSplice_DNA_max], breaks=100)

#' Keep only the gene level and most severe variant annotation
#' (remove transcript level)
variants <- variants[order(subjectID, Gene, -AbSplice_DNA_max)]
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
                         SYMBOL, simple_conseq, IMPACT, Consequence="AbSplice")], outFile)
