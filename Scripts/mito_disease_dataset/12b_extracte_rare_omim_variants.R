#'---
#' title: Produce table of all rare variants at certain MAF threshold
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/00_get_rare_variants_MAF{maf}.Rds"`'
#'  threads: 5
#'  resources:
#'   - mem_mb: 20000
#'  py:
#'  - | 
#'   with open(config["mito_processed_data"] + '/variants/prokisch_dna_ids.txt') as f:
#'    dna_ids = f.read().splitlines()
#'  input:
#'   - dna_id_list: '`sm config["mito_processed_data"] + "/variants/prokisch_dna_ids.txt"`'
#'   - sample_variants: '`sm expand("/s/project/mitoMultiOmics/raw_data/helmholtz/{dnaID}/exomicout/paired-endout/processedData/vep_anno_{dnaID}_uniq_dt.Rds", dnaID=dna_ids)`'
#'   - omim_genes: '`sm config["omim_genes"]`'
#'  output:
#'   - rare_omim_variants: '`sm config["mito_processed_data"] + "/variants/rare_omim_variants_MAF{maf}_all.Rds"`'
#'   - rare_omim_variants_no_intergenic: '`sm config["mito_processed_data"] + "/variants/rare_omim_variants_MAF{maf}_no_intergenic.Rds"`'
#'   - rare_omim_variants_no_utr: '`sm config["mito_processed_data"] + "/variants/rare_omim_variants_MAF{maf}_no_utr.Rds"`'
#'  type: script
#'---


#+ load packages
library(data.table)
library(BiocParallel)

#+ extract desired MAF threshold
MAF <- as.numeric(snakemake@wildcards$maf)

#+ read list of omim files
omim_genes <- fread(snakemake@input$omim_genes)

#+ subset to rare variants only
all_rare_omim_vars <- rbindlist (bplapply(snakemake@input$sample_variants, function(var_file){
    variants <- readRDS(var_file)
    variants[hgncid == "C19orf70", hgncid := "MICOS13"]
    variants_rare <- variants[gnomAD_MAX_AF < MAF | is.na(gnomAD_MAX_AF),]
    # variants_rare
    variants_rare_omim <- variants_rare[hgncid %in% omim_genes$gene_v29 | hgncid %in% c("MICOS13", "PTCD3", "MRPS25"),]
    return(variants_rare_omim)
}))

#+ cohort allele frequency filtering
cohort_size <- all_rare_omim_vars[, uniqueN(sample)]
all_rare_omim_vars[, nr_samples_per_var := .SD[,uniqueN(sample)], by="chr,pos,ref,alt"]
all_rare_omim_vars[, cohort_AF := nr_samples_per_var / cohort_size]
all_rare_omim_vars <- all_rare_omim_vars[nr_samples_per_var <= 3,]
# all_rare_omim_vars <- all_rare_omim_vars[cohort_AF <= 0.01,]

#+ remove intergenic and utr variants and save as different variant sets
rare_omim_vars_no_intergenic <- all_rare_omim_vars[!mstype %in% c("downstream", "upstream", "intergenic"),]
rare_omim_vars_no_utr <- rare_omim_vars_no_intergenic[!grepl("utr", mstype),]

#+ save output
saveRDS(all_rare_omim_vars, file=snakemake@output$rare_omim_variants)
saveRDS(rare_omim_vars_no_intergenic, file=snakemake@output$rare_omim_variants_no_intergenic)
saveRDS(rare_omim_vars_no_utr, file=snakemake@output$rare_omim_variants_no_utr)
