#'---
#' title: FDR only on genes with rare variant im OMIM gene
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '`sm config["log_dir"] + "/mito/00_create_DNA_IDs_list.Rds"`'
#'  threads: 1
#'  resources:
#'   - mem_mb: 4000
#'  input:
#'   - sample_anno: '`sm config["mito_sample_anno"]`'
#'   - full_sample_anno: '`sm config["mito_full_sample_anno"]`'
#'  output:
#'   - dna_id_list: '`sm config["mito_processed_data"] + "/variants/prokisch_dna_ids.txt"`'
#'  type: script
#'---


#+ load packages
library(data.table)

#+ load sample anno
sa <- fread(snakemake@input$sample_anno)
full_sa <- fread(snakemake@input$full_sample_anno)

#+ extract needed dna ids
dna_ids <- full_sa[RNA_ID %in% sa$RNA_ID, .(DNA_ID)]

#+ save dna id list
fwrite(dna_ids, file=snakemake@output$dna_id_list, col.names=F)
