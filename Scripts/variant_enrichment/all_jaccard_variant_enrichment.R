#'---
#' title: Get all jaccard variant enrichment results
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 2000
#'   input:
#'     - gg_rds_fixed_q: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/q{q}/jaccard_enrichment_plots.Rds", dataset=config["tissues_for_reproducibility"], q=config["q"])`'
#'     - gg_rds_all_q:            '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/{dPsi}__{minCoverage}__{rho}/FRASER2_all_q_enrichment_plots.Rds", dataset=config["tissues_for_reproducibility"], dPsi=config["deltaPsi"], minCoverage=config["minCoverage"], rho=config["rho"])`'
#'   output:
#'     - all_done: '`sm expand(config["DATADIR"] + "/{dataset}/jaccard_var_enrich_all_done.txt", dataset=config["tissues_for_reproducibility"])`'
#'   type: script
#'---

#+ create dummy file
file.create(snakemake@output$all_done)


