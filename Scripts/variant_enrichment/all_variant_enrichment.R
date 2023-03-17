#'---
#' title: Get variant enrichment for all datasets
#' author: Ines Scheller
#' wb:
#'   threads: 1
#'   resources:
#'     - mem_mb: 2000
#'   input:
#'     - gg_rds_fixed_q: '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/q{q}/enrichment_plots.Rds", dataset=config["dataset_name"], q=config["q"])`'
#'   output:
#'     - all_done: '`sm expand(config["DATADIR"] + "/{dataset}/var_enrich_all_done.txt", dataset=config["dataset_name"])`'
#'   type: script
#'---

# #'     - gg_rds_all_q:            '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/{dPsi}__{minCoverage}__{rho}/all_q_enrichment_plots.Rds", dataset=config["dataset_name"], dPsi=config["deltaPsi"], minCoverage=config["minCoverage"], rho=config["rho"])`'
# #'     - gg_rds_all_q_per_type:   '`sm expand(config["DATADIR"] + "/{dataset}/plot_rds/{psiType}/{dPsi}__{minCoverage}__{rho}/all_q_enrichment_plots.Rds", dataset=config["dataset_name"], dPsi=config["deltaPsi"], minCoverage=config["minCoverage"], rho=config["rho"], psiType=config["psiTypes"])`'

#+ create dummy file
file.create(snakemake@output$all_done)


