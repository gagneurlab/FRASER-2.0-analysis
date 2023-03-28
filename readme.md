# Analysis pipeline for the FRASER 2.0 manuscript

This repository contains the pipeline and scripts of the analysis described in 
the FRASER 2.0 manuscript. 

## Features

* Analysis on GTEx, comparing FRASER 2.0 to FRASER, LeafcutterMD and SPOT
* Validation on known pathgenic splice defects from a mitochondrial disease cohort

## Installation

To run this pipeline, please install FRASER 2.0 from its github first:
- `BiocManager::install("c-mertes/FRASER", ref="fraser2")`

Then install wBuild and snakemake using pip:

- `pip install wBuild`

Refer to the full documentation of wBuild in case of problems: https://wbuild.readthedocs.io

## Usage

* Clone this repository to an empty directory and navigate to it.
* Run FRASER (version 1), SPOT and LeafcutterMD separately first and update 
  the corresponding paths to their results in the wbuild.yaml file.
* Run `snakemake -n Paper_figures_all_figures_R` to run the full pipeline.

## FRASER 2.0

FRASER 2.0 is an open-source software. The source-code is available at https://github.com/gagneurlab/FRASER.

## Citation

The manuscript corresponding to this analysis is available at ...
