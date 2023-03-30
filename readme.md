# Analysis pipeline for the FRASER 2.0 manuscript

This is the accompanying analysis repository of the paper:

`Improved detection of aberrant splicing using the Intron Jaccard Index`. 

The paper can be found [on medRxiv](todo once submitted).

This repository contains the full pipeline and code to reproduce the results 
using [snakemake](https://snakemake.readthedocs.io/en/stable/) and [wBuild](https://github.com/gagneurlab/wBuild). 

## Features

* Analysis on GTEx, comparing FRASER 2.0 to FRASER, LeafcutterMD and SPOT
* Analysis of amount of splicing outliers on two rare disease datasets
* Validation on known pathgenic splice defects from a mitochondrial disease cohort

## Installation

To run this pipeline, please install FRASER 2.0 first from its github first:
- `BiocManager::install("c-mertes/FRASER", ref="fraser2")`

Then install wBuild and snakemake using pip:

- `pip install wBuild`

Refer to the full documentation of wBuild in case of problems: https://wbuild.readthedocs.io

## Usage

* Clone this repository to an empty directory and navigate to it.
* Run FRASER (version 1), SPOT and LeafcutterMD separately first and update 
  the corresponding paths to their results in the wbuild.yaml file. This pipeline 
  assumes that those results have already been generated.
* Run `snakemake -n Paper_figures_all_figures_R` to run the full pipeline.

## FRASER 2.0

FRASER 2.0 is an improved version of FRASER and an open-source software. 
The source-code is available at https://github.com/gagneurlab/FRASER.

## Citation

The manuscript corresponding to this analysis is available at [medRxiv link].
