#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")

setwd("./")

covariate <- if (is.na(args[7]) || args[7] == "" || tolower(args[7]) == "null") {
  NULL
} else {
  trimws(strsplit(args[7], ",", fixed = TRUE)[[1]])
}

gene_model <- if (is.na(args[9]) || args[9] == "" || tolower(args[9]) == "null") {
  NULL
} else {
  trimws(strsplit(args[9], ",", fixed = TRUE)[[1]])
}

holostackGP(
    projname=args[2],
    MTME=args[3],
    phenofile=args[4],
    genofile=args[5],
    metagenomefile=args[6],
    ploidy=2,
    traits=args[1],
    # for covariates: c("trait-9","traits-10")
    covariate=args[7],
    # "metagenomic", genomic","holobiont", "metagenomic+genomic"
    kernel=args[8],
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    # if full marker ata set, set to NULL
    subsample_markers=NULL,
    # "Additive", "Dominance", "metagenome or microbiome", "Full"
    gene_model = args[9],
    R_libpath=NULL
)
