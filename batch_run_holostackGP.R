#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 9)

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")
setwd("./")

covariate <- if (is.na(args[6]) || args[6] == "" || tolower(args[6]) == "null") {
  NULL
} else {
  trimws(strsplit(args[6], ",", fixed = TRUE)[[1]])
}
gene_model <- if (is.na(args[8]) || args[8] == "" || tolower(args[8]) == "null") {
  NULL
} else {
  trimws(strsplit(args[9], ",", fixed = TRUE)[[1]])
}

holostackGP(
    projname=args[2],
    phenofile=args[3],
    genofile=args[4],
    metagenomefile=args[5],
    ploidy=2,
    traits=args[1],
    # for covariates: c("trait-9","traits-10")
    covariate=covariate,
    # "metagenomic", genomic","holobiont", "metagenomic+genomic"
    kernel=args[7],
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    # if full marker ata set, set to NULL
    subsample_markers=NULL,
    # "Additive", "Dominance", "metagenome or microbiome", "Full"
    gene_model=gene_model,
    R_libpath=NULL
)
