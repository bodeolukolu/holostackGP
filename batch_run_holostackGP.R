#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")

setwd("./")

holostackGP(
    projname=args[2],
    MTME=args[3],
    phenofile=args[4],
    genofile=args[5],
    metagenomefile=args[6],
    ploidy=2,
    traits=args[1],
    covariate=args[7],             # c("trait-9","traits-10")
    kernel=args[8],       # "metagenomic", genomic","holobiont", "metagenomic+genomic"
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    subsample_markers=NULL,       # if full data set, set to NULL
    gene_model <- args[9],             # "Additive", "Dominance", "metagenome or microbiome", Full
    R_libpath=NULL
)
