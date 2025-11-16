#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")

setwd("./")

holostackGP(
    projname=args[2],
    phenofile=args[3],
    genofile=args[4],
    metagenomefile=args[5],
    ploidy=2,
    traits=args[1],
    covariate=args[6],             # c("trait-9","traits-10")
    kernel=args[7],       # "metagenomic", genomic","holobiont", "metagenomic+genomic"
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    subsample_markers=NULL,       # if full data set, set to NULL
    gene_model=c("Full"),         # "Additive", "Dominance", "metagenome or microbiome", Full
    R_libpath=NULL
)
