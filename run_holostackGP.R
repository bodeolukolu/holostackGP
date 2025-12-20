#!/usr/bin/env Rscript

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")

holostackGP(
    projname="FER",
    phenofile="Traits_FER_lsmeans.txt",
    genofile="Zm_2x_rd6_maf0.02_dose.txt",
    metagenomefile="WGS_metagenome.txt",
    ploidy=2,
    traits=c("FER"),
    covariate=c("DTA_21","Aflatoxin"),
    # "metagenomic", genomic","holobiont", "metagenomic+genomic"
    kernel="genomic",
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    # if full marker data set, set to NULL
    subsample_markers=NULL,
    # Additive, Dominance, metagenome or microbiome, Full
    gene_model=c("Full"),
    R_libpath=NULL
)
