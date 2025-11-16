#!/usr/bin/env Rscript

source("https://raw.githubusercontent.com/bodeolukolu/holostackGP/refs/heads/master/holostackGP.R")

setwd("./")

holostackGP(
    projname="FER",
    phenofile="Traits_FER_lsmeans.txt",
    genofile="Zm_2x_rd6_maf0.02_dose.txt",
    metagenomefile="WGS_metagenome.txt",
    ploidy=2,
    traits=c("FER"),
    covariate=c("DTA_21","Aflatoxin"),             # c("trait-9","traits-10")
    kernel=c("genomic"),       # "metagenomic or microbiome", genomic","holobiont", "metagenomic+genomic or microbiome+genomic"
    CVrep=100,
    k_fold=5,
    maf=0.02,
    geno_missing_rate=0.2,
    subsample_markers=NULL,       # if full data set, set to NULL
    gene_model=c("Full"),           # "Additive", "Dominance", "metagenome or microbiome", Full
    R_libpath=NULL
)
