#!/bin/bash

set -u pipefail

export OMP_NUM_THREADS=5
export MKL_NUM_THREADS=5
export OPENBLAS_NUM_THREADS=5

max_jobs=3     # each trait is run with 5 cores when model stacking is been run, so multiply each job by 5 cores.
count=0

declare -a arr=("trait-1" "trait-2")
projname="proj_GP"
phenofile="Traits.txt"
genofile="geno.txt"
metagenomefile="metagenome.txt"
# set covariate to NULL, if none
covariate="trait-10,trait-11"
# metagenomic, genomic, holobiont, metagenomic+genomic (select only one at a time or create a nested loop below)
kernel="genomic"
# Additive,Dominance,metagenome or microbiome,Full (you can specify multiple e.g "Additive","Dominance")
gene_model="Full"

# Function to count running background jobs
running_jobs() {
    jobs -pr | wc -l
}

for trait in "${arr[@]}"; do
    while (( $(running_jobs) >= max_jobs )); do
        sleep 300
    done
    echo "Launching genomic prediction for trait: ${trait}"
    Rscript batch_run_holostackGP.R \
    "${trait}" "$projname" "$phenofile" "$genofile" \
    "$metagenomefile" "$covariate" "$kernel" "$gene_model" >> "log_${trait}.out" 2>&1 &
    sleep 10
done

# Wait for all remaining jobs to finish
wait
echo "All jobs completed."
