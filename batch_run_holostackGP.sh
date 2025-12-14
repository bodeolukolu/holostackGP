#!/bin/bash
set -euo pipefail

max_jobs=3     # each trait is run with 5 cores when model stacking is been run, so multiply each job by 5 cores.
count=0

main () {

declare -a arr=("trait-1" "trait-2")
projname=proj_GP"
MTME="FALSE"
phenofile="Traits.txt"
genofile="geno.txt"
metagenomefile="metagenome.txt"
# set covariate to NULL, if none
covariate="trait-10","trait-11"
# metagenomic, genomic, holobiont, metagenomic+genomic (select only one at a time or create a nested loop below)
kernel="genomic"
# Additive,Dominance,metagenome or microbiome,Full (you can specify multiple e.g "Additive","Dominance")
gene_model="Full"


for trait in "${arr[@]}"; do
  Rscript holostackGP.R ${trait} $projname $MTME $phenofile $genofile $metagenomefile $covariate $kernel $gene_model & sleep 10
  ((count++))
  if (( count % max_jobs == 0 )); then
    wait
  fi
done
wait

}
main 2>> ./log.out
