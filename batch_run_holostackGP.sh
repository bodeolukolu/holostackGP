#!/bin/bash

cd ./
max_jobs=3     # each trait is run with 5 cores when model stacking is been run, so multiply each job by 5 cores.


main () {

declare -a arr=("trait-1" "trait-2")
projname=proj_GP
MTME=FALSE
phenofile=Traits.txt
genofile=geno.txt
metagenomefile=metagenome.txt
covariate=trait-10,trait-11
kernel=genomic       # "metagenomic", genomic","holobiont", "metagenomic+genomic"
gene_model=c("Full")         # "Additive", "Dominance", "metagenome or microbiome", "Full"


for trait in "${arr[@]}"; do
  Rscript holostackGP.R ${trait}" $projname $MTME $phenofile $genofile $metagenomefile $covariate $kernel $genomic_model & sleep 10
  ((count++))
  if (( count % max_jobs == 0 )); then
    wait
  fi
done
wait

}
main 2>> ./log.out
