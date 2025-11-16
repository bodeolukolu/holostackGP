
# Introduction
This repository provides a flexible and extensible R framework for genomic prediction, integrating multiple statistical and machine-learning models with multi-kernel (multi-omic) data. It is designed for research in quantitative genetics, plant and animal breeding, multi-omic modeling, and genomic-enabled predictions.

The pipeline implements a broad suite of genomic prediction models (GBLUP, rrBLUP, RKHS, BayesA, BayesB, BayesC, BayesL, and BRR) and combines them using a stacked ensemble. By integrating models with distinct assumptions about marker effects (linear, nonlinear, shrinkage-based, or mixture priors), the stacking approach leverages complementary patterns across learners, producing more robust and often more accurate predictions. It also models additiv and non-additive effects, including dominance, epistatic, and genome-by-metagenome interacts to improve predictive ability.

In addition, the framework supports multi-kernel integration, allowing users to incorporate genomic relationship matrices, microbiome/metagenome kernels, environmental similarity matrices, and other omic-based kernels. Any of the kernels can be run independently or in combination.

For questions, bugs, and suggestions, please contact bolukolu@gmail.com.

## Features
Implements a diverse set of prediction algorithms:
GBLUP
rrBLUP
RKHS (Reproducing Kernel Hilbert Space)
BayesA
BayesB
BayesC
BayesL / Lasso-type Bayesian model
BRR (Bayesian Ridge Regression)
Each model can be run independently or included in the stacked ensemble.

## How to run holostackGP
Download the <run_holostackGP.R>, edit parameters and run the pipeline. To submit a batch file for multiple traits and parameters, use the batch_run_holstackGP.sh (implemented with bash). Multiple traits and gene models (Additive, Dominance, metagenome or microbiome and Full or stacking of all models, kernels, gene models) can also be run with the run_holostackGP.R


## Dependencies
- mice
- data.table
- dplyr
- AGHmatrix
- vegan
- compositions
- ggcorrplot
- nlme
- lsmeans
- agricolae
- doParallel
- foreach
- parallel
- rrBLUP
- BGLR


## Select Article Referencing holostackGP


## Acknowledgment


## Troubleshooting

## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/bodeolukolu/holostackGP/blob/master/LICENSE">Apache License Version 2.0</a>
