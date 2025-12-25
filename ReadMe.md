
# Introduction
This repository provides a flexible and extensible R framework for MTME-inspired genomic prediction, integrating multiple statistical and machine-learning models with multi-kernel (multi-omic, full epistatic interactions and GxExM) data. It is designed for research in quantitative genetics, plant and animal breeding, multi-omic modeling, and genomic-enabled predictions.

The pipeline implements a broad suite of genomic prediction models (GBLUP, rrBLUP, RKHS, BayesA, BayesB, BayesC, BayesL, and BRR) and combines them using a stacked ensemble. By integrating models with distinct assumptions about marker effects (linear, nonlinear, shrinkage-based, or mixture priors), the stacking approach leverages complementary patterns across learners, producing more robust and often more accurate predictions. It also models additive and non-additive effects, including dominance, epistatic, and genome-by-metagenome interacts to improve predictive ability. Multi-trait and multi-environment (MTME)-inspired genomic prediction is implemented by model stacking (bayesian, elastic net, lasso, linear and ridge stackigng methods)

In addition, the framework supports multi-kernel integration, allowing users to incorporate genomic relationship matrices, microbiome/metagenome kernels, and other omic-based kernels. Any of the kernels can be run independently or in combination.

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
Download the <run_holostackGP.R>, edit parameters and run the pipeline. Multiple traits can be run with the run_holostackGP.R or batch_run_holostackGP.R. Additive and Dominance models can be submitted together, while metagenome/microbiome and Full (stacking of all models, kernels, gene models) models should be run independently.
On unix/Linux, the analysis can also be submitted as a batch file using batch_run_holstackGP.sh and batch_run_holostackGP.R.


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
