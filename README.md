# Synopsis

In this package we encode Constrained Instrumental Variable methods, including constrianed instrumental variable (CIV) and constrianed instrumental variable with smoothed L0 penalty (CIV_smooth), to construct valid instrumental variables and perform adjusted causal effect estimation when pleiotropic exists, focusing particularly on the situation where pleiotropic phenotypes have been measured. We also include the functions to implement some existing methods (2SLS and Allele score methods) and their variants to account for pleiotropy. More details of these methods can be found in the paper "Constrained Instruments and its Application to Mendelian Randomization with Pleiotropy".

# Code Example
 
Two datasets are included: simulation and ADNI. A detailed walk-through are provided in the vignette (CIVMR/vignette/CIVMR.pdf) to show an example of MR analysis using the functions in CIVMR.

# Installation
 
Open an R session:
 
library(devtools) ; install_github("LaiJiang/CIVMR")

# Developer
 
Lai Jiang, Department of Epidemiology, Biostatistics and Occupational Health, McGill University, 1020
Pine Avenue West, Montreal, Quebec, H3A 1A2, Canada

lai.jiang@mail.mcgill.ca
