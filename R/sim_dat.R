

#' simulation Data.
#'
#' A simulated data using the same framework from simulation series I in the CIV paper.
#' It contains 500 subjects with one phenotype (X) of interest, one pleiotropic phenotype (Z) and
#' one outcome (Y) simulated for each subject. 9 SNPs were generated and they were associated with both X and Z.
#' The true causal effect of X on Y is 1. This data file is simulated to serve as an example
#' to estimate the causal effect of X on Y while accounting for
#' potential pleiotropic effect from Z. Users can compare the performance of different MR methods on this simulation dataset
#' since the true causal effect is known.
#' @docType data
#'
#' @usage data(simulation)
#' @param simulation$Y: the simulated outcome Y. Continuous variable.
#' @param simulation$X: The simulated phenotype of interest X. Continuous variable.
#' @param simulation$Z: The potential pleiotropic phenotype Z. Continuous variable.
#' @param simulation$G: The simulated genotypes. The dosage of 9 independent SNP variants
#' were simulated with a minor allele frequency of 0.3 for all 500 subjects.
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords dataset
#'
#' @references
#'
#' @examples
#' data(simulation)
#' X <- simulation$X
#' Z <- simulation$Z
#' geno <- simulation$G
#' outcome <- simulation$Y
"simulation"

