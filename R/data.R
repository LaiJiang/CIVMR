#' ADNI Data.
#'
#' A dataset collected from ADNI project for MR analysis in the CIV paper.
#' It contains 491 subjects, whose phenotypes (\eqn{A \beta}, Ptau, Ttau, Glucose levels) and
#' Alzheimer's status were collected. 20 SNPs associated with \eqn{A \beta} were selected and
#' their dosage of the 491 subjects were recorded. This data file is extracted to serve as an example
#' to estimate the causal effect of \eqn{A \beta} on progression of Alzheimer's disease while accounting for
#' potential pleiotropic effect from other phenotypes.
#' @docType data
#'
#' @usage data(ADNI)
#' @param ADNI$Y: the Alzheimer's disease status. Continuous variable. The raw status is binary variable,
#' and we adjusted it for confounding factors such as sex, age, education ... etc.
#' @param ADNI$X: The phenotype of interest \eqn{A \beta} Continuous variable.
#' @param ADNI$Z: The potential pleiotropic phenotypes (Ptau, Ttau, Glucose levels). Continuous variables.
#' @param ADNI$G: genotypes. The adjusted dosage of 20 SNPs.
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords dataset
#'
#' @references
#'
#' @examples
#' data(ADNI)
#' X <- ADNI$X
#' Z <- ADNI$Z
#' G <- ADNI$G
#' Y <- ADNI$Y
"ADNI"
