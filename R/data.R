# R/data.R
# Documentation for datasets included in easyLSEA.

#' Example lipidomics dataset
#'
#' A synthetic dataset of 200 lipid species simulating a case vs control
#' lipidomics comparison, with known enrichment patterns built in:
#' PC and PE species are enriched in the case group, TG species are
#' depleted. Used in package examples and tests.
#'
#' @format A \code{data.frame} with 200 rows and 6 columns:
#' \describe{
#'   \item{LipidName}{Character. Lipid identifier in shorthand notation
#'     (e.g. "PC 36:2").}
#'   \item{LipidClass}{Character. Pre-assigned lipid class abbreviation.}
#'   \item{logFC}{Numeric. Log2 fold change (case / control).}
#'   \item{P.Value}{Numeric. Raw p-value from simulated differential analysis.}
#'   \item{adj.P.Val}{Numeric. Benjamini-Hochberg adjusted p-value.}
#'   \item{sig}{Integer. 1 if adj.P.Val < 0.05 and |logFC| > log2(1.25),
#'     0 otherwise.}
#' }
#'
#' @source Simulated data. See \code{data-raw/lipid_example.R} for the
#'   generation script. Seed: 2026.
"lipid_example"
