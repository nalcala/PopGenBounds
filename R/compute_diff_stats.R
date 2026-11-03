
#' Compute differentiation statistics \eqn{F_{ST}}, \eqn{G'_{ST}}, and \eqn{D}
#' for a list of data frames of allele frequencies or subclonal fractions.
#' The differentiation statistics are computed across all available subpopulations.
#'
#' @param data_subclonal A list of data frames containing allele frequencies or subclonal fractions.
#' @param col_indices Indices of columns corresponding to frequencies in subpopulations.
#'
#' @return A named list containing frequency of the most frequent allele \eqn{M}, \eqn{F_{ST}}, \eqn{G'_{ST}}, and \eqn{D}
#'   computed for each locus in each dataset.
#' @export
compute_diff_list <- function(data_subclonal, col_indices) {
  result <- lapply(names(data_subclonal), function(name) {
    df <- data_subclonal[[name]]
    list_freq <- make_popgen_input(df[, col_indices])
    lapply(list_freq, Diff)
  })
  names(result) <- names(data_subclonal)
  return(result)
}

#' Compute differentiation statistics \eqn{F_{ST}}, \eqn{G'_{ST}}, and \eqn{D} for K = 3 subpopulations.
#'
#' @param list_loci A list of \eqn{K \times 2} matrices of allele frequencies, one per locus.
#'   Each matrix contains the frequencies of the reference and alternative alleles (columns)
#'   for each subpopulation (rows).
#'
#' @return A list of four elements:
#'   - The first element is a list of matrices containing differentiation statistics
#'     \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D}, and a vector of frequencies \eqn{M}
#'     computed across all three sub-populations.
#'   - The remaining three elements are lists of matrices containing pairwise differentiation
#'     statistics for each pair of sub-populations out of the three.
#' @export
#'
#' @examples
#' # Toy example with 2 loci and 3 subpopulations
#' locus1 <- matrix(c(0.8, 0.2,   # Subpop 1
#'                    0.6, 0.4,   # Subpop 2
#'                    0.7, 0.3),  # Subpop 3
#'                  nrow = 3, byrow = TRUE)
#'
#' locus2 <- matrix(c(0.5, 0.5,
#'                    0.4, 0.6,
#'                    0.6, 0.4),
#'                  nrow = 3, byrow = TRUE)
#'
#' list_loci <- list(locus1, locus2)
#'
#' # Compute differentiation statistics
#' diff_stats <- compute_Diff_3subpop(list_loci)
#'
#' # Access 3-subpop differentiation stats for first locus
#' print(diff_stats$D_123[[1]])
#'
#' # Access pairwise differentiation (subpop 1 vs 2) for second locus
#' print(diff_stats$D_12[[2]])
compute_Diff_3subpop <- function(list_loci){
  # check whether K=3 (to do)

  # Extract couple of rows (subpops) from 3-subpop matrices
  list_loci_12 <- lapply(list_loci, function(mat) mat[c(1, 2), ])
  list_loci_23 <- lapply(list_loci, function(mat) mat[c(2, 3), ])
  list_loci_13 <- lapply(list_loci, function(mat) mat[c(1, 3), ])

  # Compute stats (F_ST, G_ST, D)
  D_loci_123 = lapply(list_loci, Diff)
  D_loci_12 = lapply(list_loci_12, Diff)
  D_loci_23 = lapply(list_loci_23, Diff)
  D_loci_13 = lapply(list_loci_13, Diff)

  return(list(D_123=D_loci_123, D_12=D_loci_12, D_13=D_loci_13, D_23=D_loci_23))

}
