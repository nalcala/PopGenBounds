#' Produces a list of allele frequency matrices for each locus, for biallelic loci.
#'
#' This function transforms a dataframe of allele frequencies in subpopulations for multiple loci
#' into a list of allele frequency matrices, one per locus, compatible with the downstream computation
#' of differentiation statistics (the `Diff()` function).
#'
#' @param data_frequencies A data frame of allele frequencies with
#'   dimensions \eqn{L \times K}, where \eqn{L} is the number of loci
#'   (rows) and \eqn{K} is the number of sub-populations (columns).
#'   Each entry corresponds to the frequency of the first allele for
#'   that locus. The second allele frequency is computed as
#'   \eqn{q = 1 - p}.
#'
#' @return A list of \eqn{K \times 2} matrices of allele frequencies,
#'   one per locus. Each matrix contains the frequencies of the
#'   reference and alternative alleles for each subpopulation.
#' @export
#'
#' @examples
#' # Create a toy dataframe of allele frequencies for 3 loci and 2 sub-populations
#' df <- data.frame(
#'   Subpop1 = c(0.8, 0.6, 0.4),
#'   Subpop2 = c(0.7, 0.5, 0.3)
#' )
#' # Use the function to generate a list of 2x2 frequency matrices
#' freq_list <- make_popgen_input(df)
make_popgen_input <- function(data_frequencies) {
  list_loci <- lapply(1:nrow(data_frequencies), function(i){
    mat <- as.matrix(data_frequencies[i,])
    second_allele <- 1 - mat[1, ]  # Compute 1-x for each element in the first row
    new_mat <- rbind(mat, second_allele)
    rownames(new_mat) <- c("First", "Second")
    return(t(new_mat))
  }
  )
  return(list_loci)
}
