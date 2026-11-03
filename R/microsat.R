#' Handles data conversion from raw microsatellite genotypes table into allele frequency matrices for each locus
# ---------------------------------------------------------------
# microsat.R
# Microsatellite genotype tables -> allele frequency matrices
#
# Standard genotype table (output of microsat_genotypes()):
#   column 1 = id, column 2 = pop (factor, file order),
#   columns 3, 4 = locus 1; columns 5, 6 = locus 2; ...
#   i.e. locus i is in columns 2*i + 1 and 2*i + 2.
#   Locus names = original header of the first column of each pair.
#
#   geno <- microsat_genotypes(raw_table)
#   lv   <- microsat_allele_levels(geno)
#   freq <- microsat_freq(geno, lv)
# ---------------------------------------------------------------


#' Prepare a microsatellite genotype table
#'
#' @param data Data frame, one row per individual, two columns per locus.
#' @param pop_col Column holding the population label.
#' @param id_col Column holding the individual ID.
#' @param allele_cols Allele columns (two consecutive columns per locus).
#' @param missing Code for missing alleles.
#' @return Data frame: id, pop, then two allele columns per locus.
#' @export
microsat_genotypes <- function(data, pop_col = 2, id_col = 1,
                               allele_cols = 4:ncol(data), missing = -9) {
  data <- as.data.frame(data)
  alleles <- data[allele_cols]
  alleles[alleles == missing] <- NA

  data.frame(id  = as.character(data[[id_col]]),
             pop = factor(data[[pop_col]], levels = unique(data[[pop_col]])),
             alleles, check.names = FALSE)
}

#' Allele set per locus
#'
#' All distinct alleles observed at each locus, sorted by size. Pass the set
#' from the original data to `microsat_freq()` so that resampled data give
#' matrices with identical columns.
#'
#' @param genotypes Table from `microsat_genotypes()`.
#' @return Named list of sorted allele vectors, one per locus.
#' @export
microsat_allele_levels <- function(genotypes) {
  n_loci <- (ncol(genotypes) - 2) / 2
  out <- lapply(seq_len(n_loci), function(i)
    sort(unique(c(genotypes[[2 * i + 1]], genotypes[[2 * i + 2]]))))   # sort drops NA
  names(out) <- names(genotypes)[seq(3, ncol(genotypes), by = 2)]
  out
}

#' Allele frequency matrices
#'
#' For each locus, the frequency of each allele in each population,
#' p_i = n_i / 2N, with missing alleles exclude,
#' with n_i - count of allele i, N - number of sampled individual in subpopultion.
#'
#' @param genotypes Table from `microsat_genotypes()`.
#' @param allele_levels Allele set per locus from `microsat_allele_levels()`.
#' @return Named list of matrices (populations x alleles), rows summing to 1.
#' @export
microsat_freq <- function(genotypes, allele_levels = microsat_allele_levels(genotypes)) {
  pop2 <- factor(rep(as.character(genotypes$pop), 2), levels = levels(genotypes$pop))

  out <- lapply(seq_along(allele_levels), function(i) {
    a   <- c(genotypes[[2 * i + 1]], genotypes[[2 * i + 2]])           # both gene copies
    cnt <- unclass(table(pop2, factor(a, levels = allele_levels[[i]])))  # NA dropped
    cnt / rowSums(cnt)                                               # p_i = n_i / 2N
  })
  names(out) <- names(allele_levels)
  out
}

