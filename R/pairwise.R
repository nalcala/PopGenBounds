# ---------------------------------------------------------------
# pairwise.R
# Pairwise differentiation statistics (K = 2) for all population pairs:
# per locus, and their arithmetic mean over loci, with bootstrap CIs.
# Output mirrors boot_stat(): res$estimate (one table per locus) and res$mean.
#
#   pw  <- pairwise_stats(microsat_freq(toad_genotypes))   # observed values
#   res <- boot_pairwise(toad_genotypes, B = 1000)          # with CIs
# ---------------------------------------------------------------


#' Pairwise differentiation statistics
#'
#' For every pair of populations, computes \eqn{M}, \eqn{F_{ST}},
#' \eqn{G'_{ST}} and \eqn{D} at each locus with the same formulas as [Diff()]
#' applied to the two rows of the pair (\eqn{K = 2}), and their arithmetic
#' mean over loci. Loci where a statistic is undefined for a pair (e.g. both
#' populations fixed for the same allele) are skipped in the mean.
#'
#' @param freq_list List of allele frequency matrices (populations x alleles),
#'   one per locus, as returned by [microsat_freq()].
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of tables with columns
#'       `pop1`, `pop2`, `statistic` and `value`, one row per pair and statistic.}
#'     \item{mean}{Table with the same columns: mean over loci of each statistic.}
#'     \item{overall}{Table with columns `statistic` and `value`: mean of each
#'       statistic over all pairs and loci.}
#'   }
#' @export
pairwise_stats <- function(freq_list) {
  pops  <- rownames(freq_list[[1]])
  up    <- which(upper.tri(diag(length(pops))), arr.ind = TRUE)     # array of index pairs: each pair of pop once
  pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])
  v     <- pairwise_values(freq_list, up)   # returns 'per_locus' (3-dim array pairs x stats x loci) and 'mean' (pairs x stats array of means over loci)

  estimate <- lapply(seq_along(freq_list), function(l)
    pairs_table(matrix(v$per_locus[, , l], nrow = nrow(up)), pairs=pairs))
  names(estimate) <- names(freq_list)

  overall <- data.frame(statistic = c("M", "FST", "G'ST", "D"),
                        value     = apply(v$per_locus, 2, mean, na.rm = TRUE))   # mean over all pairs and loci

  list(estimate = estimate, mean = pairs_table(v$mean, pairs=pairs), overall=overall)
}


#' Bootstrap confidence intervals for pairwise statistics
#'
#' Resamples individuals with replacement within each population (as in
#' [boot_stat()]); each replicate resamples all populations once, and every
#' pair is computed from it. Confidence intervals, bias and bias-corrected
#' values are computed for the per-locus statistics and for their mean over
#' loci.
#'
#' @param genotypes Table from [microsat_genotypes()].
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param type Interval type: `"percentile"` (default; the 2.5% and 97.5%
#'   quantiles of the bootstrap values, for `level = 0.95`) or `"basic"`
#'   (the same quantiles reflected around the estimate, which corrects for
#'   the bias estimated by the bootstrap).
#' @param locus_ci If `TRUE` (default), compute intervals for each locus as
#'   well as for the mean. `FALSE` gives per-locus values without intervals
#'   and uses much less memory.
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of tables with columns
#'       `pop1`, `pop2`, `statistic`, `value` (original data) and, if
#'       `locus_ci = TRUE`, `value_lower`, `value_upper`, `value_bias` and
#'       `value_corrected`.}
#'     \item{mean}{Same table for the mean over loci, always with intervals.}
#'     \item{overall}{Mean of each statistic over all pairs and loci, with
#'       `value_lower`, `value_upper`, `value_bias` and `value_corrected`,
#'       computed from the same mean within each bootstrap replicate.}
#'     \item{B}{Number of bootstrap replicates used.}
#'     \item{level}{Confidence level of the intervals.}
#'     \item{type}{Interval type used.}
#'   }
#' @export
boot_pairwise <- function(genotypes, B = 1000, level = 0.95,
                          type = c("percentile","basic"), locus_ci = TRUE) {
  type <- match.arg(type)
  genotypes$pop <- droplevels(genotypes$pop)
  lv    <- microsat_allele_levels(genotypes)
  freq  <- microsat_freq(genotypes, lv)
  pops  <- rownames(freq[[1]])
  up    <- which(upper.tri(diag(length(pops))), arr.ind = TRUE)     # each pair once
  pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])
  P <- nrow(up)                                                     # number of pairs
  L <- length(lv)                                                   # number of loci
  stat_names <- c("M", "FST", "G'ST", "D")

  obs <- pairwise_values(freq, up)                                  # original data (2 matrices: per-locus stats:
                                                                    # [P x (stats) x L] and mean stats [P x (stats)])
  overall_obs <- apply(obs$per_locus, 2, mean, na.rm = TRUE)        # mean over all pairs and loci
  # bootstrap values: pairs x statistics x (loci) x replicates
  boot_mean <- array(NA_real_, dim = c(P, 4, B))
  boot_overall_mean <- matrix(NA_real_, nrow = 4, ncol = B)              # overall means over pairs and loci (within bootstrap replicates)
  boot_loc  <- if (locus_ci) array(NA_real_, dim = c(P, 4, L, B)) else NULL

  for (b in seq_len(B)) {
    if (b %% 50 == 0) message("bootstrap ", b, "/", B)
    idx <- boot_indices(genotypes$pop)
    v   <- pairwise_values(microsat_freq(genotypes[idx, ], lv), up)
    boot_mean[, , b] <- v$mean
    boot_overall_mean[, b] <- apply(v$per_locus, 2, mean, na.rm = TRUE)
    if (locus_ci) boot_loc[, , , b] <- v$per_locus
  }
  # per-locus tables with CI
  estimate <- lapply(seq_len(L), function(l) {
    est <- matrix(obs$per_locus[, , l], nrow = P)
    ci  <- if (locus_ci) pairs_ci(est, array(boot_loc[, , l, ], dim = c(P, 4, B)), level=level, type=type) else NULL
    pairs_table(est, pairs=pairs, ci=ci)
  })
  names(estimate) <- names(lv)

  # per-pair means over loci with Ci
  mean_tab <- pairs_table(obs$mean, pairs = pairs,
                          ci = pairs_ci(obs$mean, boot_mean, level = level, type = type))

  # overall mean over all pairs and loci with CI
  overall_mean <- cbind(data.frame(statistic = stat_names, value = overall_obs),
                   boot_ci_core(overall_obs, boot_overall_mean, level = level, type = type))

  list(estimate = estimate, mean = mean_tab, overall = overall_mean,
       B = B, level = level, type = type)
}

# ---- internal helpers ---------------------------------------------


# Internal: HS, HT and M for all pairs of populations at one locus.
# f = populations x alleles frequency matrix. For pair (i, j):
#   HS = mean of the two within-population heterozygosities
#   HT = 1 - sum of squared mean frequencies of the pair
#   M  = largest mean frequency of the pair
pairwise_H <- function(f) {
  s  <- rowSums(f^2)                                          # sum of squared freqs per population
  # Compute HS from formula: Hs = 1 - S/K, with S = sum of squared frequencies over alleles and populations
  # (for each pair of populations: HS[i,j] = HS(pop_i,pop_j))).
  HS <- 1 - outer(s, s, "+") / 2
  #(for each pair of populations: HT[i,j] = HT(pop_i,pop_j))).
  HT <- 1 - (outer(s, s, "+") + 2 * tcrossprod(f)) / 4        # sum(((f_i + f_j) / 2)^2), all pairs at once.
  # Compute M - maximal mean allele frequency over populations.
  M  <- matrix(0, nrow(f), nrow(f))
  for (a in seq_len(ncol(f))) M <- pmax(M, outer(f[, a], f[, a], "+") / 2)
  list(HS = HS, HT = HT, M = M)
}

# M, FST, G'ST and D for all pairs at one locus (Diff() formulas with K = 2).
# Returns a list of four population x population matrices.
pairwise_locus <- function(f) {
  h   <- pairwise_H(f)
  mono <- h$HT < 1e-12                       # both populations fixed for the same allele (M = 1)
  FST <- ifelse(h$HT < 1e-12, 0, 1 - h$HS / h$HT)   # set FST to 0 (instead of NaN) if HT=0
  out <- list(M      = h$M,
              FST    = FST,
              "G'ST" = FST * (1 + h$HS) / (1 - h$HS),
              D      = 2 * (h$HT - h$HS) / (1 - h$HS))
  lapply(out, function(m) { m[mono] <- NA; m })   # exclude these combinations everywhere
}

# Per-locus values and their mean over loci, for the pairs in `up`.
# Returns per_locus (pairs x 4 statistics x loci) and mean (pairs x 4 statistics).
pairwise_values <- function(freq_list, up) {
  per_locus <- lapply(freq_list, function(f) {
    s <- pairwise_locus(f)
    matrix(sapply(s, function(m) m[up]), nrow = nrow(up))     # pairs x statistics
  })
  arr  <- array(unlist(per_locus), dim = c(nrow(up), 4, length(freq_list)))   # 3-dimensional array, with loci as the third dim
  avg <- matrix(apply(arr, c(1, 2), mean, na.rm = TRUE), nrow = nrow(up))
  list(per_locus = arr, mean = avg)
}

# Internal: long table (one row per pair and statistic).
# est = pairs x 4 statistics; ci = optional list of 4 data frames from boot_ci_core().
pairs_table <- function(est, pairs, ci = NULL) {
  stat_names <- c("M", "FST", "G'ST", "D")
  out <- lapply(seq_along(stat_names), function(k) {
    tab <- data.frame(pairs, statistic = stat_names[k], value = est[, k])
    if (!is.null(ci)) tab <- cbind(tab, ci[[k]])
    tab
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Internal: intervals for each statistic. est = pairs x 4; vals = pairs x 4 x B.
pairs_ci <- function(est, vals, level, type) {
  P <- nrow(est)
  lapply(seq_len(ncol(est)), function(k)
    boot_ci_core(est[, k], matrix(vals[, k, ], nrow = P), level, type))
}
