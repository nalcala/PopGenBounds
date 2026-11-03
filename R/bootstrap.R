# ---------------------------------------------------------------
# bootstrap.R
# Bootstrap confidence intervals for per-locus statistics (Fst, G'st, D)
#
# Resampling unit: individuals, within each population, with replacement.
#
#   res <- boot_stat(toad_genotypes, stat_fun = Diff, B = 1000)
#   res$estimate[[1]]   # locus 1: statistic, value, bound, ..., value_lower, value_upper
# ---------------------------------------------------------------


#' Bootstrap confidence intervals for per-locus statistics
#'
#' Resamples individuals with replacement within each population (keeping
#' the original sample sizes), recomputes the allele frequency matrices and
#' applies `stat_fun` to the matrix of each locus. Percentile confidence
#' intervals are added to the estimates of the original data. The allele set
#' of the original data is used in every replicate, so alleles not drawn get
#' frequency 0.
#'
#' @param genotypes Table from [microsat_genotypes()].
#' @param stat_fun Function applied to one allele frequency matrix
#'   (populations x alleles), returning a tibble with one row per statistic
#'   (e.g. [Diff()]). Rows must come out in the same order every time.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level of the intervals.
#' @param type Interval type: `"percentile"` (default; the 2.5% and 97.5%
#'   quantiles of the bootstrap values, for `level = 0.95`) or `"basic"`
#'   (the same quantiles reflected around the estimate, which corrects for
#'   the bias estimated by the bootstrap).
#' @param value_col Name of the column of the `stat_fun` output that holds the
#'   statistic values (default `"value"`). Confidence intervals, bias and
#'   bias-corrected values are computed for this column#' @param label_col Column of the `stat_fun` output naming the statistics.
#' @param keep_boot If `TRUE`, also return the bootstrap values.
#' @param ... Further arguments passed to `stat_fun`.
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{List (one element per locus) of the `stat_fun` output
#'       for the original data, with added columns `<value_col>_lower` and
#'       `<value_col>_upper` (confidence limits), `<value_col>_bias` (mean of
#'       bootstrap replicates minus estimate) and `<value_col>_corrected`
#'       (estimate minus bias).}
#'     \item{mean}{Data frame with the mean over loci of each statistic for the
#'       original data, with the same added columns as `estimate`, computed
#'       from the mean over loci within each bootstrap replicate.}
#'     \item{boot}{If `keep_boot = TRUE`: list (one element per locus) of the
#'       bootstrap outputs stacked into one table with a `replicate` column;
#'       otherwise `NULL`.}
#'     \item{B}{Number of bootstrap replicates used.}
#'     \item{level}{Confidence level of the intervals.}
#'     \item{type}{Interval type used: `"percentile"` or `"basic"`.}
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(2021)
#' res <- boot_stat(toad_genotypes, stat_fun = Diff, B = 20)
#' res$estimate[[1]]
#' res$mean
#' }
boot_stat <- function(genotypes, stat_fun, B = 1000, level = 0.95, type = c("percentile","basic"),
                      value_col = "value", label_col = "statistic", keep_boot = FALSE, ...) {
  type <- match.arg(type)
  genotypes$pop <- droplevels(genotypes$pop)              # only populations present (e.g. pairwise subsets)
  lv <- microsat_allele_levels(genotypes)                 # fixed allele set (original data)
  freq_orig <- microsat_freq(genotypes, lv)
  estimate  <- lapply(freq_orig, stat_fun, ...)           # stats of original data, per locus

  boot <- lapply(seq_len(B), function(b) {                # bootstrap distribution of stats
    if (b %% 100 == 0) message("bootstrap ", b, "/", B)
    idx <- boot_indices(genotypes$pop)
    freq_boot <- microsat_freq(genotypes[idx, ], lv)
    lapply(freq_boot, stat_fun, ...)
  })

  # confidence intervals, locus by locus
  result <- lapply(seq_along(estimate), function(k) {
    boot_k <- lapply(boot, function(r) r[[k]])            # locus k from every replicate
    boot_ci_df(estimate[[k]], boot_k, level=level, type=type, col=value_col)
  })
  names(result) <- names(estimate)

  # mean over loci, within each replicate, and its confidence interval
  #K <- nlevels(genotypes$pop)
  mean_obs  <- mean_over_loci(estimate, value_col, label_col)                     # original data
  mean_boot <- lapply(boot, function(r) mean_over_loci(r, value_col, label_col))  # one table per replicate
  mean_tab  <- boot_ci_df(mean_obs, mean_boot, level=level, type=type, col=value_col)

  # optional: bootstrap values, stacked per locus
  boot_tables <- NULL
  if (keep_boot) {
    boot_tables <- lapply(seq_along(estimate), function(k)
      do.call(rbind, lapply(seq_len(B), function(b) cbind(replicate = b, boot[[b]][[k]]))))
    names(boot_tables) <- names(estimate)
  }

  list(estimate = result, mean=mean_tab, boot = boot_tables, B = B, level = level, type = type)
}

# Internal: mean over loci of the value column, row by row (one value per
# statistic). tabs = list of per-locus tables with the same rows in the same
# order. Returns a table with the label column and the mean values.
mean_over_loci <- function(tabs, col="value", label_col="statistic") {
  # Ratio of averages approach for computing mean differentiation stats
  #avg <- function(s) mean(vapply(tabs, function(d) d[[col]][d[[label_col]] == s], numeric(1)))
  #M  <- avg("M")
  #HS <- avg("HS")
  #HT <- avg("HT")
  #FST  <- 1 - HS / HT
  #GpST <- FST * (K - 1 + HS) / (K - 1) / (1 - HS)
  #D    <- K / (K - 1) * (HT - HS) / (1 - HS)

  #out <- data.frame(c("M", "FST", "G'ST", "D"))
  #names(out) <- label_col
  #out[[col]] <- c(M, FST, GpST, D)
  #out

  out <- data.frame(tabs[[1]][label_col])
  vals <- matrix(unlist(lapply(tabs, function(d) d[[col]])), nrow = nrow(tabs[[1]]))
  out[[col]] <- rowMeans(vals)
  out
}

# Internal: bootstrap confidence limits, bias and bias-corrected values.
# est = vector of estimates (one per row); vals = matrix rows x B of bootstrap values.
# Returns a data frame with one row per estimate.
boot_ci_core <- function(est, vals, level = 0.95, type = "percentile") {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  q_lo  <- apply(vals, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
  q_hi  <- apply(vals, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
  bias  <- rowMeans(vals, na.rm = TRUE) - est
  if (type == "percentile") {
    lower <- q_lo; upper <- q_hi
  } else {                                        # basic: reflected around the estimate
    lower <- 2 * est - q_hi; upper <- 2 * est - q_lo
  }
  data.frame(value_lower = lower, value_upper = upper,
             value_bias = bias, value_corrected = est - bias)
}

# Internal: add bootstrap results for the value column to the estimate table
# of one locus (or of the mean over loci): value_lower, value_upper,
# value_bias and value_corrected (see boot_ci_core()).
# boot_k = list of B bootstrap tables; same row order in every table.
boot_ci_df <- function(estimate, boot_k, level, type = "percentile", col = "value") {
  vals <- matrix(unlist(lapply(boot_k, function(d) d[[col]])), nrow = nrow(estimate))   # rows x B
  ci   <- boot_ci_core(estimate[[col]], vals, level, type)
  estimate[names(ci)] <- ci
  estimate
}

# Internal: row indices of one bootstrap sample, resampling individuals
# with replacement within each population (sample sizes unchanged).
boot_indices <- function(pop) {
  unlist(lapply(split(seq_along(pop), pop),
                function(ix) ix[sample.int(length(ix), replace = TRUE)]),
         use.names = FALSE)
}
