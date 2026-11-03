# ---------------------------------------------------------------
# ccf_bootstrap.R
# Bootstrap of sequencing read sampling for biallelic subclonal-fraction
# (CCF) data -- simplest version (step 1): depth and CCF only.
#
# Each locus (row): CCF in K sub-populations and the read depth at that
# position in each sub-population. In each replicate, new CCFs are drawn as
#     CCF* = Binomial(DP, CCF) / DP
# and the differentiation statistics are recomputed.
#
#   res <- boot_ccf(data, ccf_cols, dp_cols, B = 1000)
#   res$estimate   # per locus: statistic, value, intervals
#   res$mean       # mean over loci, with intervals
#
# Requires boot_ci_core() from bootstrap.R.
# ---------------------------------------------------------------


#' Differentiation statistics for biallelic loci (vectorised)
#'
#' Same formulas as [Diff()] applied to the matrix `cbind(p, 1 - p)` of each
#' locus, computed for all loci at once.
#'
#' @param P Matrix (loci x sub-populations) of frequencies of the first allele
#'   (e.g. subclonal fractions).
#' @param exclude_fixed If `TRUE` (default), loci without variation (mean
#'   frequency 0 or 1, i.e. M = 1) get `NA` for all statistics.
#'
#' @return Matrix (loci x 4) with columns `M`, `FST`, `G'ST` and `D`.
#' @export
diff_biallelic <- function(P, exclude_fixed = TRUE) {
  K    <- ncol(P)
  pbar <- rowMeans(P)
  HS   <- rowMeans(2 * P * (1 - P))          # = mean over sub-populations of 1 - p^2 - (1 - p)^2
  HT   <- 2 * pbar * (1 - pbar)              # = 1 - pbar^2 - (1 - pbar)^2
  M    <- pmax(pbar, 1 - pbar)               # = max(colMeans(cbind(p, 1 - p)))
  FST  <- ifelse(HT < 1e-12, 0, 1 - HS / HT)
  out  <- cbind(M      = M,
                FST    = FST,
                "G'ST" = FST * (K - 1 + HS) / (K - 1) / (1 - HS),
                D      = K / (K - 1) * (HT - HS) / (1 - HS))
  if (exclude_fixed) out[HT < 1e-12, ] <- NA
  out
}


#' Bootstrap of read sampling for subclonal fractions (simple model)
#'
#' Parametric bootstrap treating the sequencing depth as the sample size of
#' each CCF: in each replicate, CCF* = Binomial(DP, CCF) / DP, independently
#' for each locus and sub-population. Purity and copy number are not modelled,
#' so the intervals are a lower bound on the read-sampling uncertainty.
#'
#' @param data Data frame, one row per locus.
#' @param ccf_cols Columns of `data` with the subclonal fractions.
#' @param dp_cols Columns of `data` with the read depths, in the same order.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param type `"percentile"` (default) or `"basic"` intervals.
#'
#' @return A list with
#'   \describe{
#'     \item{estimate}{Table with one row per locus and statistic: `locus`,
#'       `statistic`, `value` (original data), `value_lower` and `value_upper`
#'       (confidence limits), `value_bias` (mean of bootstrap values minus
#'       `value`) and `value_corrected` (`value` minus bias).}
#'     \item{mean}{Mean over loci of each statistic, with intervals computed
#'       from the mean over loci within each replicate.}
#'     \item{B}{Number of bootstrap replicates used.}
#'     \item{level}{Confidence level of the intervals.}
#'     \item{type}{Interval type used.}
#'   }
#' @export
boot_ccf <- function(data, ccf_cols, dp_cols, B = 1000, level = 0.95,
                     type = c("percentile", "basic")) {
  type <- match.arg(type)
  CCF  <- as.matrix(data[, ccf_cols])
  DP   <- as.matrix(data[, dp_cols])
  if (!identical(dim(CCF), dim(DP))) stop("'ccf_cols' and 'dp_cols' must select the same number of columns.")
  L <- nrow(CCF)         # number of loci

  # original data
  obs <- diff_biallelic(CCF)                                   # loci x 4
  stat_names <- colnames(obs)

  # bootstrap: loci x statistics x replicates
  boot_loc <- array(NA_real_, dim = c(L, 4, B))
  for (b in seq_len(B)) {
    if (b %% 100 == 0) message("bootstrap ", b, "/", B)
    boot_loc[, , b] <- diff_biallelic(ccf_resample(CCF, DP))
  }

  # ---- per-locus results: one table per statistic, then stacked ----
  boot_table_per_stat <- function(k) {
    # observed value of statistic k at each locus (vector of length L)
    observed_k <- obs[, k]
    # bootstrap values of statistic k: one row per locus, one column per replicate (L x B)
    bootstrap_k <- matrix(boot_loc[, k, ], nrow = L)
    # identification of each row and its observed value
    locus_info <- data.frame(locus     = seq_len(L),
                             statistic = stat_names[k],
                             value     = observed_k)
    # confidence limits, bias and bias-corrected value for each locus
    boot_info <- boot_ci_core(est = observed_k, vals = bootstrap_k,
                              level = 0.95, type = "percentile")
    # both parts have L rows in the same locus order: put them side by side
    cbind(locus_info, boot_info)
  }

  boot_tables <- lapply(seq_along(stat_names), boot_table_per_stat)

  # stack the four tables (M, FST, G'ST, D) into one long table of 4 x L rows
  estimate <- do.call(rbind, boot_tables)
  rownames(estimate) <- NULL

  # ---- mean over loci, with its interval ----
  # observed: mean of each statistic over loci (vector of length 4)
  mean_observed <- colMeans(obs, na.rm = TRUE)
  # bootstrap: mean over loci within each replicate (4 statistics x B replicates)
  mean_one_replicate <- function(b) colMeans(boot_loc[, , b], na.rm = TRUE)
  mean_bootstrap     <- matrix(unlist(lapply(seq_len(B), mean_one_replicate)),
                               nrow = length(stat_names))
  # observed means and their confidence limits, side by side
  mean_left  <- data.frame(statistic = stat_names, value = mean_observed)
  mean_right <- boot_ci_core(est = mean_observed, vals = mean_bootstrap,
                             level = 0.95, type = "percentile")
  mean_boot_tab   <- cbind(mean_left, mean_right)

  list(estimate = estimate, mean = mean_boot_tab, B = B, level = level, type = type)
}


# Internal: one bootstrap draw of CCFs, CCF* = Binomial(DP, CCF) / DP.
# CCF, DP: loci x sub-populations matrices.
ccf_resample <- function(CCF, DP) {
  X <- matrix(stats::rbinom(length(CCF), size = DP, prob = CCF), nrow = nrow(CCF))
  X / DP
}
