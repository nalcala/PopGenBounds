# Test microsat.R module

## Load packages -----------------------------------------------
library(readxl)

## Load microsat data ------------------------------------------

pkg_dir <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/PopGenBounds"
toad <- read_xlsx(file.path(pkg_dir, "data-raw", "ProhlEtAl2021_Bombina variegata Microsat data.xlsx"))

# load genotype data from pacakge
geno_saved <- PopGenBounds::toad_genotypes

## Testing microsat functions ----------------------------------

geno <- microsat_genotypes(toad)
lv   <- microsat_allele_levels(geno_saved)
freq <- microsat_freq(geno_saved, lv)

# plot frequency matrices

ggfreqtable(freq[[1]])


# Compute diff statistics

Diff_toad_test = lapply(freq,Diff)

## Test bootsrap -----------------------------------------------

# Test bootstrap.R module
set.seed(2021)
res <- boot_stat(geno_saved, stat_fun = Diff, B = 10)

# Results
res$estimate[["5F"]]      # one locus: statistic, value, bound, ..., value_lower, value_upper
res$estimate              # all loci

res$mean

# Save
saveRDS(res, "../output/bootstrap_toad.rds")


# Test bootstrap functions in bootstrap.R
# Internal: add <col>_lower / <col>_upper (percentile limits) to the
# estimate table of one locus. boot_k = list of B bootstrap tables of that
# locus; intervals are computed row by row (same row order in every table).
boot_ci_df <- function(estimate, boot_k, level, ci_cols) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  for (col in ci_cols) {
    # rows x B matrix: column b holds replicate b
    vals <- matrix(unlist(lapply(boot_k, function(d) d[[col]])), nrow = nrow(estimate))
    estimate[[paste0(col, "_lower")]] <- apply(vals, 1, stats::quantile, probs = probs[1], na.rm = TRUE)
    estimate[[paste0(col, "_upper")]] <- apply(vals, 1, stats::quantile, probs = probs[2], na.rm = TRUE)
  }
  estimate
}


# Internal: row indices of one bootstrap sample, resampling individuals
# with replacement within each population (sample sizes unchanged).
boot_indices <- function(pop) {
  unlist(lapply(split(seq_along(pop), pop),
                function(ix) ix[sample.int(length(ix), replace = TRUE)]),
         use.names = FALSE)
}
# tabs = list of per-locus tables with the same rows in the same order.
mean_over_loci <- function(tabs, ci_cols, label_col) {
  out <- data.frame(tabs[[1]][label_col])
  for (col in ci_cols) {
    vals <- matrix(unlist(lapply(tabs, function(d) d[[col]])), nrow = nrow(tabs[[1]]))
    out[[col]] <- rowMeans(vals)
  }
  out
}

# Compute bootstrap distribution
stat_fun <- Diff

lv <- microsat_allele_levels(geno_saved)              # fixed allele set (original data)
freq_orig <- microsat_freq(geno_saved, lv)
estimate <- lapply(freq_orig, stat_fun)   # stats of original data

B <- 200
boot <- lapply(seq_len(B), function(b) {            # bootstrap distribution of stats
  idx <- boot_indices(geno_saved$pop)
  freq_boot <- microsat_freq(geno_saved[idx, ], lv)
  lapply(freq_boot, stat_fun)
})

# confidence intervals, locus by locus
result <- lapply(seq_along(estimate), function(k) {
  boot_k <- lapply(boot, function(r) r[[k]])            # locus k from every replicate
  boot_ci_df(estimate[[k]], boot_k, level=0.95, ci_cols="value")
})
names(result) <- names(estimate)

# mean over loci, within each replicate, and its confidence interval
mean_obs  <- mean_over_loci(estimate, ci_cols="value", label_col="statistic")                     # original data
mean_boot <- lapply(boot, function(r) mean_over_loci(r, ci_cols="value", label_col="statistic"))  # one table per replicate
mean_tab  <- boot_ci_df(mean_obs, mean_boot, level=0.95, ci_cols="value")

estimate_ci <- lapply(seq_along(estimate), function(k)
  boot_ci_df(estimate[[k]], lapply(boot, function(r) r[[k]]), level = 0.95, ci_cols = "value"))
names(estimate_ci) <- names(estimate)

#test mean over loci function
tabs <- boot[[1]]
col <- "value"
vals <- matrix(unlist(lapply(tabs, function(d) d[[col]])), nrow = nrow(tabs[[1]]))

# test of confidence intervals, locus by locus

boot_1 <- lapply(boot, function(r) r[[1]])   # locus 1 from every replicate

result <- lapply(names(estimate), function(l) {
  boot_l <- lapply(boot, function(r) r[[l]])     # locus l (e.g. "5F") from every replicate
  boot_ci_df(estimate[[l]], boot_l, level, ci_cols)
})
names(result) <- names(estimate)

## Test horizontal stacking of bootstrap stat replicates to compute confidence intervals
vals <- matrix(unlist(lapply(boot_1, function(d) d[["value"]])), nrow = nrow(estimate[[1]]))

## Check boot_l content (should be B replicates for locus l)

loci <- names(estimate)
result     <- vector("list", length(loci))   # empty list, one slot per locus
boot_saved <- vector("list", length(loci))
names(result) <- names(boot_saved) <- loci

for (l in loci) {
  cat(l)
  boot_l <- lapply(boot, function(r) r[[l]])        # locus l from every replicate
  boot_saved[[l]] <- boot_l                          # save it
  result[[l]] <- boot_ci_df(estimate[[l]], boot_l, level=0.95, ci_cols="value")
}
