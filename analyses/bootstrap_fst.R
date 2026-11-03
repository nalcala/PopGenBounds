# ---------------------------------------------------------------
# Bootstrap CIs for Fst (global and pairwise)
# Resampling unit: individuals, within each population, with replacement
# Fst estimator: Nei's G_ST (multilocus: sum of numerators / sum of denominators)
# ---------------------------------------------------------------
library(readxl)

pkg_dir <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/PopGenBounds"
toad <- read_xlsx(file.path(pkg_dir, "data-raw", "ProhlEtAl2021_Bombina variegata Microsat data.xlsx"))

# ---- Prepare genotypes ------------------------------------------
pop  <- factor(toad[[2]], levels = unique(toad[[2]]))     # populations in file order
geno <- as.data.frame(toad[, 4:23])                        # 10 loci x 2 allele columns
geno[geno == -9] <- NA
n_loci     <- ncol(geno) / 2
locus_names <- colnames(toad)[seq(4, 23, 2)]

# Fixed allele set per locus (from the ORIGINAL data), so every
# bootstrap replicate gives matrices with identical columns
allele_levels <- lapply(seq_len(n_loci), function(i)
  sort(unique(na.omit(c(geno[[2 * i - 1]], geno[[2 * i]])))))

# ---- Genotypes -> list of allele frequency matrices -------------
freq_matrices <- function(geno, pop, allele_levels) {
  out <- vector("list", n_loci)
  for (i in 1:n_loci) {
    a   <- c(geno[[2 * i - 1]], geno[[2 * i]])            # both gene copies
    p   <- factor(rep(as.character(pop), 2))
    cnt <- unclass(table(p, factor(a, levels = allele_levels[[i]])))  # NA dropped
    out[[i]] <- cnt / rowSums(cnt)                         # p_i = n_i / 2N
  }
  names(out) <- locus_names
  out
}

# ---- Nei's G_ST from a list of frequency matrices ---------------
fst_nei <- function(freq_list) {
  num <- 0; den <- 0
  for (m in freq_list) {
    m <- m[stats::complete.cases(m), , drop = FALSE]      # drop pops with no data at this locus
    if (nrow(m) < 2) next
    Hs   <- mean(1 - rowSums(m^2))                         # mean within-pop heterozygosity
    Ht   <- 1 - sum(colMeans(m)^2)                         # total heterozygosity
    num  <- num + (Ht - Hs)
    den  <- den + Ht
  }
  if (den == 0) return(NA_real_)                           # all loci monomorphic
  num / den
}

# ---- Pairwise Fst matrix (reuses the per-population frequencies) -
pairwise_fst <- function(freq_list) {
  pops <- rownames(freq_list[[1]]); k <- length(pops)
  M <- matrix(0, k, k, dimnames = list(pops, pops))
  for (a in 1:(k - 1)) for (b in (a + 1):k) {
    M[a, b] <- M[b, a] <- fst_nei(lapply(freq_list, function(m) m[c(a, b), , drop = FALSE]))
  }
  M
}

# ---- One bootstrap sample: resample individuals within each pop --
boot_indices <- function(pop) {
  unlist(lapply(split(seq_along(pop), pop),
                function(ix) ix[sample.int(length(ix), replace = TRUE)]),
         use.names = FALSE)
}

# ---- Observed values --------------------------------------------
freq_obs     <- freq_matrices(geno, pop, allele_levels)
fst_obs      <- fst_nei(freq_obs)
pw_fst_obs   <- pairwise_fst(freq_obs)

# ---- Bootstrap --------------------------------------------------
set.seed(2021)
B <- 100
k <- nlevels(pop)
fst_boot    <- numeric(B)
pw_fst_boot <- array(NA_real_, dim = c(k, k, B),
                     dimnames = list(levels(pop), levels(pop), NULL))

for (b in seq_len(B)) {
  idx  <- boot_indices(pop)
  fl_b <- freq_matrices(geno[idx, , drop = FALSE], pop[idx], allele_levels)
  fst_boot[b]       <- fst_nei(fl_b)
  #pw_fst_boot[, , b] <- pairwise_fst(fl_b)
  if (b %% 100 == 0) message("bootstrap ", b, "/", B)
}

# ---- Percentile confidence intervals ----------------------------
ci_global <- quantile(fst_boot, c(0.025, 0.975), na.rm = TRUE)
cat("Global Fst =", round(fst_obs, 4),
    " 95% CI [", round(ci_global[1], 4), ",", round(ci_global[2], 4), "]\n")

pw_ci_low  <- apply(pw_fst_boot, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
pw_ci_high <- apply(pw_fst_boot, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)

# Tidy table of pairwise results
pairs <- which(upper.tri(pw_fst_obs), arr.ind = TRUE)
pw_table <- data.frame(pop1 = rownames(pw_fst_obs)[pairs[, 1]],
                       pop2 = colnames(pw_fst_obs)[pairs[, 2]],
                       fst  = pw_fst_obs[pairs],
                       ci_low  = pw_ci_low[pairs],
                       ci_high = pw_ci_high[pairs])
head(pw_table)
