

## Load package ---------------
library(PopGenBounds)

## Load data ----
dir_data <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data"
organoids_ccf_data <- get(load(file.path(dir_data, "Organoids_CCF_clusters_dataset_fracs_list.RData")))
load(file.path(dir_data, "Organoids_CNV_tables.RData"))

data_3pop <- organoids_ccf_data[c("LCNEC4", "PANEC1")]
data_2pop <- organoids_ccf_data[c("SINET8", "SINET9", "LNET6", "LNET10", "LCNEC3")]
head(data_3pop[[1]])

data_2pop <- unparse_columns(data_2pop)
data_3pop <- unparse_columns(data_3pop)

# example of CCF_clonality data for one experiment
df_LCNEC4 <- read.delim(file.path(dir_data, "small_variants_CCFs/LCNEC4_annotatedvariants_CCF_clonality.tsv"))
head(df_LCNEC4)
dim(df_LCNEC4)

df_LNET10 <- read.delim(file.path(dir_data, "small_variants_CCFs/LNET10_annotatedvariants_CCF_clonality.tsv"))
head(df_LNET10)
dim(df_LNET10)
length(unique(df_LNET10$pos))

# filtered data on LNET10
df_LNET10_filt <- read.delim(file.path(dir_data, "small_variants_CCFs/clean_data/LNET10T.tsv"))
head(df_LNET10_filt)
dim(df_LNET10_filt)
length(unique(df_LNET10_filt$pos))

df_LCNEC3 <- read.delim(file.path(dir_data, "small_variants_CCFs/LCNEC3_annotatedvariants_clonality.tsv"))
head(df_LCNEC3)
dim(df_LCNEC3)

df_wgs <- read.delim(file.path(dir_data, "small_variants_WGS_PDTOs_parents.tsv"))
head(df_wgs)
dim(df_wgs)

df_rnaseq <- read.delim(file.path(dir_data, "small_variants_drivers_RNAseq_PDTOs_parents.tsv"))
head(df_rnaseq)
dim(df_rnaseq)

# compare different data sources -----
new_LNET10 <- organoids_ccf_data[["LNET10"]]
old_LNET10 <- df_LNET10_filt

dim(new_LNET10); dim(old_LNET10)

key <- c("chr", "pos", "ref", "alt")


# extract CCF and DP column names -----
d <- filter_clonal(data_3pop[[1]], K=3)
ccf_cols <- grep("^subclonal.fractions.", names(d), value = TRUE)
dp_cols  <- sub("^subclonal.fractions.", "DP.", ccf_cols)   # same sample order as the CCF columns
ccf_cols; dp_cols                                       # check the pairing

P <- as.matrix(d[, ccf_cols])

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

K    <- ncol(P)
pbar <- rowMeans(P)
HS   <- rowMeans(2 * P * (1 - P))          # = mean over sub-populations of 1 - p^2 - (1 - p)^2
HT   <- 2 * pbar * (1 - pbar)              # = 1 - pbar^2 - (1 - pbar)^2
M    <- pmax(pbar, 1 - pbar)               # = max(colMeans(cbind(p, 1 - p)))
FST  <- ifelse(HT < 1e-12, 0, 1 - HS / HT)

## check against Diff() ----
new <- diff_biallelic(as.matrix(d[, ccf_cols]), exclude_fixed = FALSE)
old <- t(sapply(lapply(make_popgen_input(d[, ccf_cols]), Diff), function(x) x$value[1:4]))

all.equal(unname(old), unname(new))

## run bootstrap -----

df <- organoids_ccf_data[["LCNEC4"]]

#unparse column names first
LCNEC4 <- do.call(data.frame, c(as.list(df), check.names = FALSE))

#filter data
LCNEC4_clean <- filter_clonal(LCNEC4, K=3, filter_fixed=T)

# select CCF and DP columns
ccf_cols <- grep("^subclonal.fractions.", names(LCNEC4_clean), value = TRUE)
dp_cols  <- sub("^subclonal.fractions.", "DP.", ccf_cols)   # same sample order as the CCF columns

CCF  <- as.matrix(LCNEC4_clean[, ccf_cols])
DP   <- as.matrix(LCNEC4_clean[, dp_cols])

identical(dim(CCF), dim(DP))

L <- nrow(CCF)

# original data
obs <- diff_biallelic(CCF)                                   # loci x 4
stat_names <- colnames(obs)

# check that the plot is consistent with previous results
gg <- ggbounds_raw(M = obs[, "M"], FST = obs[, "FST"], GpST = obs[, "G'ST"], D = obs[, "D"],
                   K = ncol(CCF))
gg[[1]] + gg[[2]] + gg[[3]]

# draw new CCF from binomial distribution
ccf_resample <- function(CCF, DP) {
  X <- matrix(stats::rbinom(length(CCF), size = DP, prob = CCF), nrow = nrow(CCF))
  X / DP
}

# bootstrap: loci x statistics x replicates
B <- 10
boot_loc <- array(NA_real_, dim = c(L, 4, B))
for (b in seq_len(B)) {
  if (b %% 100 == 0) message("bootstrap ", b, "/", B)
  boot_loc[, , b] <- diff_biallelic(ccf_resample(CCF, DP))
}

# compute CI for bootstrap distributions
estimate <- do.call(rbind, lapply(seq_along(stat_names), function(k)
  cbind(data.frame(locus = seq_len(L), statistic = stat_names[k], value = obs[, k]),
        boot_ci_core(obs[, k], matrix(boot_loc[, k, ], nrow = L), level = level, type = type))))

## long version


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

boot_table <- lapply(seq_along(stat_names), boot_table_per_stat)

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

# Run full bootstrap --------------
res <- boot_ccf(LCNEC4_clean, ccf_cols, dp_cols, B = 100)

# plot
res$estimate





