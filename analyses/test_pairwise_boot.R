#Test pairwise.R module

## Load package ---------------
library(PopGenBounds)

## Load data ----
toad_genotypes <- PopGenBounds::toad_genotypes
#Toad_microsat_freq = PopGenBounds::Toad_microsat_freq

## Allele frequencies ----
lv   <- microsat_allele_levels(toad_genotypes)
toad_frequencies <- microsat_freq(toad_genotypes, lv)
# plot frequency table
ggfreqtable(toad_frequencies[["F22"]])

## Parallel computation of HS, HT, M for each pair of subpopulations -----
K <- 2
f <- toad_frequencies[[1]]  # one locus of allele frequencies
s <- rowSums(f^2)  # sum of squared frequencies per population
# Compute HS from formula: Hs = 1 - S/K, with S = sum of squared allele frequencies over alleles and populations
HS <- 1 - outer(s, s, "+") / K  # HS for all pairs of populations (HS[i,j] = HS(pop_i,pop_j))

# Compute HT from formula: HT = 1 - sum of squared means of allele frequencies over populations.
HT <- 1 - (outer(s, s, "+") + 2 * tcrossprod(f)) / 4  # HT for all pairs of populations (HT[i,j] = HT(pop_i,pop_j))

# Compute M - maximal mean allele frequency over populations.
M <- matrix(0, nrow(f), nrow(f))
for (a in seq_len(ncol(f))) M <- pmax(M, outer(f[, a], f[, a], "+") / 2)  # M[i,j] = max_a((p_ia + p_ja)/2)

pairwise_H <- list(HS = HS, HT = HT, M = M)

## Check with loop computation of HS, HT and M
pairwise_H_loop <- function(f) {
  K <- nrow(f)
  HS <- HT <- M <- matrix(NA, K, K)
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    f2 <- f[c(i, j), ]
    HS[i, j] <- HS[j, i] <- mean(1 - rowSums(f2^2))
    HT[i, j] <- HT[j, i] <- 1 - sum(colMeans(f2)^2)
    M[i, j]  <- M[j, i]  <- max(colMeans(f2))
  }
  list(HS = HS, HT = HT, M = M)
}

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

# check time efficiency
f <- microsat_freq(toad_genotypes)[[1]]
system.time(for (r in 1:100) pairwise_H(f))
system.time(for (r in 1:100) pairwise_H_loop(f))

#check values
a <- pairwise_H(f); b <- pairwise_H_loop(f)
up <- upper.tri(a$HS)
all.equal(a$HS[up], b$HS[up])   # TRUE
all.equal(a$HT[up], b$HT[up])   # TRUE
all.equal(a$M[up],  b$M[up])    # TRUE

## Diff for one locus ----

# M, FST, G'ST and D for all pairs at one locus (Diff() formulas with K = 2).
# Returns a list of four population x population matrices.
pairwise_locus <- function(f) {
  h   <- pairwise_H(f)
  FST <- ifelse(h$HT < 1e-12, 0, 1 - h$HS / h$HT)
  list(M      = h$M,
       FST    = FST,
       "G'ST" = FST * (1 + h$HS) / (1 - h$HS),
       D      = 2 * (h$HT - h$HS) / (1 - h$HS))
}

f <- toad_frequencies[[1]]  # one locus of allele frequencies
h <- pairwise_H(f)
head(h$HS)

FST <- 1 - h$HS / h$HT
head(FST)

pw_diff <- pairwise_locus(f)

# Test with simple matrices
f_fix <- matrix(c(0.0, 1.0, 0.0,     # P1: fixed for a2
                  0.2, 0.5, 0.3,     # P2: polymorphic
                  1.0, 0.0, 0.0),    # P3: fixed for a1
                nrow = 3, byrow = TRUE,
                dimnames = list(c("P1", "P2", "P3"), c("a1", "a2", "a3")))


# compute HT, HS, M
h <- pairwise_H(f_fix)

# compute FST, G'ST, D
pw_diff <- pairwise_locus(f_fix)
pw_diff

# automated test:
r <- pairwise_locus(f_fix)
h <- pairwise_H(f_fix)

round(h$HS, 4)        # off-diagonal: [1,2] 0.31, [1,3] 0, [2,3] 0.31
round(h$HT, 4)        # [1,2] 0.405, [1,3] 0.5, [2,3] 0.555
r$M                   # [1,2] 0.75, [1,3] 0.5, [2,3] 0.6
round(r$FST, 4)       # [1,2] 0.2346, [1,3] 1, [2,3] 0.4414
round(r$`G'ST`, 4)    # [1,2] 0.4453, [1,3] 1, [2,3] 0.8381
round(r$D, 4)         # [1,2] 0.2754, [1,3] 1, [2,3] 0.7101

# each pair against Diff
for (p in list(c(1, 2), c(1, 3), c(2, 3)))
  print(all.equal(Diff(f_fix[p, ])$value,
                  c(r$M[p[1], p[2]], r$FST[p[1], p[2]], r$`G'ST`[p[1], p[2]], r$D[p[1], p[2]])))


f2 <- matrix(c(0.5, 0.5, 0.0,
               0.0, 0.5, 0.5),
             nrow = 2, byrow = TRUE,
             dimnames = list(c("P1", "P2"), c("a1", "a2", "a3")))

h <- pairwise_H(f2);  r <- pairwise_locus(f2)
c(HS = h$HS[1, 2], HT = h$HT[1, 2], M = r$M[1, 2],
  FST = r$FST[1, 2], GpST = r$`G'ST`[1, 2], D = r$D[1, 2])
# expected: HS 0.5, HT 0.625, M 0.5, FST 0.2, GpST 0.6, D 0.5

Diff(f2)$value     # M, FST, G'ST, D from Diff: 0.5 0.2 0.6 0.5

## One locus diff matrix transform ----

up <- which(upper.tri(diag(3)), arr.ind = TRUE)     # array of index pairs: each pair of pop once
pops <- rownames(f_fix)
pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])

matrix(sapply(pw_diff, function(m) m[up]), nrow = nrow(up))     # pairs x statistics

## Stack diff for loci
pairwise_values <- function(freq_list, up) {
  per_locus <- lapply(freq_list, function(f) {
    s <- pairwise_locus(f)
    matrix(sapply(s, function(m) m[up]), nrow = nrow(up))     # pairs x statistics
  })
  arr  <- array(unlist(per_locus), dim = c(nrow(up), 4, length(freq_list)))
  mean <- matrix(apply(arr, c(1, 2), mean, na.rm = TRUE), nrow = nrow(up))
  list(per_locus = arr, mean = mean)
}

pops  <- rownames(toad_frequencies[[1]])
up    <- which(upper.tri(diag(length(pops))), arr.ind = TRUE)     # array of index pairs: each pair of pop once
pairs <- data.frame(pop1 = pops[up[, 1]], pop2 = pops[up[, 2]])

per_locus <- lapply(toad_frequencies, function(f) {
  s <- pairwise_locus(f)
  matrix(sapply(s, function(m) m[up]), nrow = nrow(up))     # pairs x statistics
})
arr  <- array(unlist(per_locus), dim = c(nrow(up), 4, length(toad_frequencies)))
mean <- matrix(apply(arr, c(1, 2), mean, na.rm = TRUE), nrow = nrow(up))

## Wrap in long tables ----

v <- pairwise_values(toad_frequencies, up)
estimate <- lapply(seq_along(toad_frequencies), function(l)
  ci_table(matrix(v$per_locus[, , l], nrow = nrow(up)), NULL, pairs))

names(estimate) <- names(toad_frequencies)

mean = ci_table(v$mean, NULL, pairs)
mean

## Test pairwise.R module -------
# compute diff statistics for each pair of populations from original data
pw  <- pairwise_stats(toad_frequencies)   # observed values

d <- pw$estimate[["5F"]]
table(d$statistic)                               # 1081 rows each for M, FST, G'ST, D
d[d$pop1 == "WUN" & d$pop2 == "LIE", ]

# run small bootstrap and compute bootstrap_CI for each pair of statistics
set.seed(1)
res <- boot_pairwise(toad_genotypes, B = 20)          # with CIs
names(res)

pts <- do.call(rbind, res$estimate)
all.equal(mean(pts$value[pts$statistic == "FST"], na.rm = TRUE),
          res$overall$value[res$overall$statistic == "FST"])   # TRUE

### Run full bootstrap and plot -----

set.seed(2021)
res_pw <- boot_pairwise(toad_genotypes, B = 1000)

# save results
out_dir <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/output"
saveRDS(res_pw, file.path(out_dir, "pairwise_bootstrap_B1000.rds"))

# check overall mean
pts <- do.call(rbind, res_pw$estimate)
mean(pts$value[pts$statistic == "FST"], na.rm = TRUE)    # mean over the 10,810 plotted points
res_pw$overall$value[res_pw$overall$statistic == "FST"]  # should be identical

sum(is.nan(pts$value))                                   # how many missing values exist at all

# summary of the results
f <- res_pw$mean[res_pw$mean$statistic == "FST", ]
summary(f$value)                          # range of pairwise FST
summary(f$value_upper - f$value_lower)    # typical CI width
mean(f$value_lower > 0)                   # share of pairs whose CI excludes 0

# a caterpillar plot
f <- f[order(f$value), ]
f$rank <- seq_len(nrow(f))
ggplot(f, aes(x = rank, y = value, ymin = value_lower, ymax = value_upper)) +
  geom_linerange(colour = "grey60") +
  geom_point(size = 0.4) +
  labs(x = "Population pairs (ranked by FST)", y = expression(italic(F[ST]))) +
  theme_bw()

caterpillar <- function(tab, stat, y_label) {
  f <- tab[tab$statistic == stat, ]
  f <- f[order(f$value), ]
  f$rank <- seq_len(nrow(f))
  ggplot(f, aes(x = rank, y = value, ymin = value_lower, ymax = value_upper)) +
    geom_linerange(colour = "grey60") +
    geom_point(size = 0.4) +
    labs(x = paste("Population pairs (ranked by", stat, ")"), y = y_label) +
    theme_bw()
}

p_fst  <- caterpillar(res_pw$mean, "FST",  expression(italic(F[ST])))
p_gpst <- caterpillar(res_pw$mean, "G'ST", expression(italic(G*minute[ST])))
p_d    <- caterpillar(res_pw$mean, "D",    expression(italic(D)))

p_fst; p_gpst; p_d
p_fst + p_gpst + p_d


# plot pairwise stats
pts <- do.call(rbind, res_pw$estimate)                 # all loci and pairs
get <- function(s) pts$value[pts$statistic == s]

ggtoad_pw = ggbounds_raw(M = get("M"), FST = get("FST"), GpST = get("G'ST"), D = get("D"),
             K = 2, ci = res_pw$overall, density_colors = TRUE)

ggtoad_pw[[1]] + ggtoad_pw[[2]] + ggtoad_pw[[3]]
## Analysis of bias ----

# global bootstrap (K=47)
set.seed(2021)
res_global <- boot_stat(toad_genotypes, stat_fun = Diff, B = 1000)

# The ratio compares the bias with the half-width of the confidence interval
ratio <- function(tab) with(tab, value_bias / ((value_upper - value_lower) / 2))

# global analysis (K = 47)
cbind(res$mean[, c("statistic", "value_bias")], ratio = ratio(res$mean))

# pairwise analysis: summaries over the 1,081 pairs
pw <- res_pw$mean
tapply(pw$value_bias, pw$statistic, median)          # typical pairwise bias
tapply(ratio(pw), pw$statistic, median)              # typical bias / half-width

## Check effciency of bootstrap commands -----
# efficiency checks
profvis::profvis(boot_pairwise(toad_genotypes, B = 50))

set.seed(2021)
system.time(res_test <- boot_pairwise(toad_genotypes, B = 50, locus_ci = TRUE))
