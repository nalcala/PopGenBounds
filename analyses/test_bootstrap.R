

# Stage 0: Setup -------------

devtools::load_all()          # run from the package root, so internal functions like boot_indices() are available
geno <- toad_genotypes        # or geno_saved
lv   <- microsat_allele_levels(geno)


# Stage 1: Resampling keeps population sizes and membership ------------------
# Each population should keep its number of toads, and every drawn toad should come from its own population.

set.seed(1)
idx <- boot_indices(geno$pop)

all(table(geno$pop[idx]) == table(geno$pop))   # TRUE: same size per population
length(idx) == nrow(geno)                       # TRUE

# Stage 2: Each toad is drawn once on average --------------------
# Over many resamples, every individual should appear about once per replicate.

set.seed(2)
draws  <- unlist(replicate(1000, boot_indices(geno$pop), simplify = FALSE))
counts <- tabulate(draws, nbins = nrow(geno)) / 1000
summary(counts)               # mean exactly 1, values roughly between 0.9 and 1.1

# Stage 3: Bootstrap allele frequencies are not shifted --------------
# Averaged over replicates, bootstrap frequencies should equal the original frequencies.

f0 <- microsat_freq(geno, lv)[[2]]                       # locus 1, original data
set.seed(3)
fb <- replicate(500, microsat_freq(geno[boot_indices(geno$pop), ], lv)[[2]],
                simplify = "array")
diffs <- apply(fb, c(1, 2), mean) - f0
summary(as.vector(diffs))     # centred on 0, small values of both signs

# Stage 4: The reported estimates are those of the original data -------------------
#The value columns must come from the original sample, not from the replicates.
set.seed(4)
res <- boot_stat(geno, stat_fun = Diff, B = 50)

# per locus: equal to Diff on the original frequency matrix
freq <- microsat_freq(geno, lv)
all.equal(res$estimate[[1]]$value, Diff(freq[[1]])$value[1:4])   # TRUE

# multilocus FST: equal to a ratio of averages computed by hand
#tabs <- lapply(freq, Diff)
#HS <- mean(sapply(tabs, function(d) d$value[d$statistic == "HS"]))
#HT <- mean(sapply(tabs, function(d) d$value[d$statistic == "HT"]))
#all.equal(res$mean$value[res$mean$statistic == "FST"], 1 - HS / HT)  # TRUE

# Stage 5: Results are reproducible -----------------
#The same seed must give exactly the same results.

set.seed(10); r1 <- boot_stat(geno, Diff, B = 30)
set.seed(10); r2 <- boot_stat(geno, Diff, B = 30)
identical(r1$mean, r2$mean)   # TRUE

# Stage 6: Basic and percentile intervals are reflections of each other --------------------
# With the same seed, both types use the same replicates, so their limits must be mirror images around the estimate,
#and the bias must be identical.

set.seed(11); rb <- boot_stat(geno, Diff, B = 50, type = "basic")
set.seed(11); rp <- boot_stat(geno, Diff, B = 50, type = "percentile")

all.equal(rb$mean$value_lower, 2 * rb$mean$value - rp$mean$value_upper)  # TRUE
all.equal(rb$mean$value_bias,  rp$mean$value_bias)                        # TRUE


# Stage 7: Increasing B stabilises the intervals ------------

set.seed(20); r_a <- boot_stat(geno, Diff, B = 1000)
set.seed(21); r_b <- boot_stat(geno, Diff, B = 1000)
cbind(r_a$mean[, c("statistic", "value_lower", "value_upper")],
      r_b$mean[, c("value_lower", "value_upper")])
# limits should agree to about the second or third decimal; if not, increase B




