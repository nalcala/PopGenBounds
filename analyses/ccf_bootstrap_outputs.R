# ---------------------------------------------------------------
# ccf_bootstrap_outputs.R
# Figures and tables from the CCF bootstrap (step 1 model).
#
# Requires, from the previous steps:
#   prep     named list (one element per experiment) from prepare_ccf():
#            $data, $ccf_cols, $dp_cols, $K
#   res_ccf  named list (same names) of boot_ccf() results:
#            $estimate, $mean, $B, $level, $type
#   out_dir  output folder
# ---------------------------------------------------------------
devtools::load_all()
library(ggplot2)
library(patchwork)
library(writexl)          # install.packages("writexl") if needed

experiments <- names(res_ccf)
stats_diff  <- c("FST", "G'ST", "D")
stat_labels <- list("FST"  = expression(bolditalic(F[ST])),
                    "G'ST" = expression(bolditalic(G*minute[ST])),
                    "D"    = expression(bolditalic(D)))


# ===============================================================
# 1. Caterpillar plots: one row per experiment, one column per statistic
# ===============================================================
caterpillar_ccf <- function(res, stat, y_label, ylim = c(0, 1)) {
  f <- res$estimate[res$estimate$statistic == stat & !is.na(res$estimate$value), ]
  f <- f[order(f$value), ]
  f$rank <- seq_len(nrow(f))
  m <- res$mean[res$mean$statistic == stat, ]

  ggplot(f, aes(x = rank, y = value, ymin = value_lower, ymax = value_upper)) +
    geom_linerange(colour = "grey70", linewidth = 0.2) +
    geom_point(size = 0.2) +
    geom_hline(yintercept = m$value, colour = "red", linetype = "dashed", linewidth = 0.4) +
    labs(x = "Loci (ranked)", y = y_label) +
    coord_cartesian(ylim = ylim) +
    theme_bw()
}

# panels in row order: all FST panels, then all G'ST panels, then all D panels
panels <- unlist(lapply(stats_diff, function(s)
  lapply(experiments, function(e) caterpillar_ccf(res_ccf[[e]], s, stat_labels[[s]]))),
  recursive = FALSE)

# experiment names as column titles, on the top row only
n <- length(experiments)
for (j in seq_len(n)) {
  e <- experiments[j]
  panels[[j]] <- panels[[j]] +
    ggtitle(e) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}

fig_caterpillar <- wrap_plots(panels, ncol = n) +
  plot_layout(axes = "collect") &
  scale_y_continuous(labels = function(x) format(x, drop0trailing = TRUE))
fig_caterpillar

ggsave(file.path(dir_out, "figures/Fig_caterpillar_organoids.pdf"), fig_caterpillar,
       width = 1.6 * n, height = 1.6 * 3)
# ===============================================================
# 2. Bounds plots (ggbounds_raw): one row per experiment
# ===============================================================
# one ggbounds_raw() call per experiment: returns list(FST panel, G'ST panel, D panel)
bounds_list <- lapply(experiments, function(e) {
  est <- res_ccf[[e]]$estimate
  get <- function(s) est$value[est$statistic == s]
  ggbounds_raw(M = get("M"), FST = get("FST"), GpST = get("G'ST"), D = get("D"),
               K = prep[[e]]$K, ci = res_ccf[[e]]$mean,
               ci_style = "cross", mean_style = "crosshair",
               ci_linewidth = 0.5, ci_cap = 0.03,
               density_colors = TRUE, show_colorbar = FALSE)
})
names(bounds_list) <- experiments

# experiment names as column titles, on the top row (FST panels) only
for (e in experiments) {
  bounds_list[[e]][[1]] <- bounds_list[[e]][[1]] +
    ggtitle(e) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}

# panels in row order: all FST panels, then all G'ST panels, then all D panels
panels <- c(lapply(experiments, function(e) bounds_list[[e]][[1]]),
            lapply(experiments, function(e) bounds_list[[e]][[2]]),
            lapply(experiments, function(e) bounds_list[[e]][[3]]))

drop_mean_point <- function(p) {
  is_point <- vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
  p$layers[[max(which(is_point))]] <- NULL
  p
}
panels <- lapply(panels, drop_mean_point)

fig_bounds <- wrap_plots(panels, ncol = length(experiments)) +
  plot_layout(guides = "collect", axes = "collect")
fig_bounds <- fig_bounds &
  scale_x_continuous(labels = function(x) format(x, drop0trailing = TRUE))&
  scale_y_continuous(labels = function(x) format(x, drop0trailing = TRUE))

fig_bounds

ggsave(file.path(dir_out, "figures/Supp_Fig_organoids_CI_cross.pdf"), fig_bounds,
       width = 1.6 * length(experiments), height = 1.6 * 3)


# ===============================================================
# 3. Supplementary table: per-locus values and intervals, all experiments
# ===============================================================
stats_all <- c("M", "FST", "G'ST", "D")
ci_cols   <- c("value", "value_lower", "value_upper", "value_bias", "value_corrected")

new_names <- c(value           = "value",
               value_lower     = "CI_lower",
               value_upper     = "CI_upper",
               value_bias      = "bias",
               value_corrected = "corrected")

per_locus_wide <- function(e) {
  est <- res_ccf[[e]]$estimate
  ids <- prep[[e]]$data
  id_cols <- intersect(c("chr", "pos", "Ref", "Alt"), names(ids))

  # one block of columns per statistic, rows in locus order
  blocks <- lapply(stats_all, function(s) {
    f <- est[est$statistic == s, ]
    f <- f[order(f$locus), intersect(ci_cols, names(f))]
    names(f) <- paste0(s, "_", new_names[names(f)])       # e.g. value_lower -> CI_lower.M
    f
  })

  locus_ids <- est$locus[est$statistic == stats_all[1]]
  locus_ids <- sort(locus_ids)

  cbind(experiment = e,
        K = prep[[e]]$K,
        ids[locus_ids, id_cols, drop = FALSE],
        do.call(cbind, blocks))
}

tab_per_locus <- do.call(rbind, lapply(names(res_ccf), per_locus_wide))
rownames(tab_per_locus) <- NULL
head(tab_per_locus)
write.csv(tab_per_locus, file.path(dir_out, "organoids/tables/ccf_bootstrap_per_locus.csv"), row.names = FALSE)


# ===============================================================
# 4. Table: mean over loci, with interval, per experiment
# ===============================================================

sem <- function(x) { x <- x[!is.na(x)]; sd(x) / sqrt(length(x)) }

tab_mean <- do.call(rbind, lapply(experiments, function(e) {
  est <- res_ccf[[e]]$estimate
  m   <- res_ccf[[e]]$mean
  data.frame(experiment      = e,
             K               = prep[[e]]$K,
             n_loci          = sum(!is.na(est$value[est$statistic == "FST"])),
             statistic       = m$statistic,
             value           = m$value,
             CI_lower        = m$value_lower,
             CI_upper        = m$value_upper,
             CI_width        = m$value_upper - m$value_lower,
             value_bias      = m$value_bias,
             value_corrected = m$value_corrected,
             SEM_loci        = unname(sapply(m$statistic, function(s) sem(est$value[est$statistic == s]))))
}))
tab_mean

write.csv(tab_mean, file.path(dir_out, "organoids/tables/ccf_bootstrap_per_sample_mean.csv"), row.names = FALSE)

# ===============================================================
# 5. Summary of per-locus intervals, per experiment and statistic
# ===============================================================
summarise_loci <- function(e) {
  est <- res_ccf[[e]]$estimate
  do.call(rbind, lapply(c("M", stats_diff), function(s) {
    f <- est[est$statistic == s & !is.na(est$value), ]
    w <- f$value_upper - f$value_lower
    data.frame(experiment       = e,
               K                = prep[[e]]$K,
               statistic        = s,
               n_loci           = nrow(f),
               value_min        = min(f$value),
               value_median     = median(f$value),
               value_max        = max(f$value),
               ci_width_median  = median(w, na.rm = TRUE),
               ci_width_mean    = mean(w, na.rm = TRUE),
               ci_width_max     = max(w, na.rm = TRUE),
               bias_median      = median(f$value_bias, na.rm = TRUE))
  }))
}
tab_summary <- do.call(rbind, lapply(experiments, summarise_loci))
rownames(tab_summary) <- NULL
tab_summary

write.csv(tab_summary, file.path(dir_out, "organoids/tables/ccf_bootstrap_summary.csv"), row.names = FALSE)

# ===============================================================
# Description sheet and export
# ===============================================================
B     <- res_ccf[[1]]$B
level <- res_ccf[[1]]$level
type  <- res_ccf[[1]]$type

tab_readme <- data.frame(
  item = c("Method", "Replicates", "Interval",
           "per_locus", "mean_over_loci", "summary_per_experiment",
           "value", "value_lower / value_upper", "value_bias", "value_corrected",
           "share_value_outside_ci"),
  description = c(
    "Parametric bootstrap of sequencing read sampling: CCF* = Binomial(DP, CCF) / DP per locus and sample. Purity and copy number not modelled (lower bound on read-sampling uncertainty).",
    paste("B =", B, "replicates per experiment"),
    paste0(level * 100, "% ", type, " bootstrap interval"),
    "Sheet: statistics and intervals for each locus (subclonal mutation) and experiment",
    "Sheet: mean over loci of each statistic, with interval computed from the mean within each replicate",
    "Sheet: summary over loci of values and per-locus interval widths, per experiment",
    "Statistic computed on the original data",
    "Confidence limits",
    "Mean of bootstrap values minus value (estimated bias)",
    "value minus value_bias (bias-corrected estimate)",
    "Proportion of loci whose observed value lies outside its own interval"))

write_xlsx(list(README                 = tab_readme,
                per_locus              = tab_per_locus,
                mean_over_loci         = tab_mean,
                summary_per_experiment = tab_summary),
           file.path(dir_out, "Table_S_ccf_bootstrap.xlsx"))

# CSV copies, if needed
write.csv(tab_per_locus, file.path(out_dir, "ccf_bootstrap_per_locus.csv"),              row.names = FALSE)
write.csv(tab_mean,      file.path(out_dir, "ccf_bootstrap_mean_over_loci.csv"),         row.names = FALSE)
write.csv(tab_summary,   file.path(out_dir, "ccf_bootstrap_summary_per_experiment.csv"), row.names = FALSE)
