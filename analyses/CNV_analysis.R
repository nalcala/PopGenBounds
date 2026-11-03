library(PopGenBounds)
library(tibble)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

## Load data ----
load("C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/Organoids_CNV_tables.RData")

diff_results_all <- list()
plots <- list()

for (sample_name in names(CNV_tables)) {

  df <- CNV_tables[[sample_name]]

  diff_results <- df %>%
    group_by(pos) %>%
    group_map(~ {
      #print(.x)
      data_freq <- .x[, 3:ncol(.x)]
      mat <- as.matrix(data_freq)
      diff_input <- t(mat)
      Diff(diff_input)
    })

  # store results for this sample
  diff_results_all[[sample_name]] <- diff_results

  plots[[sample_name]] <- plot_stats(diff_results, K=ncol(df)-3, title = sample_name)
}


## Plot 2d frequencies ----

plot_2d_freq_CNV <- function(CNV_tables){
  freq_plots <- list()
  for (sample_name in names(CNV_tables)) {
    df <- CNV_tables[[sample_name]]
    K <- ncol(df)-3
    subpop_names <- colnames(df)[4:ncol(df)]

    print(cat("Sample name: ", sample_name))
    print(cat("K: ", K))
    print(cat("Subpop name: ", subpop_names))

    x_col <- colnames(df)[4]
    y_col <- colnames(df)[5]

    g <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_point() +
      labs(x = subpop_names[[1]], y = subpop_names[[2]]) +
      theme_classic() +
      ggtitle(subpop_names[[1]])

    if (K==3){
      x_col <- colnames(df)[3]
      y_col <- colnames(df)[5]

      g1 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_point() +
        labs(x = subpop_names[[1]], y = subpop_names[[3]]) +
        theme_classic()

      x_col <- colnames(df)[4]
      y_col <- colnames(df)[5]

      g2 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
        geom_point() +
        labs(x = subpop_names[[2]], y = subpop_names[[3]]) +
        theme_classic()

      g <- g / g1 / g2 +
        plot_annotation(
          title=subpop_names[[1]],
          theme = theme(
            plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 20))
          )
        )
    }
    freq_plots[[sample_name]] <- g
  }

  return(freq_plots)
}

freq_plots_CNV <- plot_2d_freq_CNV(CNV_tables)
freq_plots_CNV


combined_plot <- wrap_plots(freq_plots_CNV, nrow = 1)  # all in one row
combined_plot

combined_2pop <- wrap_plots(freq_plots_CNV[1:5], nrow = 1)
print(combined_2pop)

combined_3pop <- wrap_plots(freq_plots_CNV[6:7], nrow = 1)
print(combined_3pop)

## Compare mean stats ----
compute_mean_stats_pop <- function(diff_loci_list, K) {
  result <- lapply(names(diff_loci_list), function(name) {
    print(name)
    diff_loci <- diff_loci_list[[name]]
    compute_mean_stats(diff_loci, K = K)
  })
  names(result) <- names(diff_loci_list)
  return(result)
}

samples_2pop <- c("LNET6", "SINET8", "SINET9", "LNET10", "LCNEC3")
samples_3pop <- c("LCNEC4", "PANEC1")

mean_stats_list_2pop <- compute_mean_stats_pop(diff_results_all[samples_2pop], 2)
mean_stats_list_3pop <- compute_mean_stats_pop(diff_results_all[samples_3pop], 3)

mean_stats_merged <- c(mean_stats_list_2pop, mean_stats_list_3pop)

sample_order <- c("LNET6", "SINET8", "SINET9", "LNET10", "LCNEC3", "LCNEC4", "PANEC1")

g_mean_stats <- plot_mean_stats(mean_stats_merged, sample_order)
g_mean_stats



SINET8_CNV <- CNV_tables$SINET8

diff_results <- SINET8_CNV %>%
  group_by(pos) %>%
  group_map(~ {
    print(.x)
    data_freq <- .x[, 3:4]
    mat <- as.matrix(data_freq)
    diff_input <- t(mat)
    Diff(diff_input)
  })

diff_results

groups_SI8 <- SINET8_CNV %>%
  group_by(pos) %>%
  group_split()

data_freq_SINET8 <- groups_SI8[[1]][4:5]
mat <- as.matrix(data_freq_SINET8)
diff_input <- t(mat)
diff_input
Diff(diff_input)
