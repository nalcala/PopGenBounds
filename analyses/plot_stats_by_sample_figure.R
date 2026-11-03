library(PopGenBounds)
library(tibble)
library(patchwork)
library(dplyr)
library(ggplot2)

# Load data ----

load("C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/Organoids_CCF_clusters_dataset_fracs_list.RData")

data_3pop <- dataset_fracs_list[c("LCNEC4", "PANEC1")]
data_2pop <- dataset_fracs_list[c("SINET8", "SINET9", "LNET6", "LNET10", "LCNEC3")]

data_2pop <- unparse_columns(data_2pop)
data_3pop <- unparse_columns(data_3pop)

# Compute stats for all the samples ----

sample_names <- c("LNET6", "SINET8", "SINET9","LNET10", "LCNEC3", "LCNEC4", "PANEC1")

compute_diff_stats <- function(data_frame_list, sample_names, K=2){
  diff_stats <- list()
  for (df in data_frame_list) {
    # Extract the sample name from the 3rd column name
    sample_name <- sub(".*\\.", "", colnames(df)[3])
    sample_name <- substr(sample_name, 1, nchar(sample_name) - 1)

    # Skip the iteration if sample not in the provided list
    if (!(sample_name %in% sample_names)){
      next
    }
    print(sample_name)
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_clonal(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    Diff_loci <- lapply(list_freq, Diff)
    diff_stats[[sample_name]] <- Diff_loci
  }
  return(diff_stats)
}

# Plot samples in one figure ----

plot_diff_stats <- function(diff_stats, K, is_normalised=FALSE){
  plots_list <- list()

  # creating the first column
  dummy1 <- ggplot() +
    labs(y = expression(bolditalic(F[ST]))) +
    theme_void() +
    theme(axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5))

  dummy2 <- ggplot() +
    labs(y = expression(bolditalic(G*"'"[ST]))) +
    theme_void() +
    theme(axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5))

  dummy3 <- ggplot() +
    labs(y = expression(bolditalic(D))) +
    theme_void() +
    theme(axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5))

  dummy_void <- ggplot() +
    theme_void()

  for (sample_name in names(diff_stats)) {
    print(paste("Sample:", sample_name))
    print(paste("K:", K))

    if (sample_name == "LNET6"){
      labels_col <- dummy1 / dummy2 / dummy3 / dummy_void

      plots_list[["labels"]] <- labels_col

    }

    d_stats <- diff_stats[[sample_name]]
    lnen.tib <- turn_diff_in_tibble(d_stats)

    if (is_normalised == FALSE){
      plot_gg <- ggbounds_raw(M=lnen.tib$M,
                               FST=lnen.tib$FST,
                               GpST=lnen.tib$GpST,
                               D=lnen.tib$D,
                               K=K)
      plot_gg <- lapply(plot_gg, function(p){
        p + scale_x_continuous(labels = function(x) {
          sapply(x, function(val) {
            if (val %in% c(0, 1)) {
              as.character(val)
            } else {
              format(round(val, 2), nsmall = 2)
            }
          })
        })
      })
    }
    else if (is_normalised == TRUE){
      plot_gg <- ggbounds_norm(M=lnen.tib$M,
                              FST=lnen.tib$FST,
                              GpST=lnen.tib$GpST,
                              D=lnen.tib$D,
                              K=K)
      plot_gg <- lapply(plot_gg, function(p) {
        # Round x-axis labels to 1 decimal place
        p + scale_x_continuous(labels = function(x) round(x, 1))
        })

    }

    # Leave y-label only in the first column and format the layout

    plot_gg <- lapply(plot_gg, function(p) {
      if (sample_name != "LNET6T") {
        p <- p + theme(
          #axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
      else {
        p <- p +
          scale_y_continuous(labels = function(x) {
              sapply(x, function(val) {
                if (val %in% c(0, 1)) {
                  as.character(val)
                } else {
                  format(round(val, 2), nsmall = 2)
                }
              })
            })
      }

      p + theme(
        axis.title = element_blank(),
        plot.margin = margin(t = 5, r = 4, b = 5, l = 4),  # reduce left/right margins
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 9),
        )
    })

    dummy_M <- ggplot() +
      labs(x = expression(italic(M))) +
      theme_void() +
      theme(axis.title.x = element_text(size = 10, hjust = 0.5))

    # Add the title directly to the first plot
    plot_gg[[1]] <- plot_gg[[1]] + ggtitle(sample_name) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 5))
      )

    plot_sample <- plot_gg[[1]]/plot_gg[[2]]/plot_gg[[3]]/dummy_M +
      patchwork::plot_annotation(theme = theme(
                                   plot.margin = margin(1, 1, 1, 1)
                                 )
      ) +
      plot_layout(axes = "collect")

    plots_list[[sample_name]] <- plot_sample
  }

  return(plots_list)
}

make_figure_1 <- function(diff_stats, samples_2pop, samples_3pop, is_normalised) {
  # plot distribution of statistics for each sample in one figure (raw values)
  plot_list_2subpop <- plot_diff_stats(diff_stats[samples_2pop], K=2, is_normalised)
  plot_list_3subpop <- plot_diff_stats(diff_stats[samples_3pop], K=3, is_normalised)

  str(plot_list_2subpop)

  plot_list_merged <- c(plot_list_2subpop, plot_list_3subpop)

  sample_names <- c(samples_2pop, samples_3pop)
  print(sample_names)
  sample_names_plot <- c("labels", sample_names)

  ordered_plots <- plot_list_merged[sample_names_plot]
  print(ordered_plots)

  final_plot <- wrap_plots(ordered_plots, nrow = 1)

  return(final_plot)
}

# compute differentiation statistics for each sample
diff_stats_2subpop <- compute_diff_stats(data_2pop,sample_names,K=2)
diff_stats_3subpop <- compute_diff_stats(data_3pop,sample_names,K=3)
diff_stats <- c(diff_stats_2subpop, diff_stats_3subpop)

# plot distribution of statistics for each sample in one figure

#sample_names <- c("LNET6T", "SINET8M", "SINET9M", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")
samples_2pop <- c("LNET6", "SINET8", "SINET9", "LNET10", "LCNEC3")
samples_3pop <- c("LCNEC4", "PANEC1")

figure_1_raw <- make_figure_1(diff_stats, samples_2pop, samples_3pop, is_normalised = FALSE)
figure_1_raw

figure_1_norm <- make_figure_1(diff_stats, samples_2pop, samples_3pop, is_normalised = TRUE)
figure_1_norm

combined_figure <- figure_1_raw / figure_1_norm
combined_figure


# save plots

output_folder <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/plots/"

ggsave(paste0(output_folder, "stats_distr_raw.svg"),
       plot = figure_1_raw,
       height = 6, width = 14)

ggsave(paste0(output_folder, "stats_distr_norm.svg"),
       plot = figure_1_norm,
       height = 7, width = 15)

ggsave(paste0(output_folder, "stats_distr_combined.svg"),
       plot = combined_figure,
       height = 12, width = 13)



## Figure for CNV data ----

# load CNV data
load("C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/Organoids_CNV_tables.RData")

compute_diff_CNV <- function(CNV_tables){
  diff_results_all <- list()

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

    #plots[[sample_name]] <- plot_stats(diff_results, K=ncol(df)-3, title = sample_name)
  }

  return(diff_results_all)
}

diff_results_CNV <- compute_diff_CNV(CNV_tables)

# plot distribution of statistics for each sample in one figure

samples_2pop <- c("LNET6", "SINET8", "SINET9", "LNET10", "LCNEC3")
samples_3pop <- c("LCNEC4", "PANEC1")

plot_list_2subpop <- plot_diff_stats(diff_results_CNV[samples_2pop], K=2, is_normalised=FALSE)

sample_names_plot <- c("labels", samples_2pop)
ordered_plots <- plot_list_2subpop[sample_names_plot]

final_plot <- wrap_plots(ordered_plots, nrow = 1)

figure_1_raw <- make_figure_1(diff_results_CNV, samples_2pop, samples_3pop, is_normalised = FALSE)
figure_1_raw

figure_1_norm <- make_figure_1(diff_results_CNV, samples_2pop, samples_3pop, is_normalised = TRUE)
figure_1_norm

combined_figure <- figure_1_raw / figure_1_norm
combined_figure

ggsave(paste0(output_folder, "stats_distr_combined_CNV.svg"),
       plot = combined_figure,
       height = 12, width = 13)


library(ggpointdensity)

dat <- bind_rows(
  tibble(x = rnorm(7000, sd = 1),
         y = rnorm(7000, sd = 10),
         group = "foo"),
  tibble(x = rnorm(3000, mean = 1, sd = .5),
         y = rnorm(3000, mean = 7, sd = 5),
         group = "bar"))

ggplot(data = dat,
       aes( x = x, y = y, color = after_stat(ndensity))) +
  geom_pointdensity( size = .25) +
  scale_color_viridis() +
  facet_wrap( ~ group) +
  labs(color = "relative\ndensity")

