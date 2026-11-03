#' Computes mean values and normalized mean valuesof FST, G'ST and D
#' statistics over all loci in a given sample.
#'
#' @param Diff_loci A list of matrices containing values of differentiation statistics FST, G'ST, D
#' and a vector of frequencies M.
#' @param K Number of subpopulations
#'
#' @return A named list with six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' }
#' @export
#'
#' @examples
#' # Example with 2 loci and 2 subpopulations
#' library(tibble)
#' library(patchwork)
#' library(ggplot2)
#' library(glue)
#' locus1 <- matrix(c(0.8, 0.2,    # Subpop 1
#'                    0.6, 0.4),   # Subpop 2
#'                  nrow = 2, byrow = TRUE)
#'
#' locus2 <- matrix(c(0.5, 0.5,
#'                    0.4, 0.6),
#'                  nrow = 2, byrow = TRUE)
#'
#' list_loci <- list(locus1, locus2)
#' Diff_loci <- lapply(list_loci, Diff)
#' mean_stats <- compute_mean_stats(Diff_loci, K=2)
compute_mean_stats <- function(Diff_loci, K){

  lcnec.tib_clean <- turn_diff_in_tibble(Diff_loci)

  FST_mean <- mean(lcnec.tib_clean$FST, na.rm = TRUE)
  FST_norm <- lcnec.tib_clean$FST / sapply(lcnec.tib_clean$M, function(m) Fup(K, m))
  FST_norm_mean <- mean(FST_norm, na.rm=T)

  GST_mean <- mean(lcnec.tib_clean$GpST, na.rm = TRUE)
  GST_norm <- lcnec.tib_clean$GpST / sapply(lcnec.tib_clean$M, function(m) Gpup(K, m))
  GST_norm_mean <- mean(GST_norm, na.rm=T)

  D_mean <- mean(lcnec.tib_clean$D, na.rm = TRUE)
  D_norm <- lcnec.tib_clean$D / sapply(lcnec.tib_clean$M, function(m) Dup(K, m))
  D_norm_mean <- mean(D_norm, na.rm=T)

  return(list(
    FST_mean = FST_mean,
    GST_mean = GST_mean,
    D_mean = D_mean,
    FST_norm_mean = FST_norm_mean,
    GST_norm_mean = GST_norm_mean,
    D_norm_mean = D_norm_mean
  ))
}

#' Computes mean and normalized mean values of FST, G'ST, D for a list of samples.
#'
#' @param data_frame_list A list of dataframes containing allele fraquencies for loci (rows)
#' in multiple subpopulations (columns)
#' @param K Number of subpopulations
#' @param sample_names List of strings containing names of samples to include in the computation.
#'
#' @returns A list of named lists containing six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' } for a provided list of samples
#'
#' @export
get_mean_stats <- function(data_frame_list, sample_names, K=2){
  mean_stats_list <- list()
  for (df in data_frame_list) {
    # Extract the sample name from the 3rd column name
    sample <- sub(".*\\.", "", colnames(df)[3])

    # Skip the iteration if sample not in the provided list
    if (!(sample %in% sample_names)){
      next
    }

    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_clonal(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])

    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      mean_stats <- compute_mean_stats(Diff_loci, K)
    }

    else if (K==3){
      D_loci = compute_Diff_3subpop(list_freq)
      mean_stats <- compute_mean_stats(D_loci$D_123, K)
    }

    mean_stats_list[[sample]] <- mean_stats
  }
  return(mean_stats_list)
}

#' Plots mean values of FST, G'ST and D and their normalised mean values for a given list of samples.
#'
#' @param mean_stats_list A list of named list containing six numeric values:
#' \describe{
#'   \item{FST_mean}{Mean raw FST value across loci}
#'   \item{GST_mean}{Mean raw G'ST value across loci}
#'   \item{D_mean}{Mean raw Jost's D value across loci}
#'   \item{FST_norm_mean}{Mean normalized FST value across loci}
#'   \item{GST_norm_mean}{Mean normalized G'ST value across loci}
#'   \item{D_norm_mean}{Mean normalized Jost's D value across loci}
#' }
#' @param sample_order Order for samples to be plotted
#' Each named list corresponds to a different data sample.
#'
#' @returns A ggplot object
#' @export
#'
plot_mean_stats_sep <- function(mean_stats_list, sample_order){
  mean_stats_df <- bind_rows(lapply(mean_stats_list, as_tibble), .id = "Sample")
  #sample_order <- c("LNET6T", "SINET8M", "SINET9M", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")
  mean_stats_df$Sample <- factor(mean_stats_df$Sample, levels = sample_order)

  # Step 1: Reshape to long format
  mean_stats_long <- mean_stats_df %>%
    select(Sample, FST_mean, GST_mean, D_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Statistic", values_to = "Value")

  mean_norm_stats_long <- mean_stats_df %>%
    select(Sample, FST_norm_mean, GST_norm_mean, D_norm_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Statistic", values_to = "Value")

  # Define the order of metrics as you want them in the legend
  metric_order <- c("FST_mean", "GST_mean", "D_mean")
  metric_norm_order <- c("FST_norm_mean", "GST_norm_mean", "D_norm_mean")

  # Set factor levels to control legend order
  mean_stats_long$Statistic <- factor(mean_stats_long$Statistic, levels = metric_order)
  mean_norm_stats_long$Statistic <- factor(mean_norm_stats_long$Statistic, levels = metric_norm_order)

  # Step 2: Plot
  g1 <- ggplot(mean_stats_long, aes(x = Sample, y = Value, color = Statistic, group = Statistic)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "",
         x = "", y = "Mean Value") +
    theme(
      plot.title = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold"),  # Make legend title bold
      legend.text = element_text(size = 11)
      ) +
    scale_color_manual(
      values = c(
        "FST_mean" = "#C93842",
        "GST_mean" = "#2678B3",
        "D_mean" = "#65A856"
      ),
      labels = c(
        expression(italic(F[ST])),
        expression(italic(G*"'"[ST])),
        expression(italic(D))
      )
    ) +
    ggtitle("(a)")

  g2 <- ggplot(mean_norm_stats_long, aes(x = Sample, y = Value, color = Statistic, group = Statistic)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "",
         x = "", y = "Mean Value") +
    theme(
      plot.title = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold"),  # Make legend title bold
      legend.text = element_text(size = 11)
    ) +
    scale_color_manual(
      values = c(
        "FST_norm_mean" = "#C93842",
        "GST_norm_mean" = "#2678B3",
        "D_norm_mean" = "#65A856"
      ),
      labels = c(
        expression(italic(bar(F)[ST])),
        expression(italic(bar(G)*"'"[ST])),
        expression(italic(bar(D)))
      )
    ) +
  ggtitle("(b)")


  combined_plot <- (g1 + g2)
    #plot_annotation(
    #  tag_levels = 'a',           # Use lowercase letters
    #  tag_prefix = '(',           # Add opening bracket
    #  tag_suffix = ')',           # Add closing bracket
    #  theme = theme(
    #    plot.tag = element_text(size = 55),  # Increase font size
    #    plot.tag.position = c(0, 1)  # x=0 (left), y=1 (top), closer to plot
    #  )
    #)
  return(combined_plot)
}


#' Plot mean values of original and normalized differentiation statistics in provided samples.
#'
#' This function creates a single plot showing the mean values of \eqn{F_{ST}},
#' \eqn{G'_{ST}}, and Jost's \eqn{D}, along with their normalized counterparts,
#' for a given list of samples.
#'
#' @param mean_stats_list A list of named lists, each corresponding to a data sample,
#'   containing six numeric values:
#'   \describe{
#'     \item{FST_mean}{Mean raw \eqn{F_{ST}} value across loci}
#'     \item{GST_mean}{Mean raw \eqn{G'_{ST}} value across loci}
#'     \item{D_mean}{Mean raw Jost's \eqn{D} value across loci}
#'     \item{FST_norm_mean}{Mean normalized \eqn{F_{ST}} value across loci}
#'     \item{GST_norm_mean}{Mean normalized \eqn{G'_{ST}} value across loci}
#'     \item{D_norm_mean}{Mean normalized Jost's \eqn{D} value across loci}
#'   }
#' @param sample_order A vector specifying the order in which samples should be plotted.
#'
#' @return A \code{ggplot} object representing the comparison of
#' mean values of the differentiation statistics.
#' @export
plot_mean_stats <- function(mean_stats_list, sample_order){
  mean_stats_df <- bind_rows(lapply(mean_stats_list, as_tibble), .id = "Sample")
  #sample_order <- c("LNET6T", "SINET8M", "SINET9M", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")
  mean_stats_df$Sample <- factor(mean_stats_df$Sample, levels = sample_order)

  # Step 1: Reshape to long format
  mean_stats_long_combined <- mean_stats_df %>%
    select(Sample, FST_mean, GST_mean, D_mean, FST_norm_mean, GST_norm_mean, D_norm_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Statistic", values_to = "Value") %>%
    mutate(
      Metric = case_when(
        grepl("^FST", Statistic) ~ "FST",
        grepl("^GST", Statistic) ~ "GST",
        grepl("^D", Statistic)   ~ "D"
      ),
      Metric = factor(Metric, levels = c("FST", "GST", "D")),
      Normalized = grepl("norm", Statistic)
    )
  print(head(mean_stats_long_combined))
  # Define the order of metrics as you want them in the legend
  #metric_order <- c("FST_mean", "GST_mean", "D_mean")
  #metric_norm_order <- c("FST_norm_mean", "GST_norm_mean", "D_norm_mean")

  # Set factor levels to control legend order
  #mean_stats_long$Statistic <- factor(mean_stats_long$Statistic, levels = metric_order)
  #mean_norm_stats_long$Statistic <- factor(mean_norm_stats_long$Statistic, levels = metric_norm_order)

  # Step 2: Plot
  g1 <- ggplot(mean_stats_long_combined, aes(x = Sample,
                                             y = Value,
                                             color = Metric,
                                             linetype = Normalized,
                                             group = interaction(Metric, Normalized))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = "",
         x = "",
         y = "Mean Value",
         linetype = "Normalized",
         color = "Statistic") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=12, color="black"),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(face = "bold"),  # Make legend title bold
      legend.text = element_text(size = 11)
    ) +
    scale_color_manual(
      values = c(
        "FST" = "#C93842",
        "GST" = "#2678B3",
        "D" = "#65A856"
      ),
      labels = c(
        expression(italic(F[ST])),
        expression(italic(G*"'"[ST])),
        expression(italic(D))
      )
    ) +
    scale_linetype_manual(values = c('YES' = "solid", 'NO' = "dotted"))


  return(g1)
}


