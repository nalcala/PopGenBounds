#' Produces a combined plot with raw and normalized differentiation statistics for a given sample.
#'
#' @param Diff_loci A list of matrices containing values of differentiation statistics \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D}
#' and frequency of the most frequent allele \eqn{M}.
#' @param K Number of subpopulations
#' @param title Title of the produced plot
#'
#' @returns A ggplot object consisting of six statistics distributions, the first column correspond
#' to the distributions of the original values of \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D}
#' and the second column to their normalized versions.
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
#' combined_plot <- plot_stats(Diff_loci, K=2, title=glue::glue("K=2: Subpop1, Subpop2"))
plot_stats <- function(Diff_loci, K=2, title="") {

  lcnec.tib_clean <- turn_diff_in_tibble(Diff_loci)

  gglcnec = ggbounds_new(M=lcnec.tib_clean$M,
                         FST=lcnec.tib_clean$FST,
                         GpST=lcnec.tib_clean$GpST,
                         D=lcnec.tib_clean$D,
                         K=K)

  combined_plot <-
    (gglcnec[[1]][[1]] + gglcnec[[2]][[1]]) /
    (gglcnec[[1]][[2]] + gglcnec[[2]][[2]]) /
    (gglcnec[[1]][[3]] + gglcnec[[2]][[3]]) +
    patchwork::plot_annotation(title = title,
                    theme = theme(
                      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 20))
                    )
    )

  return(combined_plot)
}

#' Plot differentiation statistics for a given sample in case of K = 3 subpopulations
#'
#' @param D_loci Four lists of matrices containing differentiation statistics FST, G'ST, D and a vector of frequencies M.
#' The first list corresponds to 3 subpopulations together, the remaining lists
#' correspond to pair-wise differentiation statistics for each pair of subpopulations out of 3 subpopulations.
#' @param subpop_names A string vector containing names of subpopulations to be plotted.
#'
#' @returns A stacked ggplot of four elements. The first plot corresponds to values of
#' differentiation statistics computed for \eqn{K=3} subpopulations, the next three plots
#' correspond to differentiation statistics computed across pairs of subpopulations.
#' @import glue
#' @export
#'
#' @examples
#' library(patchwork)
#' library(glue)
#' library(tibble)
#' # Initialise allele frequencies for 3 subpopulations in 2 loci
#' locus1 <- matrix(c(0.8, 0.2,   # Subpop 1
#'                    0.6, 0.4,   # Subpop 2
#'                    0.7, 0.3),  # Subpop 3
#'                  nrow = 3, byrow = TRUE)
#'
#' locus2 <- matrix(c(0.5, 0.5,
#'                    0.4, 0.6,
#'                    0.6, 0.4),
#'                  nrow = 3, byrow = TRUE)
#'
#' # Compute differentiation statistics
#' diff_stats <- compute_Diff_3subpop(list(locus1, locus2))
#'
#' # Produce plot of differentiation statistics values
#' final_plot <- plot_stats_3subpop(diff_stats, c("Subpop1", "Subpop2", "Subpop3"))
plot_stats_3subpop <- function(D_loci, subpop_names){
  plot123 <- plot_stats(D_loci$D_123, K=3, title=glue("K=3: {subpop_names[1]}, {subpop_names[2]}, {subpop_names[3]}"))
  plot12 <- plot_stats(D_loci$D_12, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
  plot23 <- plot_stats(D_loci$D_23, K=2, title=glue("K=2: {subpop_names[2]}, {subpop_names[3]}"))
  plot13 <- plot_stats(D_loci$D_13, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[3]}"))

  final_plot <- wrap_elements(plot123) /
    wrap_elements(plot12) /
    wrap_elements(plot23) /
    wrap_elements(plot13)

  return(final_plot)
}

#' Pipeline to compute distributions of differentiation statistics for LNEN data.
#'
#' @param data_frame_list A list of dataframes with LNEN samples.
#' @param K Number of subpopulations. K should be the same for all samples in data_frame_list.
#' @param output_folder An output directory in str format.
#'
#' @export
run_lnen_plotting <- function(data_frame_list, K=2, output_folder) {
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_clonal(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      combined_plot <- plot_stats(Diff_loci, K, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
    }
    else if (K==3){
      D_loci = compute_Diff_3subpop(list_freq)
      combined_plot <- plot_stats_3subpop(D_loci, subpop_names)
      combined_plot
    }
    save_plots(subpop_names, combined_plot, output_folder, type="stats", K)
  }
}


#' Make plots of distributions of differentiation statistics \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D}
#' and their normalized versions.
#'
#' Returns a list of plots of \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D} values versus
#' the mean frequency of the most frequent allele \eqn{M}.
#'
#' @param diff_stats Named list of differentiation statistics values for each sample.
#' @param K The number of subpopulations
#' @param is_normalised Boolean indicator whether to plot original statistics or their normalized versions.
#'
#' @returns A named list of stacked plots, each element correspond to a stacked plot of 3 subplots,
#' corresponding for \eqn{F_{ST}}, \eqn{G'_{ST}}, \eqn{D} values plotted versus \eqn{M}.
#' @export
plot_diff_over_samples <- function(diff_stats, K, is_normalised=FALSE){
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

#' Function to produce a joint plot of differentiation statistics over samples.
#'
#' @param diff_stats Named list of differentiation statistics values for each sample.
#' @param samples_2pop
#' @param samples_3pop
#' @param is_normalised Bool value indicating whether the statistics values were normalized by the corresponding upper bound
#'
#' @returns ggplot object
#' @export
#'
make_figure_sup <- function(diff_stats, samples_2pop, samples_3pop, is_normalised) {
  # plot distribution of statistics for each sample in one figure (raw values)
  plot_list_2subpop <- plot_diff_over_samples(diff_stats[samples_2pop], K=2, is_normalised)
  plot_list_3subpop <- plot_diff_over_samples(diff_stats[samples_3pop], K=3, is_normalised)

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
