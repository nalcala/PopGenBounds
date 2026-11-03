#' Plot a 2d-scatter plot of allele frequencies for each pair of subpopulations in a sample
#'
#' @param df A dataframe with loci in rows and allele freqencies in subpopulations in columns
#' @param K Number of subpopulations
#'
#' @returns A ggplot object
#' @export
#'
plot_2d_freq <- function(df, K){

  subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
  print(cat("Subpop names: ", subpop_names))

  x_col <- colnames(df)[3]
  y_col <- colnames(df)[4]

  g <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color=.data[["Cluster"]])) +
    geom_point() +
    labs(x = subpop_names[[1]], y = subpop_names[[2]]) +
    theme_classic() +
    ggtitle(subpop_names[[1]])

  if (K==3){
    x_col <- colnames(df)[3]
    y_col <- colnames(df)[5]

    g1 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color=.data[["Cluster"]])) +
      geom_point() +
      labs(x = subpop_names[[1]], y = subpop_names[[3]]) +
      theme_classic()

    x_col <- colnames(df)[4]
    y_col <- colnames(df)[5]

    g2 <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color=.data[["Cluster"]])) +
      geom_point() +
      labs(x = subpop_names[[2]], y = subpop_names[[3]]) +
      theme_classic()

    g <- g + g1 + g2 +
      plot_annotation(
        title=subpop_names[[1]],
        theme = theme(
          plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 20))
        )
      )
  }

  return(g)
}

#' Run a pipeline for plotting pair-wise allele frequencies of subpopulations for a list of samples
#'
#' @param data_frame_list A list of dataframes containing allele fraquencies for loci (rows)
#' in multiple subpopulations (columns)
#' @param K Number of subpopulations
#' @param output_folder An output directory in str format.
#'
#' @export
run_2d_freq_plotting <- function(data_frame_list, K=2, output_folder){
  for (df in data_frame_list) {
    data_clean <- filter_clonal(df, K, filter_fixed=F)
    freq_plot <- plot_2d_freq(data_clean, K)
    save_plots(df, freq_plot, output_folder, type="freq", K)
  }
}
