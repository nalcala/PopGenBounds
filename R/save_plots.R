#' Save to file plots of differentiation statistics or of allele frequency.
#'
#' @param combined_plot A ggplot object to plot.
#' @param type A type of the plot: differentiation statistics (type = "stats"),
#' or allele frequencies (type = "freqs")
#' @param K Number of subpopulations.
#' @param subpop_names A string vector of subpopulations names to plot.
#' @param output_folder An output directory in str format.
#'
#' @returns None
#' @export
#'
save_plots <- function(subpop_names, combined_plot, output_folder, type="stats", K=2) {
  #output_folder <- "../../plots/"
  if (K == 2) {
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]

    filename <- paste0(output_folder, type, "_", pop1, "_vs_", pop2, ".svg")
    print(filename)
    if (type == "stats"){
      ggsave(filename,
             plot = combined_plot,
             height = 9, width = 8)
    }
    else if (type == "freqs"){
      ggsave(filename,
             plot = combined_plot,
             height = 3, width = 4)
    }
  }
  else if (K == 3) {
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]
    pop3 <- subpop_names[3]

    filename <- paste0(output_folder, type, "_", pop1, "_vs_", pop2, "_vs_", pop3, ".svg")
    print(filename)
    if (type=="stats"){
      ggsave(filename,
             plot = combined_plot,
             height = 36, width = 8)
    }
    else if (type == "freqs"){
      ggsave(filename,
             plot = combined_plot,
             height = 3, width = 12)
    }
  }
}
