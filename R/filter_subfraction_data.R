#' Filter subclonal fractions dataframe to exclude clonal mutations.
#'
#' Keep only subclonal mutations ("Clonal" = FALSE) and exclude loci with subfraction frequency = 0.
#'
#' @param df A dataframe with loci in rows and different metrics in columns.
#' @param K Number of subpopulations, K=2 by default.
#' @param filter_fixed bool value, defines whether to filter out loci with 0 frequencies in both subpopulations.
#'
#' @returns A dataframe containing only rows with "Clonal" = FALSE.
#' @export
#'
#' @examples
#' #' library(dplyr)
#' # Simulated example dataset
#' df <- tibble::tibble(
#'   chr = c("1", "1", "1"),
#'   pos = c(12345, 12346, 12347),
#'   sub.frac.Sample1 = c(0.2, 0, 0.5),
#'   sub.frac.Sample2 = c(0.3, 0, 0.5),
#'   Clonal = c("FALSE", "TRUE", "Uncertain")
#' )
#' filtered_data <- filter_clonal(df, K=2)
filter_clonal <- function(df, K=2, filter_fixed=T){
  cat("Initial data dimension: ", dim(df), "\n")
  # Include only subclonal loci
  data_subclonal <- subset(df, Clonal == "FALSE")
  cat("After filtering clonal loci: ", dim(data_subclonal), "\n")

  # Filter loci with 0 frequency in both subclones
  if (filter_fixed == T) {
    if (K == 2) {
      data_subclonal <- data_subclonal[!(data_subclonal[[3]] == 0 & data_subclonal[[4]] == 0),]
    }
    else if (K == 3) {
      data_subclonal <- data_subclonal[
        !(data_subclonal[[3]] == 0 & data_subclonal[[4]] == 0 & data_subclonal[[5]] == 0),]
    }

    cat("After filtering loci with no mutations: ", dim(data_subclonal), "\n")

  # Count the number of loci with fixed alleles
  if (K == 2) {
    cat("Loci having at least one fixed subpop", sum(data_subclonal[[3]] == 0 |
                                                       data_subclonal[[4]] == 0),
        "\n\n")
  }

  else if (K == 3) {
    cat("Loci having at least one fixed subpop", sum(data_subclonal[[3]] == 0 |
                                                       data_subclonal[[4]] == 0 |
                                                       data_subclonal[[5]] == 0), "\n\n")
    }
  }

  return(data_subclonal)
}



#' Filters out loci from other samples, duplicated rows, Clonal loci and loci with 0 in both subpopulations
#'
#' @param df A dataframe with loci in rows and different metrics in columns
#' @param K Number of subpopulations, K=2 by default
#' @param filter_fixed bool value, defines whether to filter out loci with 0 frequencies in both subpopulations
#'
#' @returns A cleaned dataframe for further pipeline
#' @export
#'
#' @examples
#' library(dplyr)
#' # Simulated example dataset
#' df <- tibble::tibble(
#'   chr = c("1", "1", "1"),
#'   pos = c(12345, 12346, 12347),
#'   sub.frac.Sample1 = c(0.2, 0, 0.5),
#'   sub.frac.Sample2 = c(0.3, 0, 0.5),
#'   Clonal = c(FALSE, TRUE, FALSE),
#'   Sample = c("Sample1", "Sample1", "Sample2")
#' )
#' filtered_data <- filter_data_old(df, K=2)
filter_data_old <- function(df, K=2, filter_fixed=T) {
  # filter out loci form other samples
  if (K == 2){
    subpop_names <- sub(".*\\.", "", colnames(df)[3:4])
  }

  else if (K == 3) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:5])
  }
  cat("Subpopulation: ", subpop_names, "\n")
  cat("Initial data dimension: ", dim(df), "\n")

  data_filtered <- df[df$Sample %in% subpop_names,]
  cat("After sample name filtering: ", dim(data_filtered), "\n")

  # Clean duplicated rows
  # number of rows for each chr position (locus)
  #plot(as.vector(table(data_filtered$pos)))
  pos_before <- names(table(data_filtered$pos))

  # keep only one of rows in identical locus
  data_deduped <- data_filtered %>%
    distinct(pos, .keep_all = TRUE)

  # check number of rows for each chr position (locus) after the deduplication
  #plot(as.vector(table(data_deduped$pos)))
  pos_after <- names(table(data_deduped$pos))
  cat("After de-doubling: ", dim(data_deduped), "\n")

  # check whether after de-duplication we didn't lose any loci
  are_equal <- setequal(pos_after, pos_before)
  cat("No locus is lost: ", are_equal, "\n")

  ## Filter out Clonal loci
  data_filtered <- subset(data_deduped, Clonal == "FALSE")
  cat("After filtering clonal loci: ", dim(data_filtered), "\n")

  if (filter_fixed == T) {
    ## Filter loci with 0 frequency in both subclones
    if (K == 2) {
      data_filtered <- data_filtered[!(data_filtered[[3]] == 0 & data_filtered[[4]] == 0),]
    }
    else if (K == 3) {
      data_filtered <- data_filtered[
        !(data_filtered[[3]] == 0 & data_filtered[[4]] == 0 & data_filtered[[5]] == 0),]
    }

    cat("After filtering loci with no mutations", dim(data_filtered), "\n")

    ## Count the number of loci with fixed alleles
    if (K == 2) {
      cat("Loci having at least one fixed subpop", sum(data_filtered[[3]] == 0 |
                                                         data_filtered[[4]] == 0),
          "\n\n")
    }

    else if (K == 3) {
      cat("Loci having at least one fixed subpop", sum(data_filtered[[3]] == 0 |
                                                         data_filtered[[4]] == 0 |
                                                         data_filtered[[5]] == 0), "\n\n")
    }
  }

  return(data_filtered)
}

