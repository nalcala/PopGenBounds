#' Transform a list of matrices containing differentiation statistics into a tibble.
#' Exclude rows containing nan values.
#'
#' @param Diff_loci A list of matrices containing M, FST, G'ST and D values
#' for each loci in a sample. The values expected to be in the second column.
#'
#' @returns A tibble object without nan values.
#' @export
#'
#' @examples
#' library(tibble)
#' # Create example list of 2 matrices with 4 rows (M, FST, G'ST, D)
#' Diff_loci <- list(
#'   matrix(c("M", 0.85, NA, NA,
#'            "FST", 0.01, NA, NA,
#'            "G'ST", 0.02, NA, NA,
#'            "D", 0.015, NA, NA),
#'          ncol = 4, byrow = TRUE,
#'          dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness"))),
#'
#'   matrix(c("M", 0.5, NA, NA,
#'            "FST", 0.02, NA, NA,
#'            "G'ST", 0.03, NA, NA,
#'            "D", 0.025, NA, NA),
#'          ncol = 4, byrow = TRUE,
#'          dimnames = list(NULL, c("statistic", "value", "bound", "bound_closeness")))
#' )
#'
#' # Convert to tibble and remove rows with NaNs
#' clean_tibble <- turn_diff_in_tibble(Diff_loci)
#' print(clean_tibble)
turn_diff_in_tibble <- function(Diff_loci){
  diff.tib = tibble(M=unlist(sapply(Diff_loci,function(x){x[1,2]})),
                     FST=unlist(sapply(Diff_loci,function(x){x[2,2]})),
                     GpST=unlist(sapply(Diff_loci,function(x){x[3,2]})),
                     D=unlist(sapply(Diff_loci,function(x){x[4,2]})))

  diff.tib_clean <- diff.tib[!apply(is.nan(as.matrix(diff.tib)), 1, any), ]

  return(diff.tib_clean)
}

#' Swap two columns in a dataframe without changing other columns
#'
#' @param df A dataframe with at least two columns.
#' @param col1 Name of the first columns to swap.
#' @param col2 Name of the second column to swap.
#'
#' @returns A datafrme with specified two columns swapped.
#' @export
#'
#' @examples
#' # Create a sample data frame
#' df <- data.frame(A = 1:3, B = 4:6, C = 7:9)
#' # Swap columns B and C
#' swapped_df <- swap_columns(df, "B", "C")
swap_columns <- function(df, col1, col2) {
  cols <- colnames(df)
  idx1 <- which(cols == col1)
  idx2 <- which(cols == col2)
  cols[c(idx1, idx2)] <- cols[c(idx2, idx1)]
  df[, cols]
}

#' Transform a combined column in multiple separate columns
#'
#' If a dataframe column consist of two columns, replaces the combined column by the two separate columns.
#'
#' @param data_list A list of dataframes
#' @param idx Index of the column to replace
#' @export
#'
#' @returns List of modified dataframes
unparse_columns <- function(data_list, idx=3){
  print("Subpopulations:")
  for (i in seq_along(data_list)){
    df <- data_list[[i]]

    subpop_names <- colnames(df[[idx]])
    print(subpop_names)

    # Insert each subclonal fraction column with a generalized name
    for (j in seq_along(subpop_names)) {
      name_j <- subpop_names[j]
      col_vector <- df$subclonal.fractions[, name_j]

      # Insert at position 2 + (j - 1) to keep the new columns adjacent
      df <- add_column(df, !!paste0("subclonal.fractions.", name_j) := col_vector, .after = 1 + j)
    }

    df$subclonal.fractions <- NULL

    # Reassign back to the list
    data_list[[i]] <- df
  }
  return(data_list)
}
