#' Split combined (matrix or data frame) columns into separate columns
#'
#' Replaces each combined column by its sub-columns, named "<column>.<subcolumn>",
#' and keeps them at the same position in the dataframe.
#'
#' @param data_list A list of dataframes
#' @param cols Names of the combined columns to split. If NULL, all combined columns are split.
#' @export
#'
#' @returns List of modified dataframes
unparse_columns <- function(data_list, cols = NULL) {
  lapply(data_list, function(df) {
    # Find combined columns automatically if none are given
    if (is.null(cols)) {
      cols <- names(df)[vapply(df, function(x) !is.null(dim(x)), logical(1))]
    }

    for (col in cols) {
      m <- df[[col]]
      if (is.null(dim(m))) next   # skip columns that are missing or already simple

      sub_names <- colnames(m)
      if (is.null(sub_names)) sub_names <- seq_len(ncol(m))
      message(col, ": ", paste(sub_names, collapse = ", "))

      new <- as.data.frame(m, check.names = FALSE)
      names(new) <- paste0(col, ".", sub_names)

      pos <- match(col, names(df))
      df <- dplyr::bind_cols(df[seq_len(pos - 1)], new, df[-seq_len(pos)])
    }
    df
  })
}

d_unparsed <- do.call(data.frame, c(as.list(d), check.names = FALSE))

# manuscript data
dir_data <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data"
organoids_ccf_data <- get(load(file.path(dir_data, "Organoids_CCF_clusters_dataset_fracs_list.RData")))

data_3pop <- organoids_ccf_data[c("LCNEC4", "PANEC1")]
data_2pop <- organoids_ccf_data[c("SINET8", "SINET9", "LNET6", "LNET10", "LCNEC3")]

LCNEC4_used <- data_3pop[[1]]
head(LCNEC4_used)
dim(LCNEC4_used)

LCNEC4_used_unparsed <- unparse_columns(data_3pop)[[1]]
head(LCNEC4_used_unparsed)
dim(LCNEC4_used_unparsed)
class(LCNEC4_used_unparsed)

LCNEC4_used_base_unparsed <- do.call(data.frame, c(as.list(LCNEC4_used), check.names = FALSE))
head(LCNEC4_used_base_unparsed)
dim(LCNEC4_used_base_unparsed)
class(LCNEC4_used_base_unparsed)

all.equal(LCNEC4_used_unparsed, LCNEC4_used_base_unparsed)


