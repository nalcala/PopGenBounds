library(PopGenBounds)
library(tibble)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)



# Load data ----

load("C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/Organoids_CCF_clusters_dataset_fracs_list.RData")

data_3pop <- dataset_fracs_list[c("LCNEC4", "PANEC1")]
data_2pop <- dataset_fracs_list[c("SINET8", "SINET9", "LNET6", "LNET10", "LCNEC3")]

#data_list <- c(data_2pop, data_3pop)

unparse_columns <- function(data_list){
  for (i in seq_along(data_list)){
    df <- data_list[[i]]

    subpop_names <- colnames(df[[3]])
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

data_2pop <- unparse_columns(data_2pop)
data_3pop <- unparse_columns(data_3pop)

output_folder <- output_folder <- "../plots/new/"

run_lnen_plotting(data_2pop, K=2, output_folder)
run_lnen_plotting(data_3pop, K=3, output_folder)

run_2d_freq_plotting(data_2pop, K=2, output_folder)
run_2d_freq_plotting(data_3pop, K=3, output_folder)

# Plot figure with mean values comparison over samples and normalization -----
sample_names <- c("LNET6T", "SINET8M", "SINET9M", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")

# compute mean values of differentiation statistics per sample
mean_stats_list_2pop <- get_mean_stats(data_2pop, sample_names, K=2)
mean_stats_list_3pop <- get_mean_stats(data_3pop, sample_names, K=3)

# merge samples containing 2 and 3 subpopulations
mean_stats_merged <- c(mean_stats_list_2pop, mean_stats_list_3pop)

g_mean_stats <- plot_mean_stats(mean_stats_merged)
g_mean_stats

SINET8 <- data_2pop[[1]]

subclonal.fractions.SINET8M <- SINET8$subclonal.fractions[, "SINET8M"]
subclonal.fractions.SINET8Mp2 <- SINET8$subclonal.fractions[, "SINET8Mp2"]

SINET8 <- add_column(SINET8, subclonal.fractions.SINET8M, .after = 2)
SINET8 <- add_column(SINET8, subclonal.fractions.SINET8Mp2, .after = 3)
SINET8$subclonal.fractions <- NULL


#SINET8$subclonal.fractions.SINET8M <- SINET8$subclonal.fractions[, "SINET8M"]
#SINET8$subclonal.fractions.SINET8Mp2 <- SINET8$subclonal.fractions[, "SINET8Mp2"]

data_subclonal <- subset(SINET8, Clonal == "FALSE")
data_subclonal <- data_subclonal[!(data_subclonal[[3]] == 0 & data_subclonal[[4]] == 0),]

SINET8_subclonal <- filter_clonal(SINET8, K=2)

SINET8_freq_plot <- plot_2d_freq(SINET8_subclonal, K=2)
SINET8_freq_plot

SINET8_subpop_names <- sub(".*\\.", "", colnames(SINET8)[3:4])
SINET8_subpop_names

SINET8_list_freq <- make_popgen_input(SINET8_subclonal[,3:4])
SINET8_list_freq[[1]]

SINET8_Diff_loci <- lapply(SINET8_list_freq, Diff)
SINET8_Diff_loci[[1]]

SINET8_combined_plot <- plot_stats(SINET8_Diff_loci, K=2,
                                    title=glue("K=2: {SINET8_subpop_names[1]},{SINET8_subpop_names[2]}"))
SINET8_combined_plot
