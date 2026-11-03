library(PopGenBounds)

# load data ----

path_to_data <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/small_variants_CCFs"
setwd(path_to_data)

data_frame_names <- list.files(pattern = "*.tsv")       # Get all file names

file_names_3_subpop <- c("PANEC1_annotatedvariants_CCF_clonality.tsv", "LCNEC4_annotatedvariants_CCF_clonality.tsv")
file_names_2_subpop <- data_frame_names[!data_frame_names %in% file_names_3_subpop]

data_frame_list_2pop <- lapply(file_names_2_subpop, read.delim)  # Read all data frames
data_frame_list_3pop <- lapply(file_names_3_subpop, read.delim)  # Read all data frames

# clean the data ----

## exclude other samples ----
# filter out loci form other samples

clean_samples <- function(df, K){
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
  
  return(data_filtered)
}

## exclude duplicated rows ----

clean_duplicated <- function(df){
  # check whether same loci have the same values in allele frequencies, Cluster and Clonal columns
  df_subset <- df[, c("pos", names(df)[3], names(df)[4], "Cluster", "Clonal")]
  ## Find rows with duplicated `pos`
  dups <- df_subset %>%
    group_by(pos) %>%
    filter(n() > 1) %>%
    ungroup()
  
  ## Check for variation in the other columns
  # This returns only the problematic positions
  inconsistent_dups <- dups %>% 
    group_by(pos) %>% 
    filter(n_distinct(across(everything())) > 1) %>%
    arrange(pos)
  
  cat("Inconsistent rows: ", dim(inconsistent_dups)[1] > 0, "\n")
  
  cat("Dimensions before de-doubling: ", dim(df), "\n")
  pos_before <- names(table(df$pos))
  
  # keep only one of rows in identical locus
  data_deduped <- df %>%
    distinct(pos, .keep_all = TRUE)
  
  # check number of rows for each chr position (locus) after the deduplication
  #plot(as.vector(table(data_deduped$pos)))
  pos_after <- names(table(data_deduped$pos))
  cat("Dimensions after de-doubling: ", dim(data_deduped), "\n")
  
  # check whether after de-duplication we didn't lose any loci
  are_equal <- setequal(pos_after, pos_before)
  cat("No locus is lost: ", are_equal, "\n")
  
  return(data_deduped)
}

data_clean_sample <- clean_samples(data_frame_list_2pop[[1]], K=2)
data_clean_duplic <- clean_duplicated(data_clean_sample)

## change the order of rows in samples with K=3

# print the order of subpopulations in samples with K=3

PANEC1 <- data_frame_list_3pop[[1]]
LCNEC4 <- data_frame_list_3pop[[2]]

print(subpop_names <- sub(".*\\.", "", colnames(PANEC1)[3:5]))
print(subpop_names <- sub(".*\\.", "", colnames(LCNEC4)[3:5]))

PopGenBounds::swap_columns(PANEC1, col1=colnames(PANEC1)[4], col2=colnames(PANEC1)[5])
PopGenBounds::swap_columns(LCNEC4, col1=colnames(LCNEC4)[4], col2=colnames(LCNEC4)[5])







