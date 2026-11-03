# Apply PopGenBounds functions to LNEN SNV samples

library(PopGenBounds)
library(tibble)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glue)


# Load data --------------------------------------

path_to_data <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/MS_panNEN_organoids/data/small_variants_CCFs"
setwd(path_to_data)

data_frame_names <- list.files(pattern = "*.tsv")       # Get all file names
print(data_frame_names)

file_names_3_subpop <- c("PANEC1_annotatedvariants_CCF_clonality.tsv", "LCNEC4_annotatedvariants_CCF_clonality.tsv")
file_names_2_subpop <- data_frame_names[!data_frame_names %in% file_names_3_subpop]
print(file_names_2_subpop)

data_frame_list_2pop <- lapply(file_names_2_subpop, read.delim)  # Read all data frames
data_frame_list_3pop <- lapply(file_names_3_subpop, read.delim)  # Read all data frames

print(file_names_3_subpop)

# Functions to process data ------------------------------------------------

filter_data <- function(df, K=2, filter_fixed=T) {
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

make_popgen_input <- function(data_frequencies) {

  # make list of allele frequency matrices for each locus
  list_loci <- lapply(1:nrow(data_frequencies), function(i){
    mat <- as.matrix(data_frequencies[i,])
    second_allele <- 1 - mat[1, ]  # Compute 1-x for each element in the first row
    new_mat <- rbind(mat, second_allele)
    rownames(new_mat) <- c("First", "Second")
    return(t(new_mat))
  }
  )
  return(list_loci)
}

compute_Diff_3subpop <- function(list_loci){
  # check whether K=3 (to do)

  # Extract couple of rows (subpops) from 3-subpop matrices
  list_loci_12 <- lapply(list_loci, function(mat) mat[c(1, 2), ])
  list_loci_23 <- lapply(list_loci, function(mat) mat[c(2, 3), ])
  list_loci_13 <- lapply(list_loci, function(mat) mat[c(1, 3), ])

  # Compute stats (F_ST, G_ST, D)
  D_loci_123 = lapply(list_loci, Diff)
  D_loci_12 = lapply(list_loci_12, Diff)
  D_loci_23 = lapply(list_loci_23, Diff)
  D_loci_13 = lapply(list_loci_13, Diff)

  return(list(D_123=D_loci_123, D_12=D_loci_12, D_13=D_loci_13, D_23=D_loci_23))

}

plot_stats <- function(Diff_loci, K=2, title="") {
  lcnec.tib = tibble(M=unlist( sapply(Diff_loci,function(x){x[1,2]})),
                     FST=unlist( sapply(Diff_loci,function(x){x[2,2]})),
                     GpST=unlist( sapply(Diff_loci,function(x){x[3,2]})),
                     D=unlist( sapply(Diff_loci,function(x){x[4,2]})))

  lcnec.tib[apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]

  lcnec.tib_clean <- lcnec.tib[!apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]

  gglcnec = ggbounds_new(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=K)

  #combined_plot <- gglcnec[[1]] + gglcnec[[2]]+ gglcnec[[3]] +
  #  plot_annotation(title = title,
  #                  theme = theme(
  #                  plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 20))
  #                  )
  #  )
  combined_plot <-
    (gglcnec[[1]][[1]] + gglcnec[[2]][[1]]) /
    (gglcnec[[1]][[2]] + gglcnec[[2]][[2]]) /
    (gglcnec[[1]][[3]] + gglcnec[[2]][[3]]) +
    plot_annotation(title = title,
                    theme = theme(
                      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 20))
                      )
                    )

  return(combined_plot)
}

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

save_plots <- function(df, combined_plot, type="stats", K=2) {
  if (K == 2) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:4])
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]

    output_folder <- "../../../plots/"
    filename <- paste0(output_folder, type, "_", pop1, "_vs_", pop2, ".svg")
    print(filename)
    if (type=="stats"){
      ggsave(filename,
             plot = combined_plot,
             height = 9, width = 8)
    }
    else{
      ggsave(filename,
             plot = combined_plot,
             height = 3, width = 4)
    }

  }
  else if (K == 3) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:5])
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]
    pop3 <- subpop_names[3]

    output_folder <- "../../../plots/"
    filename <- paste0(output_folder, type, "_", pop1, "_vs_", pop2, "_vs_", pop3, ".svg")
    print(filename)
    if (type=="stats"){
      print("Starting saving the file")
      ggsave(filename,
             plot = combined_plot,
             height = 36, width = 8)
    }
    else{
      ggsave(filename,
             plot = combined_plot,
             height = 3, width = 12)
    }
  }
}

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

compute_mean_stats <- function(Diff_loci, K){
  # make matrix with diff stats values in columns
  lcnec.tib = tibble(M=unlist( sapply(Diff_loci,function(x){x[1,2]})),
                     FST=unlist( sapply(Diff_loci,function(x){x[2,2]})),
                     GpST=unlist( sapply(Diff_loci,function(x){x[3,2]})),
                     D=unlist( sapply(Diff_loci,function(x){x[4,2]})))

  lcnec.tib[apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]

  lcnec.tib_clean <- lcnec.tib[!apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]

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

get_mean_stats <- function(data_frame_list, K=2){
  mean_stats_list <- list()
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    sample <- subpop_names[[1]]
    data_clean <- filter_data(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      mean_stats <- compute_mean_stats(Diff_loci, K)
    }
    else if (K==3){
      print("Computing Diff for all pairs of subpopulations")
      D_loci = compute_Diff_3subpop(list_freq)
      mean_stats <- compute_mean_stats(D_loci$D_123, K)
    }
    mean_stats_list[[sample]] <- mean_stats
  }
  return(mean_stats_list)
}

plot_mean_stats <- function(mean_stats_list){
  mean_stats_df <- bind_rows(lapply(mean_stats_list, as_tibble), .id = "Sample")
  sample_order <- c("SINET8M", "SINET9M", "LNET6T", "LNET10T", "LCNEC3T", "LCNEC4T", "PANEC1T")
  mean_stats_df$Sample <- factor(mean_stats_df$Sample, levels = sample_order)

  # Step 1: Reshape to long format
  mean_stats_long <- mean_stats_df %>%
    select(Sample, FST_mean, GST_mean, D_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Metric", values_to = "Value")

  mean_norm_stats_long <- mean_stats_df %>%
    select(Sample, FST_norm_mean, GST_norm_mean, D_norm_mean) %>%
    pivot_longer(cols = -Sample, names_to = "Metric", values_to = "Value")

  # Define the order of metrics as you want them in the legend
  metric_order <- c("FST_mean", "GST_mean", "D_mean")
  metric_norm_order <- c("FST_norm_mean", "GST_norm_mean", "D_norm_mean")

  # Set factor levels to control legend order
  mean_stats_long$Metric <- factor(mean_stats_long$Metric, levels = metric_order)
  mean_norm_stats_long$Metric <- factor(mean_norm_stats_long$Metric, levels = metric_norm_order)

  # Step 2: Plot
  g1 <- ggplot(mean_stats_long, aes(x = Sample, y = Value, color = Metric, group = Metric)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "Mean Diff Statistics per Sample",
         x = "Sample", y = "Mean Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  g2 <- ggplot(mean_norm_stats_long, aes(x = Sample, y = Value, color = Metric, group = Metric)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    theme_bw() +
    labs(title = "Mean Normalised Diff Statistics per Sample",
         x = "Sample", y = "Mean Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  combined_plot <- g1 + g2
  return(combined_plot)
}

run_lnen_plotting <- function(data_frame_list, K=2) {
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_data(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    print(list_freq[[1]])
    if (K==2){
      Diff_loci <- lapply(list_freq, Diff)
      combined_plot <- plot_stats(Diff_loci, K, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
    }
    else if (K==3){
      print("Computing Diff for all pairs of subpopulations")
      D_loci = compute_Diff_3subpop(list_freq)
      print(D_loci$D_123[1])
      print("Plotting stats")
      combined_plot <- plot_stats_3subpop(D_loci, subpop_names)
      combined_plot
      print("Combined plot is successfully computed")
    }
    print("Saving plot")
    save_plots(subpop_names, combined_plot, type="stats", K)
  }
}

run_2d_freq_plotting <- function(data_frame_list, K=2){
  for (df in data_frame_list) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
    data_clean <- filter_data(df, K, filter_fixed=F)
    freq_plot <- plot_2d_freq(data_clean, K)
    save_plots(subpop_names, freq_plot, type="freq", K)
  }
}

# plot mean values of differentiation statistics per sample
mean_stats_list_2pop <- get_mean_stats(data_frame_list_2pop, K=2)
mean_stats_list_3pop <- get_mean_stats(data_frame_list_3pop, K=3)
mean_stats_merged <- c(mean_stats_list_2pop, mean_stats_list_3pop)

g_mean_stats <- plot_mean_stats(mean_stats_merged)

output_folder <- "../../../plots/"
filename <- paste0(output_folder, "mean stats.svg")
ggsave(filename,
       plot = g_mean_stats,
       height = 4, width = 12)

# plot distributions of diff statistics across loci for each sample
run_lcnec_plotting(data_frame_list_2pop, K=2)

run_lcnec_plotting(data_frame_list_3pop, K=3)

run_2d_freq_plotting(data_frame_list_2pop, K=2)

run_2d_freq_plotting(data_frame_list_3pop, K=3)


# Example on one dataset ------------------------------------

LCNEC3 <- data_frame_list_2pop[[1]] ##LCNEC3
K <- 2
dim(LCNEC3)

LNET6 <- data_frame_list_2pop[[5]] ##LNET6
View(LNET6)

PANEC1 <- data_frame_list_3pop[[1]]
K <- 3

## Filter dataset --------------------------------------------------------

df <- PANEC1

# keep only the rows corresponding to the current sample
subpop_names <- sub(".*\\.", "", colnames(df)[3:(2+K)])
cat("Subpopulation: ", subpop_names, "\n")

cat("Initial data dimension: ", dim(df), "\n")

data_filtered <- df[df$Sample %in% subpop_names,]
cat("After sample name filtering: ", dim(data_filtered), "\n")

# Clean duplicated rows
# number of rows for each chr position (locus)
plot(as.vector(table(data_filtered$pos)))
pos_before <- names(table(data_filtered$pos))

# keep only one of rows in identical locus
data_deduped <- data_filtered %>%
  distinct(pos, .keep_all = TRUE)

# check number of rows for each chr position (locus) after the deduplication
plot(as.vector(table(data_deduped$pos)))
pos_after <- names(table(data_deduped$pos))
cat("After de-doubling: ", dim(data_deduped), "\n")

# check whether after de-duplication we didn't lose any loci
are_equal <- setequal(pos_after, pos_before)
cat("No locus is lost: ", are_equal, "\n")

## Filter out Clonal loci
data_filtered <- subset(data_deduped, Clonal == "FALSE")
cat("After filtering clonal loci: ", dim(data_filtered), "\n")

x_col <- colnames(data_filtered)[3]
y_col <- colnames(data_filtered)[4]
color_col <- colnames(data_filtered)[5]

subpop_names <- sub(".*\\.", "", colnames(data_filtered)[3:4])

ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color=.data[["Cluster"]])) +
  geom_point() +
  labs(x = subpop_names[[1]], y = subpop_names[[2]], title = subpop_names[1]) +
  theme_classic()

## Filter loci with 0 frequency in both subclones

data_filtered <- data_filtered[!(data_filtered[[3]] == 0 & data_filtered[[4]] == 0 & data_filtered[[5]] == 0),]
cat("After filtering loci with no mutations", dim(data_filtered), "\n")

## Count the number of loci with fixed alleles
cat("Loci having at least one fixed subpop", sum(data_filtered[[3]] == 0 |
                                                   data_filtered[[4]] == 0|
                                                   data_filtered[[5]] == 0), "\n\n")



## Format as allele frequency matrix for PopGenBound package -----

## Subset allele frequency columns
data_frequencies <- data_filtered[,3:(K+2)]
head(data_frequencies, 10)

list_loci <- lapply(1:nrow(data_frequencies), function(i){
  mat <- as.matrix(data_frequencies[i,])
  second_allele <- 1 - mat[1, ]  # Compute 1-x for each element in the first row
  new_mat <- rbind(mat, second_allele)
  rownames(new_mat) <- c("First", "Second")
  return(t(new_mat))
}
)

list_loci[[2]]

## Compute statistics with Diff() for the list of loci -------

loci <- list_loci[[3]]
print(loci)

subpop_12 <- loci[c(1,2),]
subpop_13 <- loci[c(1,3),]
subpop_23 <- loci[c(2,3),]

Diff(subpop_12)
Diff(subpop_13)
Diff(subpop_23)

# Extract rows 1 & 2
list_12 <- lapply(list_loci, function(mat) mat[c(1, 2), ])

# Extract rows 2 & 3
list_23 <- lapply(list_loci, function(mat) mat[c(2, 3), ])

# Extract rows 1 & 3
list_13 <- lapply(list_loci, function(mat) mat[c(1, 3), ])

Diff(loci)

D_loci = compute_Diff_3subpop(list_loci)

plot123 <- plot_stats(D_loci$D_123, K=3, title=glue("K=3: {subpop_names[1]}, {subpop_names[2]}, {subpop_names[3]}"))
plot12 <- plot_stats(D_loci$D_12, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[2]}"))
plot23 <- plot_stats(D_loci$D_23, K=2, title=glue("K=2: {subpop_names[2]}, {subpop_names[3]}"))
plot13 <- plot_stats(D_loci$D_13, K=2, title=glue("K=2: {subpop_names[1]}, {subpop_names[3]}"))

final_plot <- final_plot <- wrap_elements(plot123) /
                            wrap_elements(plot12) /
                            wrap_elements(plot23) /
                            wrap_elements(plot13)
final_plot

ggbounds(M=Diff(loci)$value[1],FST=Diff(loci)$value[2],GpST = Diff(loci)$value[3],D=Diff(loci)$value[4],K=nrow(loci))

Diff_loci <- lapply(list_loci, Diff)

compute_mean_stats(Diff_loci)
## Plotting --------------------------------------------------

lcnec.tib = tibble(M=unlist( sapply(Diff_loci,function(x){x[1,2]})),
                   FST=unlist( sapply(Diff_loci,function(x){x[2,2]})),
                   GpST=unlist( sapply(Diff_loci,function(x){x[3,2]})),
                   D=unlist( sapply(Diff_loci,function(x){x[4,2]})))

#lcnec.tib[apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]

lcnec.tib_clean <- lcnec.tib[!apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]



FST_mean <- mean(lcnec.tib_clean$FST, na.rm = TRUE)

FST_norm <- lcnec.tib_clean$FST / sapply(lcnec.tib_clean$M, function(m) Fup(K, m))
FST_norm_mean <- mean(FST_norm, na.rm=T)
#summary(FST_norm)

summary()# normalise stats by their max values for a given locus
# to do (or not?)

gglcnec = ggbounds_loc(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=2)

gglcnec_norm = ggbounds_loc_norm(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=2)

gg_all <- ggbounds_new(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=2)

(gg_all[[1]][[1]] + gg_all[[2]][[1]]) /
(gg_all[[1]][[2]] + gg_all[[2]][[2]]) /
(gg_all[[1]][[3]] + gg_all[[2]][[3]])

(gglcnec[[1]] + gglcnec_norm[[1]]) /
(gglcnec[[2]] + gglcnec_norm[[2]]) /
(gglcnec[[3]] + gglcnec_norm[[3]])

MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                      FST= Fup(K,seq(0.001,1-0.001,0.001)))

dplyr::tibble(M=mean(lcnec.tib_clean$M,na.rm=T),FST=mean(lcnec.tib_clean$FST,na.rm=T))


#mean_FST <- mean(lcnec.tib_clean$FST, na.rm = TRUE)

FST_norm <- lcnec.tib_clean$FST / sapply(lcnec.tib_clean$M, function(m) Fup(K, m))

MF_ST_tib <- dplyr::tibble(M=lcnec.tib_clean$M,FST_n=FST_norm)

mean_M <- mean(lcnec.tib_clean$M, na.rm = TRUE)
mean_FST_n <- mean(FST_norm, na.rm=T)

plot_FST_norm <-
  ggplot2::ggplot(MF_ST_tib, ggplot2::aes(x = M, y = FST_n)) +
  ggpointdensity::geom_pointdensity() +
  ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST_n),
                        ggplot2::aes(x = M, xend = M, y = 0, yend = 1),
                        col = "red", size = 0.8, linetype = "dashed") +
  ggplot2::geom_point(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                      col="red", pch=16,size=3,stroke=2) +

  ggplot2::coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
  ggplot2::xlab(expression(italic(M))) +
  ggplot2::ylab(expression(italic(F[ST]))) +
  ggplot2::theme_bw()


plot_FST_norm

ggbounds_loc = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M,na.rm=T)<0.5)*0.16-(mean(M,na.rm=T)>=0.5)*0.25
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                         GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )

  mean_M <- mean(M, na.rm = TRUE)
  mean_FST <- mean(FST, na.rm = TRUE)

  plotFST <-
    ggplot2::ggplot(dplyr::tibble(M=M,FST=FST), ggplot2::aes(x = M, y = FST)) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::geom_line(data = MFtmp, ggplot2::aes(x = M, y = FST)) +
    ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST),
                          ggplot2::aes(x = M, xend = M, y = 0, yend = 1),
                          col = "red", size = 0.8, linetype = "dashed") +
    ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, FST_mean = mean_FST),
                          ggplot2::aes(x = 0, xend = 1, y = FST_mean, yend = FST_mean),
                          col = "red", size = 0.8, linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
    ggplot2::xlab(expression(italic(M))) +
    ggplot2::ylab(expression(italic(F[ST]))) +
    ggplot2::theme_bw()


  if(!is.null(GpST)){
    mean_GpST <- mean(GpST, na.rm = TRUE)

    plotGpST <-
      ggplot2::ggplot(dplyr::tibble(M=M,GpST=GpST),ggplot2::aes(x=M,y=GpST)) +
      ggpointdensity::geom_pointdensity() +
      ggplot2::geom_line(data=MGptmp,ggplot2::aes(x=M,y=GpST)) +
      ggplot2::geom_segment(data=dplyr::tibble(M=mean_M,GpST=mean_GpST),
                            ggplot2::aes(x=M,xend=M,y=0,yend=1),
                            col="red", size=0.8,linetype = "dashed") +
      ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, GpST_mean = mean_GpST),
                            ggplot2::aes(x = 0, xend = 1, y = GpST_mean, yend = GpST_mean),
                            col = "red", size = 0.8, linetype = "dashed") +
      ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
      ggplot2::xlab(expression(italic(M))) +
      ggplot2::ylab(expression(italic(G[ST]))) +
      ggplot2::theme_bw()
  }
  else{plotGpST = NULL}

  if(!is.null(D)){
    mean_D <- mean(D, na.rm = TRUE)

    plotD <-
      ggplot2::ggplot(dplyr::tibble(M=M,D=D),ggplot2::aes(x=M,y=D)) +
      ggpointdensity::geom_pointdensity() +
      ggplot2::geom_line(data=MDtmp,ggplot2::aes(x=M,y=D)) +
      ggplot2::geom_segment(data=dplyr::tibble(M=mean_M,D=mean_D),
                            ggplot2::aes(x=M,xend=M,y=0,yend=1),
                            col="red", size=0.8,linetype = "dashed") +
      ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, D_mean = mean_D),
                            ggplot2::aes(x = 0, xend = 1, y = D_mean, yend = D_mean),
                            col = "red", size = 0.8, linetype = "dashed") +
      ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
      ggplot2::xlab(expression(italic(M))) +
      ggplot2::ylab(expression(italic(D))) +
      ggplot2::theme_bw()
  }
  else{plotD=NULL}

  if(length(M)>2) plot <- plot + ggplot2::scale_color_viridis_b()
  return(list(plotFST,plotGpST,plotD) )
}

ggbounds_loc_norm = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M,na.rm=T)<0.75)*0.16-(mean(M,na.rm=T)>=0.75)*0.16
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                         GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )
  # compute normalised FST
  FST_norm <- FST / sapply(M, function(m) Fup(K, m))
  MF_ST_tib <- dplyr::tibble(M=M,FST_n=FST_norm)

  mean_M <- mean(M, na.rm = TRUE)
  mean_FST_n <- mean(FST_norm, na.rm=T)

  plotFST_norm <-
    ggplot2::ggplot(MF_ST_tib, ggplot2::aes(x = M, y = FST_n)) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST_n),
                          ggplot2::aes(x = M, xend = M, y = 0, yend = 1),
                          col = "red", size = 0.8, linetype = "dashed") +
    ggplot2::geom_point(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                        col="red", pch=16,size=3,stroke=2) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                        ggplot2::aes(x=M,y=FST_n,
                                     label = paste0("mean FST=",
                                                    format(FST_n,digits=2))),
                        nudge_x = nudge,nudge_y=0,
                        col="red", size=3) +
    ggplot2::coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
    ggplot2::xlab(expression(italic(M))) +
    ggplot2::ylab(expression(italic(F[ST]))) +
    ggplot2::theme_bw()

  if(!is.null(GpST)){
    GST_norm <- GpST / sapply(M, function(m) Gpup(K, m))
    MG_ST_tib <- dplyr::tibble(M=M,GST_n=GST_norm)
    mean_GST_n <- mean(GST_norm, na.rm=T)

    plotGpST_norm <-
      ggplot2::ggplot(MG_ST_tib, ggplot2::aes(x = M, y = GST_n)) +
      ggpointdensity::geom_pointdensity() +
      ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, GST_n = mean_GST_n),
                            ggplot2::aes(x = M, xend = M, y = 0, yend = 1),
                            col = "red", size = 0.8, linetype = "dashed") +
      ggplot2::geom_point(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                          col="red", pch=16,size=3,stroke=2) +
      ggplot2::geom_label(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                          ggplot2::aes(x=M,y=GST_n,
                            label = paste0("mean GpST=",
                                      format(GST_n,digits=2))),
                            nudge_x = nudge,nudge_y=0,
                            col="red", size=3) +
      ggplot2::coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
      ggplot2::xlab(expression(italic(M))) +
      ggplot2::ylab(expression(italic(G[ST]))) +
      ggplot2::theme_bw()
  }
  else{
    plotGpST_norm = NULL
  }

  if(!is.null(D)){
    D_norm <- D / sapply(M, function(m) Dup(K, m))
    MD_tib <- dplyr::tibble(M=M,D_n=D_norm)
    mean_D_n <- mean(D_norm, na.rm=T)

    plotD_norm <-
      ggplot2::ggplot(MD_tib, ggplot2::aes(x = M, y = D_n)) +
      ggpointdensity::geom_pointdensity() +
      ggplot2::geom_segment(data = dplyr::tibble(M = mean_M, D_n = mean_D_n),
                            ggplot2::aes(x = M, xend = M, y = 0, yend = 1),
                            col = "red", size = 0.8, linetype = "dashed") +
      ggplot2::geom_point(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                          col="red", pch=16,size=3,stroke=2) +
      ggplot2::geom_label(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                          ggplot2::aes(x=M,y=D_n,
                                       label = paste0("mean D=",
                                                      format(D_n,digits=2))),
                          nudge_x = nudge,nudge_y=0,
                          col="red", size=3) +
      ggplot2::coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
      ggplot2::xlab(expression(italic(M))) +
      ggplot2::ylab(expression(italic(D))) +
      ggplot2::theme_bw()
  }
  else{
    plotD_norm=NULL
  }

  if(length(M)>2) plot <- plot + ggplot2::scale_color_viridis_b()
  return(list(plotFST_norm,plotGpST_norm,plotD_norm) )
}


ggbounds_new <- function(M,FST,GpST=NULL,D=NULL,K=2){
  gg1 <- ggbounds_loc(M, FST, GpST, D, K)
  gg2 <- ggbounds_loc_norm(M, FST, GpST, D, K)
  return(list(gg1, gg2))
}

