# Functions to process data ------------------------------------------------

filter_data <- function(df, K=2) {
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

plot_stats <- function(Diff_loci, K=2) {
  lcnec.tib = tibble(M=unlist( sapply(Diff_loci,function(x){x[1,2]})),
                     FST=unlist( sapply(Diff_loci,function(x){x[2,2]})),
                     GpST=unlist( sapply(Diff_loci,function(x){x[3,2]})),
                     D=unlist( sapply(Diff_loci,function(x){x[4,2]})))
  
  lcnec.tib[apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]
  
  lcnec.tib_clean <- lcnec.tib[!apply(is.nan(as.matrix(lcnec.tib)), 1, any), ]
  
  gglcnec = ggbounds(M=lcnec.tib_clean$M,FST=lcnec.tib_clean$FST,GpST=lcnec.tib_clean$GpST,D=lcnec.tib_clean$D,K=K)
  
  combined_plot <- gglcnec[[1]] + gglcnec[[2]]+ gglcnec[[3]]
  
  return(combined_plot)
}

save_plots <- function(df, combined_plot, K=2) {
  if (K == 2) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:4])
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]
    
    output_folder <- "../../../plots/"
    filename <- paste0(output_folder, pop1, "_vs_", pop2, ".svg")
    print(filename)
    ggsave(filename,
           plot = combined_plot,
           height = 3, width = 12)
  }
  else if (K == 3) {
    subpop_names <- sub(".*\\.", "", colnames(df)[3:5])
    pop1 <- subpop_names[1]
    pop2 <- subpop_names[2]
    pop3 <- subpop_names[3]
    
    output_folder <- "../../../plots/"
    filename <- paste0(output_folder, pop1, "_vs_", pop2, "_vs_", pop3, ".svg")
    print(filename)
    ggsave(filename,
           plot = combined_plot,
           height = 3, width = 12)
  }
}

run_lcnec_plotting <- function(data_frame_list, K=2) {
  for (df in data_frame_list) {
    data_clean <- filter_data(df, K)
    list_freq <- make_popgen_input(data_clean[,3:(2+K)])
    Diff_loci <- lapply(list_freq, Diff)
    combined_plot <- plot_stats(Diff_loci, K)
    save_plots(df, combined_plot, K)
  }
}
