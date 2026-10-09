library(tidyverse)
library(related)
library(pheatmap)

allele_frequencies <- function(file) {
  
  # Read the data
  dat <- read_csv(file)
  
  # Find the "a" copy of each locus
  allele_a <- names(dat)[str_ends(names(dat), "a")]
  
  # Calculate allele frequencies for each locus
  results <- map_dfr(allele_a, function(col_a) {
    
    # Find the matching "b" column
    col_b <- str_replace(col_a, "a$", "b")
    
    # Get the locus name
    locus <- str_remove(col_a, "a$")
    
    dat %>%
      select(all_of(c(col_a, col_b))) %>%
      
      # Put both allele copies into one column
      pivot_longer(
        cols = everything(),
        values_to = "allele"
      ) %>%
      
      # Remove missing alleles
      filter(!is.na(allele)) %>%
      
      # Count each allele
      count(allele) %>%
      
      # Calculate allele frequency
      mutate(
        frequency = n / sum(n),
        locus = locus
      ) %>%
      
      select(locus, allele, n, frequency)
  })
  
  return(results)
}

locus_statistics <- function(file) {
  
  # Read the data
  dat <- read_csv(file)
  
  # Find the "a" copy of each locus
  allele_a <- names(dat)[str_ends(names(dat), "a")]
  
  # Calculate statistics for each locus
  results <- map_dfr(allele_a, function(col_a) {
    
    # Find the matching "b" column
    col_b <- str_replace(col_a, "a$", "b")
    
    # Get the locus name
    locus <- str_remove(col_a, "a$")
    
    # Keep individuals with complete genotypes
    genotypes <- dat %>%
      select(all_of(c(col_a, col_b))) %>%
      filter(!is.na(.data[[col_a]]),
             !is.na(.data[[col_b]]))
    
    # Number of samples
    N <- nrow(genotypes)
    
    # Put both allele copies into one vector
    alleles <- c(genotypes[[col_a]], genotypes[[col_b]])
    
    # Count alleles and calculate allele frequencies
    allele_counts <- table(alleles)
    allele_freq <- allele_counts / sum(allele_counts)
    
    # Number of different alleles
    num_alleles <- length(allele_freq)
    
    # Effective number of alleles
    effective_alleles <- 1 / sum(allele_freq^2)
    
    # Shannon's information index
    shannon <- -sum(allele_freq * log(allele_freq))
    
    # Observed heterozygosity
    Ho <- mean(genotypes[[col_a]] != genotypes[[col_b]])
    
    # Expected heterozygosity
    He <- 1 - sum(allele_freq^2)
    
    # Unbiased expected heterozygosity
    # N = number of diploid individuals
    # Therefore, 2N = number of allele copies
    if (N > 0 && He > 0) {
      uHe <- (2 * N / (2 * N - 1)) * He
    } else {
      uHe <- NA
    }
    
    # Fixation index
    if (He > 0) {
      Fis <- (He - Ho) / He
    } else {
      Fis <- NA
    }
    
    # Return statistics for this locus
    tibble(
      locus = locus,
      samples = N,
      different_alleles = num_alleles,
      effective_alleles = effective_alleles,
      shannons_index = shannon,
      observed_heterozygosity = Ho,
      expected_heterozygosity = He,
      unbiased_heterozygosity = uHe,
      fixation_index = Fis
    )
  })
  
  return(results)
}

pairwise_relatedness <- function(file) {
  
  # Read the original CSV
  dat <- read_csv(file)
  
  # Find the genotype columns
  # They should end in "a" or "b"
  genotype_cols <- names(dat)[
    str_ends(names(dat), "a") | str_ends(names(dat), "b")
  ]
  
  # Keep the individual ID and genotype data
  genotype_data <- dat %>%
    select("individual ID", all_of(genotype_cols))
  
  # Make sure individual IDs are characters
  genotype_data$"individual ID" <- as.character(
    genotype_data$"individual ID"
  )
  
  # Replace missing values with 0
  genotype_data <- genotype_data %>%
    mutate(
      across(
        all_of(genotype_cols),
        ~ replace_na(.x, 0)
      )
    )
  
  # related requires the genotype data to be in a temporary
  # text file with no column names
  temp_file <- tempfile(fileext = ".txt")
  
  write.table(
    genotype_data,
    file = temp_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  # Calculate Queller & Goodnight relatedness
  results <- coancestry(
    temp_file,
    quellergt = 1
  )
  
  # Remove the temporary file
  unlink(temp_file)
  
  # Get the pairwise relatedness results
  relatedness <- results$relatedness
  
  # The Queller & Goodnight estimate is column 10
  relatedness <- relatedness %>%
    select(
      ID1 = 2,
      ID2 = 3,
      relatedness = 10
    )
  
  # Get all individual IDs
  IDs <- genotype_data$"individual ID"
  
  # Create an empty relatedness matrix
  relatedness_matrix <- matrix(
    NA,
    nrow = length(IDs),
    ncol = length(IDs),
    dimnames = list(IDs, IDs)
  )
  
  # Put the pairwise estimates into the matrix
  for (i in 1:nrow(relatedness)) {
    
    id1 <- as.character(relatedness$ID1[i])
    id2 <- as.character(relatedness$ID2[i])
    value <- relatedness$relatedness[i]
    
    relatedness_matrix[id1, id2] <- value
    relatedness_matrix[id2, id1] <- value
  }
  
  # Relatedness of an individual with itself
  diag(relatedness_matrix) <- 1
  
  # Define a contrasting color scale
max_abs <- max(abs(relatedness_matrix), na.rm = TRUE)

my_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(100)

my_breaks <- seq(-max_abs, max_abs, length.out = 101)

# Cluster individuals based on their relatedness profiles
dist_matrix <- as.dist(
  1 - relatedness_matrix
)

hc <- hclust(
  dist_matrix,
  method = "average"
)

# Plot with clustered individuals
pheatmap(
  relatedness_matrix,
  color = my_colors,
  breaks = my_breaks,
  cluster_rows = hc,
  cluster_cols = hc,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  display_numbers = FALSE,
  main = "Pairwise Relatedness (Queller & Goodnight)",
  fontsize = 10
)
  
  # Return the matrix
  return(relatedness_matrix)
}