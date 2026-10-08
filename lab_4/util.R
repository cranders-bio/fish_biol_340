library(tidyverse)

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