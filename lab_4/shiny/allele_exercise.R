library(shiny)

# ---- Simulation -------------------------------------------------------------

# n: number of individuals sampled
# k: number of alleles segregating in the population
# f: inbreeding coefficient (0 = Hardy-Weinberg proportions)
simulate_genotypes <- function(n, k, f = 0) {
  alleles <- LETTERS[seq_len(k)]
  
  # Random allele frequencies that sum to 1
  p <- rgamma(k, shape = 1)
  p <- p / sum(p)
  
  a1 <- sample(alleles, n, replace = TRUE, prob = p)
  a2 <- sample(alleles, n, replace = TRUE, prob = p)
  
  # With probability f, the two alleles are identical by descent
  ibd <- runif(n) < f
  a2[ibd] <- a1[ibd]
  
  # Write each genotype in alphabetical order (A/C, never C/A)
  first  <- pmin(a1, a2)
  second <- pmax(a1, a2)
  
  data.frame(
    Individual = seq_len(n),
    Allele1    = first,
    Allele2    = second,
    Genotype   = paste(first, second, sep = "/"),
    stringsAsFactors = FALSE
  )
}

n_alleles <- function(pop) {
  length(unique(c(pop$Allele1, pop$Allele2)))
}

obs_het <- function(pop) {
  mean(pop$Allele1 != pop$Allele2)
}

# Frequency of allele A in the sampled chromosomes
allele_freq_A <- function(pop) {
  sum(c(pop$Allele1, pop$Allele2) == "A") / (2 * nrow(pop))
}

# Answers within this distance of the true value count as correct
# (allows for rounding to two decimal places)
ANSWER_TOLERANCE <- 0.005

# ---- UI ---------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Genotypes, alleles and heterozygosity"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Simulate a population sample"),
      
      sliderInput(
        "n",
        "Number of individuals sampled",
        min = 5, max = 100, value = 20, step = 1
      ),
      
      sliderInput(
        "k",
        "Alleles segregating in the population",
        min = 2, max = 8, value = 4, step = 1
      ),
      
      sliderInput(
        "f",
        "Inbreeding coefficient (F)",
        min = 0, max = 1, value = 0, step = 0.05
      ),
      
      actionButton("new", "New sample", class = "btn-primary"),
      
      hr(),
      
      h4("2. Enter your answers"),
      
      numericInput(
        "ans_alleles",
        "Number of alleles observed in the sample",
        value = NA,
        min = 1,
        step = 1
      ),
      
      numericInput(
        "ans_het",
        "Observed heterozygosity (proportion, 2 decimal places)",
        value = NA,
        min = 0,
        max = 1,
        step = 0.01
      ),
      
      numericInput(
        "ans_freq_A",
        "Allele frequency of A (proportion, 2 decimal places)",
        value = NA,
        min = 0,
        max = 1,
        step = 0.01
      ),
      
      actionButton("check", "Check answers", class = "btn-success"),
      actionButton("reveal", "Show solution"),
      
      hr(),
      
      uiOutput("feedback")
    ),
    
    mainPanel(
      h4("Sampled genotypes"),
      
      p(
        "Each individual is diploid, so it carries two alleles at this locus.",
        "Count the different alleles present, work out the proportion of",
        "individuals that are heterozygous, and calculate the frequency of",
        "allele A among all sampled chromosomes."
      ),
      
      checkboxInput(
        "show_counts",
        "Hint: show genotype counts",
        value = FALSE
      ),
      
      conditionalPanel(
        "input.show_counts",
        tableOutput("genotype_counts")
      ),
      
      div(
        style = "max-height: 500px; overflow-y: auto;",
        tableOutput("genotypes")
      )
    )
  )
)

# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {
  
  pop <- reactiveVal(NULL)
  feedback <- reactiveVal(NULL)
  
  # Draw a new sample at start-up and whenever the button is pressed
  observeEvent(input$new, {
    pop(simulate_genotypes(input$n, input$k, input$f))
    feedback(NULL)
    
    updateNumericInput(
      session,
      "ans_alleles",
      value = NA
    )
    
    updateNumericInput(
      session,
      "ans_het",
      value = NA
    )
    
    updateNumericInput(
      session,
      "ans_freq_A",
      value = NA
    )
  }, ignoreNULL = FALSE)
  
  output$genotypes <- renderTable({
    req(pop())
    
    pop()[, c("Individual", "Genotype")]
  }, digits = 0, striped = TRUE)
  
  output$genotype_counts <- renderTable({
    req(pop())
    
    counts <- as.data.frame(
      table(Genotype = pop()$Genotype),
      responseName = "Count",
      stringsAsFactors = FALSE
    )
    
    counts
  }, digits = 0)
  
  # Check the student's answers
  observeEvent(input$check, {
    req(pop())
    
    true_alleles <- n_alleles(pop())
    true_het <- obs_het(pop())
    true_freq_A <- allele_freq_A(pop())
    
    result_line <- function(ok, label, hint) {
      if (ok) {
        tags$p(
          style = "color: #2e7d32;",
          strong("\u2713 "),
          label,
          ": correct."
        )
      } else {
        tags$p(
          style = "color: #c62828;",
          strong("\u2717 "),
          label,
          ": not correct. ",
          hint
        )
      }
    }
    
    # Number of alleles
    if (is.na(input$ans_alleles)) {
      alleles_msg <- tags$p(
        style = "color: #555;",
        "Number of alleles: please enter an answer."
      )
    } else {
      alleles_msg <- result_line(
        isTRUE(input$ans_alleles == true_alleles),
        "Number of alleles",
        "Count each different letter once, across all genotypes."
      )
    }
    
    # Observed heterozygosity
    if (is.na(input$ans_het)) {
      het_msg <- tags$p(
        style = "color: #555;",
        "Observed heterozygosity: please enter an answer."
      )
    } else if (input$ans_het > 1) {
      het_msg <- tags$p(
        style = "color: #c62828;",
        strong("\u2717 "),
        "Observed heterozygosity is a proportion, so it must",
        " lie between 0 and 1."
      )
    } else {
      het_msg <- result_line(
        abs(input$ans_het - true_het) <= ANSWER_TOLERANCE + 1e-9,
        "Observed heterozygosity",
        "Divide the number of heterozygous individuals by the total number of individuals."
      )
    }
    
    # Allele frequency of A
    if (is.na(input$ans_freq_A)) {
      freq_A_msg <- tags$p(
        style = "color: #555;",
        "Allele frequency of A: please enter an answer."
      )
    } else if (input$ans_freq_A < 0 || input$ans_freq_A > 1) {
      freq_A_msg <- tags$p(
        style = "color: #c62828;",
        strong("\u2717 "),
        "Allele frequency must lie between 0 and 1."
      )
    } else {
      freq_A_msg <- result_line(
        abs(input$ans_freq_A - true_freq_A) <= ANSWER_TOLERANCE + 1e-9,
        "Allele frequency of A",
        "Count the number of A alleles and divide by the total number of sampled chromosomes (2 × number of individuals)."
      )
    }
    
    feedback(
      tagList(
        alleles_msg,
        het_msg,
        freq_A_msg
      )
    )
  })
  
  # Reveal the worked solution
  observeEvent(input$reveal, {
    req(pop())
    
    n_het <- sum(pop()$Allele1 != pop()$Allele2)
    n_ind <- nrow(pop())
    
    present <- sort(
      unique(c(pop()$Allele1, pop()$Allele2))
    )
    
    n_A <- sum(
      c(pop()$Allele1, pop()$Allele2) == "A"
    )
    
    n_chromosomes <- 2 * n_ind
    freq_A <- n_A / n_chromosomes
    
    feedback(
      tagList(
        tags$p(strong("Solution")),
        
        tags$p(
          "Alleles present: ",
          paste(present, collapse = ", "),
          " \u2192 ",
          strong(length(present)),
          " alleles."
        ),
        
        tags$p(
          "Heterozygotes: ",
          n_het,
          " of ",
          n_ind,
          " individuals \u2192 ",
          "observed heterozygosity = ",
          n_het,
          " / ",
          n_ind,
          " = ",
          strong(sprintf("%.2f", n_het / n_ind)),
          "."
        ),
        
        tags$p(
          "Allele A: ",
          n_A,
          " of ",
          n_chromosomes,
          " sampled chromosomes \u2192 ",
          "allele frequency of A = ",
          n_A,
          " / ",
          n_chromosomes,
          " = ",
          strong(sprintf("%.2f", freq_A)),
          "."
        )
      )
    )
  })
  
  output$feedback <- renderUI(feedback())
}

shinyApp(ui, server)