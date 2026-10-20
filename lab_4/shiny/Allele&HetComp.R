# Population genetics practice app
# Simulates diploid genotypes at one locus; students work out the number of
# alleles, the allele frequencies and the observed heterozygosity, and the app
# checks their answers.
#
# Run with: shiny::runApp("app.R")

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

  list(
    alleles = alleles,   # all alleles in the population (A, B, ... up to k)
    data = data.frame(
      Individual = seq_len(n),
      Allele1    = first,
      Allele2    = second,
      Genotype   = paste(first, second, sep = "/"),
      stringsAsFactors = FALSE
    )
  )
}

n_alleles <- function(pop) length(unique(c(pop$data$Allele1, pop$data$Allele2)))
obs_het   <- function(pop) mean(pop$data$Allele1 != pop$data$Allele2)

# Number of copies of each allele in the sample (homozygotes contribute two)
allele_counts <- function(pop) {
  counts <- table(factor(c(pop$data$Allele1, pop$data$Allele2),
                         levels = pop$alleles))
  setNames(as.vector(counts), pop$alleles)
}

# Allele frequency = copies of the allele / (2 x number of individuals)
allele_freqs <- function(pop) allele_counts(pop) / (2 * nrow(pop$data))

# Answers within this distance of the true value count as correct
# (allows for rounding to two decimal places)
TOLERANCE <- 0.005

within_tol <- function(answer, truth) abs(answer - truth) <= TOLERANCE + 1e-9

# ---- UI ---------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Genotypes, alleles and heterozygosity"),

  sidebarLayout(
     sidebarPanel(
    #   h4("1. Simulate a population sample"),
    #   sliderInput("n", "Number of individuals sampled",
    #               min = 5, max = 100, value = 20, step = 1),
    #   sliderInput("k", "Alleles segregating in the population",
    #               min = 2, max = 8, value = 4, step = 1),
    #   sliderInput("f", "Inbreeding coefficient (F)",
    #               min = 0, max = 1, value = 0, step = 0.05),
    #   actionButton("new", "New sample", class = "btn-primary"),

      hr(),

      h4("2. Enter your answers"),
      numericInput("ans_alleles",
                   "Number of alleles observed in the sample",
                   value = NA, min = 1, step = 1),

      tags$label("Allele frequencies (2 decimal places)"),
      helpText("Enter 0 for any allele that does not appear in the sample."),
      uiOutput("freq_inputs"),

      numericInput("ans_het",
                   "Observed heterozygosity (proportion, 2 decimal places)",
                   value = NA, min = 0, max = 1, step = 0.01),
      actionButton("check", "Check answers", class = "btn-success"),
      actionButton("reveal", "Show solution"),

      hr(),
      uiOutput("feedback")
    ),

    mainPanel(
      h4("Sampled genotypes"),
      p("Each individual is diploid, so it carries two alleles at this locus.",
        "Count the different alleles present, work out the frequency of each",
        "allele, and find the proportion of individuals that are heterozygous."),
      checkboxInput("show_counts", "Hint: show genotype counts", value = FALSE),
      conditionalPanel(
        "input.show_counts",
        tableOutput("genotype_counts")
      ),
      div(style = "max-height: 500px; overflow-y: auto;",
          tableOutput("genotypes"))
    )
  )
)

# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {

  pop      <- reactiveVal(NULL)
  feedback <- reactiveVal(NULL)

  # Draw a new sample at start-up and whenever the button is pressed
  observeEvent(input$new, {
    pop(simulate_genotypes(20, sample(2:5, 1), 0))
    feedback(NULL)
    updateNumericInput(session, "ans_alleles", value = NA)
    updateNumericInput(session, "ans_het", value = NA)
  }, ignoreNULL = FALSE)

  # One frequency box per allele in the population. Rebuilt (and so cleared)
  # every time a new sample is drawn.
  output$freq_inputs <- renderUI({
    req(pop())
    boxes <- lapply(pop()$alleles, function(a) {
      div(style = "display: inline-block; width: 80px; margin-right: 6px;",
          numericInput(paste0("freq_", a), label = a,
                       value = NA, min = 0, max = 1, step = 0.01))
    })
    div(boxes)
  })

  output$genotypes <- renderTable({
    req(pop())
    pop()$data[, c("Individual", "Genotype")]
  }, digits = 0, striped = TRUE)

  output$genotype_counts <- renderTable({
    req(pop())
    as.data.frame(table(Genotype = pop()$data$Genotype),
                  responseName = "Count",
                  stringsAsFactors = FALSE)
  }, digits = 0)

  # Check the student's answers: a table of each answer next to the correct
  # answer, with a tick or cross
  observeEvent(input$check, {
    req(pop())

    true_alleles <- n_alleles(pop())
    true_freqs   <- allele_freqs(pop())
    true_het     <- obs_het(pop())

    ans_freqs <- sapply(pop()$alleles, function(a) {
      v <- input[[paste0("freq_", a)]]
      if (is.null(v)) NA_real_ else v
    })

    # One row per question
    rows <- list(
      list(label   = "Number of alleles",
           answer  = input$ans_alleles,
           correct = true_alleles,
           ok      = isTRUE(input$ans_alleles == true_alleles),
           digits  = 0)
    )
    for (a in pop()$alleles) {
      rows[[length(rows) + 1]] <- list(
        label   = paste("Frequency of", a),
        answer  = ans_freqs[[a]],
        correct = true_freqs[[a]],
        ok      = isTRUE(within_tol(ans_freqs[[a]], true_freqs[[a]])),
        digits  = 2)
    }
    rows[[length(rows) + 1]] <- list(
      label   = "Observed heterozygosity",
      answer  = input$ans_het,
      correct = true_het,
      ok      = isTRUE(within_tol(input$ans_het, true_het)),
      digits  = 2)

    fmt <- function(x, digits) {
      if (is.null(x) || is.na(x)) "—" else formatC(x, format = "f", digits = digits)
    }
    cell <- function(..., style = "") {
      tags$td(style = paste("padding: 3px 8px; border-bottom: 1px solid #ddd;", style), ...)
    }
    head_cell <- function(x) {
      tags$th(style = "padding: 3px 8px; border-bottom: 2px solid #999;", x)
    }

    results <- tags$table(
      style = "width: 100%; margin-bottom: 10px;",
      tags$tr(head_cell(""), head_cell("Your answer"),
              head_cell("Correct answer"), head_cell("")),
      lapply(rows, function(r) {
        blank <- is.null(r$answer) || is.na(r$answer)
        mark <- if (blank) {
          span(style = "color: #888;", "no answer")
        } else if (r$ok) {
          span(style = "color: #2e7d32; font-weight: bold; font-size: 1.2em;", "✓")
        } else {
          span(style = "color: #c62828; font-weight: bold; font-size: 1.2em;", "✗")
        }
        tags$tr(cell(r$label),
                cell(fmt(r$answer, r$digits)),
                cell(fmt(r$correct, r$digits)),
                cell(mark, style = "text-align: center;"))
      })
    )

    n_right <- sum(vapply(rows, function(r) r$ok, logical(1)))
    score <- tags$p(strong(sprintf("Score: %d / %d", n_right, length(rows))))

    # Frequencies must add up to 1; a useful cue when they don't
    sum_note <- NULL
    if (!any(is.na(ans_freqs)) && abs(sum(ans_freqs) - 1) > 0.02) {
      sum_note <- tags$p(style = "color: #c62828;",
        sprintf("Your allele frequencies add up to %.2f; they should sum to 1.",
                sum(ans_freqs)))
    }

    feedback(tagList(
      score,
      results,
      sum_note,
      helpText("Answers within 0.005 of the correct value count as correct.",
               "Click \"Show solution\" for the working.")
    ))
  })

  # Reveal the worked solution
  observeEvent(input$reveal, {
    req(pop())
    d       <- pop()$data
    n_het   <- sum(d$Allele1 != d$Allele2)
    n_ind   <- nrow(d)
    present <- sort(unique(c(d$Allele1, d$Allele2)))
    counts  <- allele_counts(pop())
    freqs   <- allele_freqs(pop())

    cell <- function(...) tags$td(style = "padding: 2px 10px;", ...)
    freq_table <- tags$table(
      style = "margin-bottom: 10px;",
      tags$tr(tags$th("Allele"), tags$th("Copies"), tags$th("Frequency")),
      lapply(names(counts), function(a) {
        tags$tr(cell(a),
                cell(counts[[a]], " / ", 2 * n_ind),
                cell(sprintf("%.2f", freqs[[a]])))
      })
    )

    feedback(tagList(
      tags$p(strong("Solution")),
      tags$p("Alleles present: ", paste(present, collapse = ", "),
             " → ", strong(length(present)), " alleles."),
      tags$p("Allele frequencies: copies of each allele divided by 2N = ",
             2 * n_ind, " allele copies."),
      freq_table,
      tags$p("Heterozygotes: ", n_het, " of ", n_ind, " individuals → ",
             "observed heterozygosity = ", n_het, " / ", n_ind, " = ",
             strong(sprintf("%.2f", n_het / n_ind)), ".")
    ))
  })

  output$feedback <- renderUI(feedback())
}

shinyApp(ui, server)
