if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("shinyjs", quietly = TRUE)) install.packages("shinyjs")


library(shiny)
library(shinyjs)
library(tidyverse)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Population-Specific Genetics Parameter Validator"),
  sidebarLayout(
    sidebarPanel(
      h3("Population 1 Genotypes"),
      fluidRow(
        column(6, numericInput("p1_AA", "AA:", 0, min = 0)),
        column(6, numericInput("p1_Aa", "Aa:", 0, min = 0)),
        column(6, numericInput("p1_aa", "aa:", 0, min = 0))
      ),
      
      h3("Population 2 Genotypes"),
      fluidRow(
        column(6, numericInput("p2_AA", "AA:", 0, min = 0)),
        column(6, numericInput("p2_Aa", "Aa:", 0, min = 0)),
        column(6, numericInput("p2_aa", "aa:", 0, min = 0))
      ),
      
      h3("Enter Your Calculations"),
      fluidRow(
        column(6,
               h4("Population 1"),
               numericInput("user_p1_Hobs", "Hobs:", 0, min = 0, max = 1, step = 0.01),
               numericInput("user_p1_Hexp", "Hexp:", 0, min = 0, max = 1, step = 0.01),
               numericInput("user_p1_FIS", "FIS:", 0, min = -1, max = 1, step = 0.01)
        ),
        column(6,
               h4("Population 2"),
               numericInput("user_p2_Hobs", "Hobs:", 0, min = 0, max = 1, step = 0.01),
               numericInput("user_p2_Hexp", "Hexp:", 0, min = 0, max = 1, step = 0.01),
               numericInput("user_p2_FIS", "FIS:", 0, min = -1, max = 1, step = 0.01)
        )
      ),
      numericInput("user_FST", "FST:", 0, min = 0, max = 1, step = 0.01),
      actionButton("validate", "Validate Parameters")
    ),
    
    mainPanel(
      h3("Validation Results"),
      tableOutput("resultsTable"),
      htmlOutput("validationMessage")
    )
  )
)

server <- function(input, output) {
  calculate_params <- reactive({
    # Population 1 calculations
    p1_total <- input$p1_AA + input$p1_Aa + input$p1_aa
    p1_p <- ifelse(p1_total > 0, (2*input$p1_AA + input$p1_Aa)/(2*p1_total), NA)
    p1_Hobs <- ifelse(p1_total > 0, input$p1_Aa/p1_total, NA)
    p1_Hexp <- ifelse(!is.na(p1_p), 2 * p1_p * (1 - p1_p), NA)
    p1_FIS <- ifelse(!is.na(p1_Hexp) && p1_Hexp > 0, (p1_Hexp - p1_Hobs)/p1_Hexp, NA)
    
    # Population 2 calculations
    p2_total <- input$p2_AA + input$p2_Aa + input$p2_aa
    p2_p <- ifelse(p2_total > 0, (2*input$p2_AA + input$p2_Aa)/(2*p2_total), NA)
    p2_Hobs <- ifelse(p2_total > 0, input$p2_Aa/p2_total, NA)
    p2_Hexp <- ifelse(!is.na(p2_p), 2 * p2_p * (1 - p2_p), NA)
    p2_FIS <- ifelse(!is.na(p2_Hexp) && p2_Hexp > 0, (p2_Hexp - p2_Hobs)/p2_Hexp, NA)
    
    # FST calculation
    total_AA <- input$p1_AA + input$p2_AA
    total_Aa <- input$p1_Aa + input$p2_Aa
    total_aa <- input$p1_aa + input$p2_aa
    total_indiv <- p1_total + p2_total
    
    p_total <- ifelse(total_indiv > 0, 
                      (2*total_AA + total_Aa)/(2*total_indiv), 
                      NA)
    
    HT <- ifelse(!is.na(p_total), 2 * p_total * (1 - p_total), NA)
    HS <- ifelse(!is.na(p1_Hexp) && !is.na(p2_Hexp), 
                 (p1_Hexp * p1_total + p2_Hexp * p2_total)/(p1_total + p2_total), 
                 NA)
    FST <- ifelse(!is.na(HT) && !is.na(HS) && HT > 0, (HT - HS)/HT, NA)
    
    list(
      p1 = list(Hobs = p1_Hobs, Hexp = p1_Hexp, FIS = p1_FIS),
      p2 = list(Hobs = p2_Hobs, Hexp = p2_Hexp, FIS = p2_FIS),
      FST = FST
    )
  })
  
  observeEvent(input$validate, {
    actual <- calculate_params()
    user_entries <- c(
      input$user_p1_Hobs, input$user_p1_Hexp, input$user_p1_FIS,
      input$user_p2_Hobs, input$user_p2_Hexp, input$user_p2_FIS,
      input$user_FST
    )
    
    calculated <- c(
      actual$p1$Hobs, actual$p1$Hexp, actual$p1$FIS,
      actual$p2$Hobs, actual$p2$Hexp, actual$p2$FIS,
      actual$FST
    )
    
    results <- data.frame(
      Parameter = c("Pop1 Hobs", "Pop1 Hexp", "Pop1 FIS",
                    "Pop2 Hobs", "Pop2 Hexp", "Pop2 FIS",
                    "FST"),
      Your_Value = user_entries,
      Calculated_Value = sapply(calculated, function(x) ifelse(is.na(x), "NA", round(x, 3))),
      Match = mapply(function(x, y) {
        if(is.na(x) || is.na(y)) return(FALSE)
        abs(x - y) < 0.01
      }, user_entries, calculated)
    )
    
    output$resultsTable <- renderTable({
      results
    })
    
    validation_text <- ifelse(all(results$Match, na.rm = TRUE),
                              "<span style='color:green;font-weight:bold'>All parameters correct! ✔</span>",
                              "<span style='color:red;font-weight:bold'>Some parameters incorrect ✘</span>"
    )
    
    output$validationMessage <- renderUI({
      HTML(validation_text)
    })
  })
}

shinyApp(ui, server)



#######################################

# This is the code for the shiny app that I have created. The app is a genetics parameter validator that takes in user inputs for the genotypes of two populations and their respective genetic parameters (Hobs, Hexp, FIS) and FST. The app then calculates the genetic parameters and FST based on the user inputs and compares them to the calculated values. The app then displays a table with the user inputs, calculated values, and whether they match or not, along with a validation message indicating if all parameters are correct or not.

library(shiny)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Population Genetics Calculator with Random Data"),
  sidebarLayout(
    sidebarPanel(
      actionButton("generate", "Generate New Populations"),
      h4("Sample Data (100 individuals per population)"),
      tableOutput("genoTable"),
      hr(),
      h3("Enter Your Calculations:"),
      fluidRow(
        column(6,
               h4("Population 1"),
               numericInput("p1_Hobs", "Observed Het (Hobs)", 0, step = 0.01),
               numericInput("p1_Hexp", "Expected Het (Hexp)", 0, step = 0.01),
               numericInput("p1_FIS", "FIS", 0, step = 0.01)
        ),
        column(6,
               h4("Population 2"),
               numericInput("p2_Hobs", "Observed Het (Hobs)", 0, step = 0.01),
               numericInput("p2_Hexp", "Expected Het (Hexp)", 0, step = 0.01),
               numericInput("p2_FIS", "FIS", 0, step = 0.01)
        )
      ),
      fluidRow(
        column(6, numericInput("user_Hs", "Hs (Within-pop)", 0, step = 0.01)),
        column(6, numericInput("user_Ht", "Ht (Total)", 0, step = 0.01))
      ),
      numericInput("user_FST", "FST", 0, step = 0.01),
      actionButton("validate", "Check Answers")
    ),
    mainPanel(
      h3("Validation Results"),
      tableOutput("validationTable"),
      htmlOutput("summary")
    )
  )
)

server <- function(input, output) {
  rv <- reactiveValues(
    data = NULL,
    trueValues = NULL
  )
  
  observeEvent(input$generate, {
    # Generate random genotype data
    generate_pop <- function() {
      p <- runif(1, 0.1, 0.9)
      FIS <- runif(1, -0.5, 0.5)
      Hexp <- 2 * p * (1 - p)
      Hobs <- Hexp * (1 - FIS)
      
      # Ensure positive counts
      AA <- max(0, round((p^2 + FIS*p*(1-p)) * 100))
      Aa <- max(0, round(Hobs * 100))
      aa <- max(0, 100 - AA - Aa)
      
      # Adjust counts if necessary to ensure sum is 100
      total <- AA + Aa + aa
      if (total != 100) {
        diff <- 100 - total
        if (diff > 0) {
          # Add to the largest count
          largest <- which.max(c(AA, Aa, aa))
          if (largest == 1) AA <- AA + diff
          else if (largest == 2) Aa <- Aa + diff
          else aa <- aa + diff
        } else {
          # Subtract from the largest count
          largest <- which.max(c(AA, Aa, aa))
          if (largest == 1) AA <- max(0, AA + diff)
          else if (largest == 2) Aa <- max(0, Aa + diff)
          else aa <- max(0, aa + diff)
        }
      }
      
      # Recalculate parameters based on adjusted counts
      total <- AA + Aa + aa
      p_adj <- (2*AA + Aa) / (2*total)
      Hobs_adj <- Aa / total
      Hexp_adj <- 2 * p_adj * (1 - p_adj)
      FIS_adj <- ifelse(Hexp_adj > 0, (Hexp_adj - Hobs_adj) / Hexp_adj, 0)
      
      list(
        counts = c(AA = AA, Aa = Aa, aa = aa),
        p = p_adj,
        Hobs = Hobs_adj,
        Hexp = Hexp_adj,
        FIS = FIS_adj
      )
    }
    
    pop1 <- generate_pop()
    pop2 <- generate_pop()
    
    # Calculate global parameters
    p_total <- (2*(pop1$counts["AA"] + pop2$counts["AA"]) + 
                  pop1$counts["Aa"] + pop2$counts["Aa"]) / 400
    Ht <- 2 * p_total * (1 - p_total)
    Hs <- (pop1$Hexp + pop2$Hexp)/2
    FST <- (Ht - Hs)/Ht
    
    rv$data <- list(pop1 = pop1, pop2 = pop2)
    rv$trueValues <- list(
      p1_Hobs = pop1$Hobs,
      p1_Hexp = pop1$Hexp,
      p1_FIS = pop1$FIS,
      p2_Hobs = pop2$Hobs,
      p2_Hexp = pop2$Hexp,
      p2_FIS = pop2$FIS,
      Hs = Hs,
      Ht = Ht,
      FST = FST
    )
  })
  
  output$genoTable <- renderTable({
    if (!is.null(rv$data)) {
      data.frame(
        Population = c("1", "2"),
        AA = c(rv$data$pop1$counts["AA"], rv$data$pop2$counts["AA"]),
        Aa = c(rv$data$pop1$counts["Aa"], rv$data$pop2$counts["Aa"]),
        aa = c(rv$data$pop1$counts["aa"], rv$data$pop2$counts["aa"])
      )
    }
  }, digits = 0)
  
  observeEvent(input$validate, {
    req(rv$trueValues)
    user <- list(
      p1_Hobs = input$p1_Hobs,
      p1_Hexp = input$p1_Hexp,
      p1_FIS = input$p1_FIS,
      p2_Hobs = input$p2_Hobs,
      p2_Hexp = input$p2_Hexp,
      p2_FIS = input$p2_FIS,
      Hs = input$user_Hs,
      Ht = input$user_Ht,
      FST = input$user_FST
    )
    
    results <- data.frame(
      Parameter = c("Pop1 Hobs", "Pop1 Hexp", "Pop1 FIS",
                    "Pop2 Hobs", "Pop2 Hexp", "Pop2 FIS",
                    "Hs", "Ht", "FST"),
      Your_Value = unlist(user),
      True_Value = unlist(rv$trueValues),
      Match = abs(unlist(user) - unlist(rv$trueValues)) < 0.01
    )
    
    output$validationTable <- renderTable({
      results %>% mutate(
        True_Value = round(True_Value, 3),
        Match = ifelse(Match, "✓", "✗")
      )
    })
    
    correct <- sum(results$Match)
    total <- nrow(results)
    output$summary <- renderUI({
      HTML(paste0(
        "<div style='font-size:16px; margin-top:20px'>",
        "<span style='color:", ifelse(correct == total, "green", "red"), "'>",
        "Score: ", correct, "/", total, "</span>",
        "</div>"
      ))
    })
  })
}

shinyApp(ui, server)
