library(shiny)
library(shinyjs)
library(dplyr)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Population Genetics Validation: Heterozygosity, F-statistics, and χ² Test"),
  sidebarLayout(
    sidebarPanel(
      actionButton("generate", "Generate New Populations"),
      h4("Sample Data (100 individuals per population)"),
      tableOutput("genoTable"),
      hr(),
      h3("Estimating Deviations from HWE:"),
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
      h3("Chi-squared Test for HWE:"),
      fluidRow(
        column(6,
               numericInput("p1_chi", "χ² Value", value = NA, step = 0.1),
               numericInput("p1_df", "Degrees of Freedom", value = NA, min = 1, max = 10, step = 1),
               selectInput("p1_pval", "P-value (lookup)", 
                           choices = c("<0.001","<0.01","<0.05","<0.10",">0.10"))
        ),
        column(6,
               numericInput("p2_chi", "χ² Value", value = NA, step = 0.1),
               numericInput("p2_df", "Degrees of Freedom", value = NA, min = 1, max = 10, step = 1),
               selectInput("p2_pval", "P-value (lookup)", 
                           choices = c("<0.001","<0.01","<0.05","<0.10",">0.10"))
        )
      ),
      h3("Estimating Population Differentiation:"),
      fluidRow(
        column(6, numericInput("user_Hs", "Hs (Within-pop)", 0, step = 0.01)),
        column(6, numericInput("user_Ht", "Ht (Total)", 0, step = 0.01))
      ),
      numericInput("user_FST", "FST", 0, step = 0.01),
      hr(),
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
  rv <- reactiveValues(data=NULL, trueValues=NULL)
  
  observeEvent(input$generate, {
    shinyjs::enable("validate")
    
    generate_pop <- function() {
      p <- runif(1, 0.1, 0.9)
      FIS <- runif(1, -0.3, 0.3)
      Hexp <- 2 * p * (1 - p)
      Hobs <- Hexp * (1 - FIS)
      
      AA <- max(0, round((p^2 + FIS*p*(1-p)) * 100))
      Aa <- max(0, round(Hobs * 100))
      aa <- max(0, 100 - AA - Aa)
      total <- AA + Aa + aa
      
      if (total != 100) {
        diff <- 100 - total
        largest <- which.max(c(AA, Aa, aa))
        if (largest == 1) AA <- AA + diff else if (largest == 2) Aa <- Aa + diff else aa <- aa + diff
      }
      
      total <- AA + Aa + aa
      p_adj <- (2*AA + Aa)/(2*total)
      Hobs_adj <- Aa/total
      Hexp_adj <- 2 * p_adj * (1 - p_adj)
      FIS_adj <- ifelse(Hexp_adj > 0, (Hexp_adj - Hobs_adj)/Hexp_adj, 0)
      
      # Chi-squared test components
      exp_AA <- p_adj^2 * total
      exp_Aa <- 2*p_adj*(1-p_adj)*total
      exp_aa <- (1 - p_adj)^2 * total
      chi_sq <- sum((c(AA, Aa, aa) - c(exp_AA, exp_Aa, exp_aa))^2 / c(exp_AA, exp_Aa, exp_aa))
      
      df_true <- 1  # Always 1 for HWE with 3 genotypes and 2 alleles
      p_val <- pchisq(chi_sq, df = df_true, lower.tail = FALSE)
      label_p <- if (p_val < 0.001) "<0.001"
      else if (p_val < 0.01) "<0.01"
      else if (p_val < 0.05) "<0.05"
      else if (p_val < 0.10) "<0.10"
      else ">0.10"
      
      list(counts=c(AA=AA, Aa=Aa, aa=aa),
           p=p_adj, Hobs=Hobs_adj, Hexp=Hexp_adj, FIS=FIS_adj,
           ChiSq=chi_sq, df=df_true, P_label=label_p)
    }
    
    pop1 <- generate_pop()
    pop2 <- generate_pop()
    
    p_total <- (2*(pop1$counts["AA"]+pop2$counts["AA"]) +
                  (pop1$counts["Aa"]+pop2$counts["Aa"])) / 400
    Ht <- 2*p_total*(1-p_total)
    Hs <- (pop1$Hexp + pop2$Hexp)/2
    FST <- (Ht - Hs)/Ht
    
    rv$data <- list(pop1=pop1, pop2=pop2)
    rv$trueValues <- list(
      p1_Hobs=pop1$Hobs, p1_Hexp=pop1$Hexp, p1_FIS=pop1$FIS,
      p1_ChiSq=pop1$ChiSq, p1_df=pop1$df, p1_P=pop1$P_label,
      p2_Hobs=pop2$Hobs, p2_Hexp=pop2$Hexp, p2_FIS=pop2$FIS,
      p2_ChiSq=pop2$ChiSq, p2_df=pop2$df, p2_P=pop2$P_label,
      Hs=Hs, Ht=Ht, FST=FST
    )
  })
  
  output$genoTable <- renderTable({
    if (!is.null(rv$data))
      data.frame(
        Population=c("1","2"),
        AA=c(rv$data$pop1$counts["AA"], rv$data$pop2$counts["AA"]),
        Aa=c(rv$data$pop1$counts["Aa"], rv$data$pop2$counts["Aa"]),
        aa=c(rv$data$pop1$counts["aa"], rv$data$pop2$counts["aa"]))
  }, digits=0)
  
  observeEvent(input$validate, {
    shinyjs::disable("validate")
    req(rv$trueValues)
    
    user <- list(
      p1_Hobs=input$p1_Hobs, p1_Hexp=input$p1_Hexp, p1_FIS=input$p1_FIS, 
      p1_ChiSq=input$p1_chi, p1_df=input$p1_df, p1_P=input$p1_pval,
      p2_Hobs=input$p2_Hobs, p2_Hexp=input$p2_Hexp, p2_FIS=input$p2_FIS, 
      p2_ChiSq=input$p2_chi, p2_df=input$p2_df, p2_P=input$p2_pval,
      Hs=input$user_Hs, Ht=input$user_Ht, FST=input$user_FST
    )
    
    results <- data.frame(
      Parameter=c("Pop1 Hobs","Pop1 Hexp","Pop1 FIS","Pop1 χ²","Pop1 df","Pop1 P-value",
                  "Pop2 Hobs","Pop2 Hexp","Pop2 FIS","Pop2 χ²","Pop2 df","Pop2 P-value",
                  "Hs","Ht","FST"),
      Your_Value=unlist(user),
      True_Value=unlist(rv$trueValues),
      Match=mapply(function(a,b){
        if (is.numeric(a) && is.numeric(b)) return(abs(a-b)<0.01)
        else return(a==b)
      }, user, rv$trueValues)
    )
    
    output$validationTable <- renderTable({
      results %>% mutate(
        True_Value=ifelse(sapply(True_Value,is.numeric),
                          round(as.numeric(True_Value),3), True_Value),
        Match=ifelse(Match,"✓","✗")
      )
    })
    
    score <- sum(results$Match)
    total <- nrow(results)
    output$summary <- renderUI({
      HTML(paste0("<div style='font-size:16px; margin-top:20px'>",
                  "<span style='color:", ifelse(score==total,"green","red"), "'>",
                  "Score: ", score, "/", total, "</span></div>"))
    })
  })
}

shinyApp(ui, server)
