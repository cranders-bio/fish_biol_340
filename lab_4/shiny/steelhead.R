library(shiny)
library(ggplot2)

# -------------------------
# User interface
# -------------------------

ui <- fluidPage(
  titlePanel("Length and Reproductive Success"),

  p("Enter data for 18 individuals, then generate the regression plots."),
  p("Length is the explanatory variable; number of offspring is the response."),

  tags$style(HTML("
    .record-row {
      border-bottom: 1px solid #eeeeee;
      padding-top: 5px;
      padding-bottom: 5px;
    }
    .record-header {
      font-weight: bold;
      padding-bottom: 8px;
    }
  ")),

  fluidRow(
    column(1, strong("Record")),
    column(3, strong("Sex")),
    column(4, strong("Number of offspring")),
    column(4, strong("Length"))
  ),

  # Create 18 input rows
  lapply(seq_len(18), function(i) {
    fluidRow(
      class = "record-row",

      column(
        1,
        div(style = "padding-top: 10px;", i)
      ),

      column(
        3,
        selectInput(
          inputId = paste0("sex", i),
          label = NULL,
          choices = c(
            "Select sex" = "",
            "Male (m)" = "m",
            "Female (f)" = "f"
          ),
          selected = ""
        )
      ),

      column(
        4,
        numericInput(
          inputId = paste0("offspring", i),
          label = NULL,
          value = NA,
          min = 0,
          step = 1
        )
      ),

      column(
        4,
        numericInput(
          inputId = paste0("length", i),
          label = NULL,
          value = NA,
          min = 0
        )
      )
    )
  }),

  br(),

  actionButton(
    "generate",
    "Generate plots",
    class = "btn-primary"
  ),

  br(),
  br(),

  textOutput("status"),

  fluidRow(
    column(
      6,
      h3("Males"),
      plotOutput("male_plot", height = "400px")
    ),
    column(
      6,
      h3("Females"),
      plotOutput("female_plot", height = "400px")
    )
  )
)


# -------------------------
# Server
# -------------------------

server <- function(input, output, session) {

  # Collect data only when the button is clicked
  submitted_data <- eventReactive(input$generate, {

    sex <- vapply(
      seq_len(18),
      function(i) input[[paste0("sex", i)]],
      character(1)
    )

    offspring <- vapply(
      seq_len(18),
      function(i) {
        x <- input[[paste0("offspring", i)]]
        if (is.null(x)) NA_real_ else x
      },
      numeric(1)
    )

    length <- vapply(
      seq_len(18),
      function(i) {
        x <- input[[paste0("length", i)]]
        if (is.null(x)) NA_real_ else x
      },
      numeric(1)
    )

    # Check that all 18 records are complete
    validate(
      need(
        all(sex %in% c("m", "f")),
        "Please select a sex for all 18 individuals."
      ),
      need(
        all(is.finite(offspring)) &&
          all(offspring >= 0) &&
          all(offspring == floor(offspring)),
        "Offspring counts must be non-negative integers."
      ),
      need(
        all(is.finite(length)) &&
          all(length >= 0),
        "Please enter a valid, non-negative length for every individual."
      )
    )

    dat <- data.frame(
      sex = factor(sex, levels = c("m", "f")),
      offspring = offspring,
      length = length
    )

    # Each sex needs at least two observations and
    # at least two distinct lengths for regression.
    for (s in c("m", "f")) {
      group <- dat[dat$sex == s, ]

      validate(
        need(
          nrow(group) >= 2,
          paste(
            "Enter at least two individuals for sex",
            s,
            "to fit a regression."
          )
        ),
        need(
          length(unique(group$length)) >= 2,
          paste(
            "The length values for sex",
            s,
            "must not all be identical."
          )
        )
      )
    }

    dat
  })

  # Display status after successful submission
  output$status <- renderText({
    dat <- submitted_data()

    paste(
      "Successfully submitted", nrow(dat),
      "individuals:",
      sum(dat$sex == "m"), "males and",
      sum(dat$sex == "f"), "females."
    )
  })

  # Reusable plotting function
  make_plot <- function(dat, sex_value, sex_label) {
    
    group <- dat[dat$sex == sex_value, ]
    
    # Fit linear regression
    model <- lm(offspring ~ length, data = group)
    
    # Extract regression coefficients and R-squared
    intercept <- coef(model)[1]
    slope <- coef(model)[2]
    r_squared <- summary(model)$r.squared
    
    # Format the regression equation and R-squared
    equation_label <- sprintf(
      "Offspring = %.2f %s %.2f × Length",
      intercept,
      ifelse(slope < 0, "-", "+"),
      abs(slope)
    )
    
    r_squared_label <- sprintf(
      "R² = %.3f",
      r_squared
    )
    
    # Create plot
    ggplot(
      group,
      aes(x = length, y = offspring)
    ) +
      geom_point(
        color = "#2878B5",
        size = 3,
        alpha = 0.8
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "#C0392B",
        fill = "#E6A09A"
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste(
          equation_label,
          r_squared_label,
          sep = "\n"
        ),
        hjust = 1.1,
        vjust = 1.5,
        size = 4
      ) +
      labs(
        title = paste(sex_label, "Length vs. Offspring"),
        x = "Length",
        y = "Number of offspring"
      ) +
      theme_minimal(base_size = 13)
  }

  # Male regression plot
  output$male_plot <- renderPlot({
    dat <- submitted_data()
    make_plot(dat, "m", "Male")
  })

  # Female regression plot
  output$female_plot <- renderPlot({
    dat <- submitted_data()
    make_plot(dat, "f", "Female")
  })
}

shinyApp(ui = ui, server = server)