# Genotype-calling cutoffs and what they do to population genetics
# -----------------------------------------------------------------------
# A teaching app for introductory genetics.
#
# Two simulated fish populations are genotyped with GTseq. Every fish is
# genotyped twice (two independent libraries), so we can measure how often
# replicate genotypes disagree. Students change the sequencing error rate,
# read depth, genotype-calling cutoffs and filters, and see:
#   * missing data vs. genotyping errors (replicate discordance)
#   * how both distort a PCA and the estimate of FST between populations
#
# Run it:
#   install.packages(c("shiny", "ggplot2"))   # once
#   shiny::runApp("gtseq_cutoffs_app")        # folder that holds this file
# -----------------------------------------------------------------------

library(shiny)
library(ggplot2)

call_cols <- c("A/A" = "#0072B2", "A/B" = "#009E73", "B/B" = "#E69F00",
               "No call" = "grey70")
pop_cols  <- c("Population 1" = "#0072B2", "Population 2" = "#D55E00")

# ---------- Simulate two populations ----------
# Allele frequencies in the two populations are drawn around a shared
# ancestral frequency (Balding-Nichols model) with the chosen true FST.
# Everything that does not depend on the error rate is drawn first, so
# moving the error slider changes the reads, not the fish.
simulate_data <- function(n_per, n_loc, mean_depth, err, imbalance, fst, seed,
                          sd_fish = 0.8, sd_loc = 0.5, af_range = "wide",
                          mono_frac = 0.15) {
  set.seed(seed)
  n   <- 2 * n_per
  pop <- rep(c("Population 1", "Population 2"), each = n_per)
  # Ancestral allele frequencies: either a wide range that includes rare
  # alleles (U-shaped), or common alleles only
  p <- if (af_range == "common") runif(n_loc, 0.1, 0.9) else
         pmin(pmax(rbeta(n_loc, 0.5, 0.5), 0.02), 0.98)
  a   <- p * (1 - fst) / fst
  b   <- (1 - p) * (1 - fst) / fst
  p1  <- rbeta(n_loc, a, b)
  p2  <- rbeta(n_loc, a, b)
  # Some panel loci do not vary in these populations (monomorphic): every
  # fish is homozygous, so any minor-allele reads are sequencing errors
  mono  <- runif(n_loc) < mono_frac
  fixed <- ifelse(runif(n_loc) < 0.5, 0, 1)
  p1[mono] <- fixed[mono]; p2[mono] <- fixed[mono]
  freq <- rbind(matrix(p1, n_per, n_loc, byrow = TRUE),
                matrix(p2, n_per, n_loc, byrow = TRUE))
  g <- matrix(rbinom(n * n_loc, 2, freq), nrow = n)   # true genotype = no. of A alleles

  # Sequencing success varies among fish (DNA quality) and among loci
  # (amplification). Log-normal with mean 1; larger SD = more uneven.
  dna_quality <- exp(rnorm(n, 0, sd_fish) - sd_fish^2 / 2)
  locus_eff   <- exp(rnorm(n_loc, 0, sd_loc) - sd_loc^2 / 2)

  d <- expand.grid(rep = 1:2, ind = 1:n, loc = 1:n_loc)
  d$true  <- g[cbind(d$ind, d$loc)]
  mu      <- mean_depth * dna_quality[d$ind] * locus_eff[d$loc]
  d$depth <- rnbinom(nrow(d), size = 3, mu = mu)

  if (imbalance > 0) {
    k <- 0.25 / imbalance^2 - 1
    q <- rbeta(nrow(d), k * 0.5, k * 0.5)
  } else {
    q <- rep(0.5, nrow(d))
  }
  frac_A   <- ifelse(d$true == 2, 1, ifelse(d$true == 0, 0, q))
  p_read_A <- frac_A * (1 - err) + (1 - frac_A) * err   # an error flips a read to the other allele
  d$A <- rbinom(nrow(d), d$depth, p_read_A)
  d$B <- d$depth - d$A
  list(reads = d, true = g, pop = pop, n = n, n_loc = n_loc)
}

# ---------- Call genotypes from allele counts ----------
# ratio = A reads / B reads
#   A/A  if ratio >= hom          (A has at least `hom` times the reads of B)
#   B/B  if ratio <= 1/hom
#   A/B  if 1/het <= ratio <= het (neither allele has more than `het` times the other)
#   no call otherwise, or with fewer than min_depth reads
call_genotypes <- function(A, B, min_depth, hom, het) {
  het   <- min(het, hom)
  n     <- A + B
  ratio <- ifelse(B == 0, Inf, A / pmax(B, 1e-12))
  out <- rep(NA_integer_, length(A))
  out[ratio >= hom]     <- 2L
  out[ratio <= 1 / hom] <- 0L
  out[ratio <= het & ratio >= 1 / het] <- 1L
  out[n < min_depth | n == 0] <- NA_integer_
  out
}

geno_label <- function(x) ifelse(is.na(x), "No call", c("B/B", "A/B", "A/A")[x + 1])

# Plain-language version of the rules for one example read depth
describe_rules <- function(n, min_depth, hom, het) {
  if (n < min_depth)
    return(paste0("With ", n, " reads: no call (fewer than the minimum of ", min_depth, ")."))
  B <- 0:n
  lab <- geno_label(call_genotypes(n - B, B, min_depth, hom, het))
  runs <- rle(lab)
  ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  parts <- mapply(function(l, s, e) {
    rng <- if (s == e) B[s] else paste0(B[s], "-", B[e])
    paste0(rng, " B reads = ", l)
  }, runs$values, starts, ends)
  paste0("Example with ", n, " reads (A + B): ", paste(parts, collapse = "; "), ".")
}

# ---------- Summaries ----------
summarise_calls <- function(d, calls) {
  r1 <- calls[d$rep == 1]; r2 <- calls[d$rep == 2]
  both <- !is.na(r1) & !is.na(r2)
  called <- !is.na(calls)
  list(missing   = mean(is.na(calls)),
       discord   = if (any(both)) mean(r1[both] != r2[both]) else NA,
       n_pairs   = sum(both),
       n_discord = sum(r1[both] != r2[both]),
       true_err  = if (any(called)) mean(calls[called] != d$true[called]) else NA)
}

# Genotype matrix (fish x loci) from replicate 1
geno_matrix <- function(sim, calls) {
  matrix(calls[sim$reads$rep == 1], nrow = sim$n)
}

# Minor allele frequency per locus from the called genotypes
called_maf <- function(M) {
  q <- colMeans(M, na.rm = TRUE) / 2
  q[is.na(q)] <- 0
  pmin(q, 1 - q)
}

# Remove loci (too much missing data or minor allele frequency too low),
# then fish with too much missing data
apply_filters <- function(M, loc_min, fish_min, maf_min = 0) {
  keep_loc  <- colMeans(!is.na(M)) >= loc_min & called_maf(M) >= maf_min
  M2 <- M[, keep_loc, drop = FALSE]
  keep_fish <- if (ncol(M2) > 0) rowMeans(!is.na(M2)) >= fish_min else rep(FALSE, nrow(M))
  list(M = M2[keep_fish, , drop = FALSE], keep_fish = keep_fish, keep_loc = keep_loc)
}

# Hudson's FST (ratio of averages; Bhatia et al. 2013)
fst_hudson <- function(M, pop) {
  num <- 0; den <- 0
  for (j in seq_len(ncol(M))) {
    x1 <- M[pop == "Population 1", j]; x1 <- x1[!is.na(x1)]
    x2 <- M[pop == "Population 2", j]; x2 <- x2[!is.na(x2)]
    if (length(x1) < 2 || length(x2) < 2) next
    q1 <- mean(x1) / 2; q2 <- mean(x2) / 2
    n1 <- 2 * length(x1); n2 <- 2 * length(x2)
    num <- num + (q1 - q2)^2 - q1 * (1 - q1) / (n1 - 1) - q2 * (1 - q2) / (n2 - 1)
    den <- den + q1 * (1 - q2) + q2 * (1 - q1)
  }
  if (den == 0) NA else num / den
}

# Ho and FIS within populations (averaged over the two populations)
ho_fis <- function(M, pop) {
  res <- sapply(unique(pop), function(pp) {
    X  <- M[pop == pp, , drop = FALSE]
    ok <- !is.na(X)
    Ho <- mean(X[ok] == 1)
    q  <- colMeans(X, na.rm = TRUE) / 2
    He <- mean(2 * q * (1 - q), na.rm = TRUE)
    c(Ho = Ho, He = He)
  })
  Ho <- mean(res["Ho", ]); He <- mean(res["He", ])
  c(Ho = Ho, FIS = 1 - Ho / He)
}

# PCA with missing genotypes replaced by the locus mean
pca_scores <- function(M) {
  q    <- colMeans(M, na.rm = TRUE) / 2
  keep <- !is.na(q) & q > 0 & q < 1
  if (sum(keep) < 2 || nrow(M) < 3) return(NULL)
  X <- M[, keep, drop = FALSE]; q <- q[keep]
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- 2 * q[j]
  Z  <- sweep(sweep(X, 2, 2 * q), 2, sqrt(2 * q * (1 - q)), "/")
  pc <- prcomp(Z, center = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  list(x = pc$x[, 1:2, drop = FALSE], var = ve[1:2], n_loci = ncol(X), n_fish = nrow(X))
}

pct <- function(x, digits = 1) {
  if (is.na(x)) "-" else paste0(formatC(100 * x, format = "f", digits = digits), "%")
}
num3 <- function(x) if (is.na(x)) "-" else formatC(x, format = "f", digits = 3)

stat_box <- function(label, value, sub = NULL, color = "#222") {
  div(style = "flex:1 1 150px; background:#f4f6f7; border-radius:8px; padding:10px 14px;",
      div(style = "font-size:13px; color:#555;", label),
      div(style = paste0("font-size:26px; font-weight:600; color:", color, ";"), value),
      if (!is.null(sub)) div(style = "font-size:12px; color:#666;", sub))
}

# Default settings. The calling rules match the GTseq genotyper
# (Campbell et al. 2015): at least 10 reads, allele ratio of 10 or more
# for a homozygote, 5 or less (0.2 to 5) for a heterozygote.
defaults <- list(err = 0.01, depth = 60, fold_fish = 3, fold_loc = 2.75, mono_frac = 0.15,
                 min_depth = 10, hom = 10, het = 5,
                 loc_min = 0, maf_min = 0, fish_min = 0, fst_true = 0.05, n_per = 20,
                 n_loc = 80, imb = 0.08)

base_theme <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# ======================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 15px; }
    .well { background: #fafbfb; }
    .irs--shiny .irs-bar { background: #0072B2; border-color: #0072B2; }
    h4 { margin-top: 18px; }
    .rule { font-size: 13px; color: #444; background: #eef3f6; border-radius: 6px;
            padding: 8px 10px; margin: -6px 0 12px; }
  "))),
  titlePanel("Why genotype-calling cutoffs matter"),
  sidebarLayout(
    sidebarPanel(width = 4,
      actionButton("reset", "Reset to GTseq genotyper defaults",
                   style = "width:100%; margin-bottom:8px;"),
      h4("The sequencing run"),
      sliderInput("err", "Sequencing error (chance a read shows the wrong allele)",
                  min = 0, max = 0.10, value = 0.01, step = 0.005),
      sliderInput("depth", "Average reads per fish per locus",
                  min = 5, max = 150, value = 60, step = 5),

      h4("Genotype-calling rules"),
      sliderInput("min_depth", "Minimum total reads: fewer reads than this gives no call",
                  min = 0, max = 60, value = 10, step = 1),
      sliderInput("hom", "Homozygote cutoff: call A/A (or B/B) when one allele has at least this many times the reads of the other",
                  min = 1.5, max = 40, value = 10, step = 0.5),
      sliderInput("het", "Heterozygote cutoff: call A/B when neither allele has more than this many times the reads of the other",
                  min = 1, max = 10, value = 5, step = 0.5),
      div(class = "rule", textOutput("rule_text")),
      uiOutput("cutoff_warning"),

      h4("Filtering after calling"),
      sliderInput("fish_min", "Remove fish genotyped at fewer than this share of loci",
                  min = 0, max = 1, value = 0, step = 0.05),
      sliderInput("loc_min", "Remove loci genotyped in fewer than this share of fish",
                  min = 0, max = 1, value = 0, step = 0.05),
      sliderInput("maf_min", "Remove loci with a minor allele frequency below",
                  min = 0, max = 0.25, value = 0, step = 0.01),

      tags$details(
        tags$summary("More settings"),
        sliderInput("fold_fish", "Unevenness in reads among fish (fold spread; 1 = every fish the same)",
                    min = 1, max = 6, value = 3, step = 0.25),
        div(class = "rule", textOutput("fish_text")),
        sliderInput("fold_loc", "Unevenness in reads among loci (fold spread; 1 = every locus the same)",
                    min = 1, max = 6, value = 2.75, step = 0.25),
        sliderInput("mono_frac", "Share of loci that do not vary (monomorphic) in these populations",
                    min = 0, max = 0.5, value = 0.15, step = 0.05),
        selectInput("af_range", "Allele frequencies across loci",
                    c("Wide range, including rare alleles" = "wide",
                      "Common alleles only (0.1 to 0.9)" = "common")),
        sliderInput("fst_true", "True FST between the populations",
                    min = 0.01, max = 0.3, value = 0.05, step = 0.01),
        sliderInput("n_per", "Fish per population", min = 10, max = 50, value = 20, step = 5),
        sliderInput("n_loc", "Number of loci", min = 20, max = 200, value = 80, step = 10),
        sliderInput("imb", "Allele-balance noise in heterozygotes (SD of A fraction)",
                    min = 0, max = 0.2, value = 0.08, step = 0.01),
        actionButton("reseed", "Simulate a new data set")
      )
    ),
    mainPanel(width = 8,
      tabsetPanel(
        tabPanel("Genotype quality",
          br(),
          uiOutput("boxes"),
          h4("Every replicate genotype"),
          helpText("Each point is one fish at one locus in one replicate. Black rings",
                   "mark replicate pairs that disagree. Dashed lines are the homozygote",
                   "cutoff, solid lines the heterozygote cutoff, and the dotted curve at",
                   "the lower left the minimum read count."),
          plotOutput("scatter", height = "540px"),
          h4("The trade-off"),
          fluidRow(
            column(6, selectInput("sweep", "Vary this setting:",
                                  c("Homozygote cutoff" = "hom",
                                    "Heterozygote cutoff" = "het",
                                    "Minimum total reads" = "min_depth"))),
            column(6, br(), helpText("Other settings stay as set on the left.",
                                     "The dashed line marks your current value."))
          ),
          plotOutput("sweep_plot", height = "380px")
        ),
        tabPanel("Two populations",
          br(),
          p("The fish come from two populations that differ a little in allele",
            "frequencies. Compare what the true genotypes show with what the",
            "called genotypes show after your cutoffs and filters."),
          uiOutput("pop_boxes"),
          h4("Which fish and loci get filtered?"),
          helpText("Share of genotypes called for each fish and locus, before filtering.",
                   "Dashed lines are your filter thresholds; crosses are removed. Red loci",
                   "do not really vary: if they show a minor allele frequency above 0,",
                   "sequencing errors made them look like SNPs (false SNPs)."),
          fluidRow(column(6, plotOutput("callrate_fish", height = "300px")),
                   column(6, plotOutput("callrate_loci", height = "300px"))),
          h4("Genetic diversity"),
          helpText("Expected heterozygosity (He), averaged over loci, in each population:",
                   "true genotypes of all fish and loci vs. called genotypes after your",
                   "cutoffs and filters. Watch the called bars as you move the minor allele",
                   "frequency filter."),
          plotOutput("he_bar", height = "260px"),
          h4("PCA and FST"),
          helpText("Each point is a fish. Panel titles give the number of fish and loci used.",
                   "In the called-genotype PCA, larger points have more missing genotypes;",
                   "missing genotypes are filled with the locus average, which pulls those",
                   "fish toward the middle. The bars in the middle compare FST from the true",
                   "and the called genotypes."),
          fluidRow(column(5, plotOutput("pca_true", height = "380px")),
                   column(2, plotOutput("fst_bar", height = "380px")),
                   column(5, plotOutput("pca_called", height = "380px")))
        )
      )
    )
  )
)

# ======================================================================
server <- function(input, output, session) {

  seed <- reactiveVal(2026)
  observeEvent(input$reseed, seed(sample.int(1e6, 1)))

  observeEvent(input$reset, {
    for (id in names(defaults)) updateSliderInput(session, id, value = defaults[[id]])
    updateSelectInput(session, "af_range", selected = "wide")
    seed(2026)
  })

  sim <- reactive({
    simulate_data(input$n_per, input$n_loc, input$depth, input$err,
                  input$imb, input$fst_true, seed(),
                  sd_fish = log(input$fold_fish), sd_loc = log(input$fold_loc), af_range = input$af_range,
                  mono_frac = input$mono_frac)
  })

  calls <- reactive({
    d <- sim()$reads
    call_genotypes(d$A, d$B, input$min_depth, input$hom, input$het)
  })

  output$fish_text <- renderText({
    sdf <- log(input$fold_fish)
    q   <- function(z) round(input$depth * exp(z * sdf - sdf^2 / 2))
    if (sdf == 0) return(paste("Every fish averages", input$depth, "reads per locus."))
    paste0("Two-thirds of fish average ", q(-1), " to ", q(1), " reads per locus; ",
           "the poorest 1 in 10 fish average fewer than ", q(qnorm(0.1)), ".")
  })

  output$rule_text <- renderText(describe_rules(40, input$min_depth, input$hom, input$het))

  output$cutoff_warning <- renderUI({
    if (input$het > input$hom)
      div(style = "color:#b3541e; font-size:13px;",
          "The heterozygote cutoff is larger than the homozygote cutoff,",
          "so it is treated as equal to the homozygote cutoff.")
  })

  # ===================== Genotype quality =====================
  output$boxes <- renderUI({
    s <- summarise_calls(sim()$reads, calls())
    div(style = "display:flex; flex-wrap:wrap; gap:10px;",
      stat_box("Missing data", pct(s$missing), "genotypes with no call", "#555"),
      stat_box("Replicate discordance", pct(s$discord, 2),
               paste(s$n_discord, "of", s$n_pairs, "pairs disagree"), "#C0392B"),
      stat_box("Error vs. the truth", pct(s$true_err, 2),
               "only knowable in a simulation", "#7F4FA0"))
  })

  output$scatter <- renderPlot({
    d <- sim()$reads; cl <- calls()
    d$call <- factor(geno_label(cl), names(call_cols))
    r1 <- cl[d$rep == 1]; r2 <- cl[d$rep == 2]
    d$discordant <- rep(!is.na(r1) & !is.na(r2) & r1 != r2, each = 2)

    mx  <- max(c(d$A, d$B, 10))
    xs  <- exp(seq(0, log(mx + 1), length.out = 200)) - 1
    het <- min(input$het, input$hom)
    ln  <- function(m, type, grp) data.frame(x = xs, y = xs / m, type = type, grp = grp)
    lines <- rbind(ln(input$hom, "Homozygote cutoff", 1), ln(1 / input$hom, "Homozygote cutoff", 2),
                   ln(het, "Heterozygote cutoff", 3), ln(1 / het, "Heterozygote cutoff", 4))
    lines <- lines[lines$y <= mx, ]
    md <- input$min_depth
    depth_line <- data.frame(x = seq(0, md, length.out = 100))
    depth_line$y <- md - depth_line$x

    ggplot(d, aes(A, B)) +
      geom_line(data = lines, aes(x, y, group = grp, linetype = type), colour = "grey35") +
      { if (md > 0) geom_line(data = depth_line, aes(x, y), colour = "grey35", linetype = "dotted") } +
      geom_point(aes(colour = call), alpha = 0.7, size = 2) +
      geom_point(data = d[d$discordant, ], shape = 21, size = 4, stroke = 1.1,
                 colour = "black", fill = NA) +
      scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
      scale_colour_manual(values = call_cols, drop = FALSE) +
      scale_linetype_manual(values = c("Homozygote cutoff" = "dashed",
                                       "Heterozygote cutoff" = "solid")) +
      coord_equal() +
      labs(x = "Reads with allele A", y = "Reads with allele B", colour = "Call", linetype = NULL) +
      base_theme + theme(legend.position = "right")
  })

  output$sweep_plot <- renderPlot({
    d <- sim()$reads
    vals <- switch(input$sweep, hom = seq(1.5, 40, by = 0.5),
                   het = seq(1, 10, by = 0.25), min_depth = 0:60)
    label <- switch(input$sweep, hom = "Homozygote cutoff",
                    het = "Heterozygote cutoff", min_depth = "Minimum total reads")
    res <- do.call(rbind, lapply(vals, function(v) {
      args <- list(min_depth = input$min_depth, hom = input$hom, het = input$het)
      args[[input$sweep]] <- v
      s <- summarise_calls(d, call_genotypes(d$A, d$B, args$min_depth, args$hom, args$het))
      data.frame(value = v,
                 panel = c("Missing data (%)", "Errors (%)", "Errors (%)"),
                 measure = c("Missing data", "Replicate discordance", "Error vs. the truth"),
                 y = 100 * c(s$missing, s$discord, s$true_err))
    }))
    res$panel   <- factor(res$panel, c("Missing data (%)", "Errors (%)"))
    res$measure <- factor(res$measure, c("Missing data", "Replicate discordance", "Error vs. the truth"))
    ggplot(res, aes(value, y, colour = measure, linetype = measure)) +
      geom_line(linewidth = 1.1, na.rm = TRUE) +
      geom_vline(xintercept = input[[input$sweep]], linetype = "dashed", colour = "grey30") +
      facet_wrap(~panel, ncol = 1, scales = "free_y") +
      scale_colour_manual(values = c("Missing data" = "grey40",
                                     "Replicate discordance" = "#C0392B",
                                     "Error vs. the truth" = "#7F4FA0")) +
      scale_linetype_manual(values = c("solid", "solid", "22")) +
      expand_limits(y = 0) +
      labs(x = label, y = NULL, colour = NULL, linetype = NULL) +
      base_theme + theme(strip.text = element_text(face = "bold", hjust = 0))
  })

  # ===================== Two populations =====================
  filtered <- reactive({
    s <- sim()
    M <- geno_matrix(s, calls())
    f <- apply_filters(M, input$loc_min, input$fish_min, input$maf_min)
    f$pop <- s$pop[f$keep_fish]
    f$missing_fish <- rowMeans(is.na(M[, f$keep_loc, drop = FALSE]))[f$keep_fish]
    f
  })

  callrates <- reactive({
    s <- sim(); M <- geno_matrix(s, calls()); f <- filtered()
    list(fish = data.frame(rate = rowMeans(!is.na(M)), pop = s$pop,
                           status = ifelse(f$keep_fish, "Kept", "Removed")),
         loci = data.frame(rate = colMeans(!is.na(M)), maf = called_maf(M),
                           type = ifelse(colMeans(s$true) %in% c(0, 2),
                                         "Monomorphic: no real SNP", "Real SNP"),
                           status = ifelse(f$keep_loc, "Kept", "Removed")))
  })

  output$callrate_fish <- renderPlot({
    d <- callrates()$fish
    d <- d[order(d$rate), ]; d$rank <- seq_len(nrow(d))
    ggplot(d, aes(rank, rate, colour = pop, shape = status)) +
      geom_hline(yintercept = input$fish_min, linetype = "dashed", colour = "grey40") +
      geom_point(size = 3) +
      scale_colour_manual(values = pop_cols) +
      scale_shape_manual(values = c(Kept = 16, Removed = 4), drop = FALSE) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(title = "Fish", x = "Fish, sorted by call rate", y = "Loci genotyped",
           colour = NULL, shape = NULL) + base_theme
  })

  output$callrate_loci <- renderPlot({
    d <- callrates()$loci
    ggplot(d, aes(maf, rate, shape = status, colour = type)) +
      geom_hline(yintercept = input$loc_min, linetype = "dashed", colour = "grey40") +
      geom_point(size = 3) +
      scale_colour_manual(values = c("Real SNP" = "#555555",
                                     "Monomorphic: no real SNP" = "#C0392B")) +
      scale_shape_manual(values = c(Kept = 16, Removed = 4), drop = FALSE) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      geom_vline(xintercept = input$maf_min, linetype = "dashed", colour = "grey40") +
      labs(title = "Loci", x = "Minor allele frequency (called genotypes)",
           y = "Fish genotyped", shape = NULL, colour = NULL) + base_theme +
      guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  })

  bar_fill <- c("True genotypes" = "#A6A6A6", "Called genotypes" = "#0072B2")

  output$he_bar <- renderPlot({
    s <- sim(); f <- filtered()
    he <- function(X) {
      q <- colMeans(X, na.rm = TRUE) / 2
      mean(2 * q * (1 - q), na.rm = TRUE)
    }
    res <- do.call(rbind, lapply(names(pop_cols), function(pp) {
      rows <- f$pop == pp
      data.frame(pop = pp, type = c("True genotypes", "Called genotypes"),
                 He = c(he(s$true[s$pop == pp, , drop = FALSE]),
                        if (sum(rows) >= 2) he(f$M[rows, , drop = FALSE]) else NA))
    }))
    res$type <- factor(res$type, names(bar_fill))
    ggplot(res, aes(pop, He, fill = type)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7, na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(He), "", sprintf("%.3f", He))),
                position = position_dodge(width = 0.75), vjust = -0.4, size = 4.5) +
      scale_fill_manual(values = bar_fill) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      expand_limits(y = 0.5) +
      labs(x = NULL, y = "Average He per locus", fill = NULL) + base_theme
  })

  output$pop_boxes <- renderUI({
    s <- sim(); f <- filtered()
    t_hf <- ho_fis(s$true, s$pop)
    o_hf <- if (nrow(f$M) > 0) ho_fis(f$M, f$pop) else c(Ho = NA, FIS = NA)
    M_all <- geno_matrix(s, calls())
    true_mono  <- colMeans(s$true) %in% c(0, 2)
    looks_var  <- called_maf(M_all) > 0
    false_snps <- sum(true_mono & looks_var & f$keep_loc)
    div(style = "display:flex; flex-wrap:wrap; gap:10px;",
      stat_box("Observed heterozygosity", num3(o_hf["Ho"]),
               paste("true:", num3(t_hf["Ho"]))),
      stat_box("FIS within populations", num3(o_hf["FIS"]),
               paste("true:", num3(t_hf["FIS"]))),
      stat_box("False SNPs left after filtering", false_snps,
               paste0("of ", sum(true_mono), " monomorphic loci; ",
                      sum(true_mono & looks_var), " look variable before filtering"),
               if (false_snps > 0) "#C0392B" else "#222"),
      stat_box("Data kept after filtering",
               paste0(sum(f$keep_fish), " fish"),
               paste0(sum(f$keep_loc), " of ", length(f$keep_loc), " loci")))
  })

  pcas <- reactive({
    s <- sim(); f <- filtered()
    pt <- pca_scores(s$true)
    po <- if (nrow(f$M) >= 3) pca_scores(f$M) else NULL
    if (!is.null(pt) && !is.null(po)) {
      # PCA axes have arbitrary signs: flip the called PCA to match the true one
      tx <- pt$x[f$keep_fish, , drop = FALSE]
      for (k in 1:2) if (isTRUE(cor(po$x[, k], tx[, k]) < 0)) po$x[, k] <- -po$x[, k]
    }
    list(true = pt, called = po)
  })

  pca_panel <- function(pc, pop, missing, title) {
    d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], pop = pop, missing = missing)
    ggplot(d, aes(PC1, PC2, colour = pop, size = missing)) +
      geom_point(alpha = 0.8) +
      scale_colour_manual(values = pop_cols) +
      scale_size_continuous(range = c(2, 7), limits = c(0, 1), guide = "none") +
      labs(title = title,
           subtitle = sprintf("%d fish, %d loci", pc$n_fish, pc$n_loci),
           x = sprintf("PC1 (%.0f%%)", pc$var[1]), y = sprintf("PC2 (%.0f%%)", pc$var[2]),
           colour = NULL) +
      base_theme + theme(legend.position = "bottom")
  }

  output$pca_true <- renderPlot({
    pt <- pcas()$true
    validate(need(!is.null(pt), "Not enough variable loci for a PCA."))
    pca_panel(pt, sim()$pop, 0, "True genotypes")
  })

  output$pca_called <- renderPlot({
    po <- pcas()$called; f <- filtered()
    validate(need(!is.null(po), "Too few fish or loci left for a PCA."))
    pca_panel(po, f$pop, f$missing_fish, "Called genotypes")
  })

  output$fst_bar <- renderPlot({
    s <- sim(); f <- filtered()
    d <- data.frame(type = factor(names(bar_fill), names(bar_fill)),
                    FST = c(fst_hudson(s$true, s$pop),
                            if (nrow(f$M) > 3) fst_hudson(f$M, f$pop) else NA))
    ggplot(d, aes(type, FST, fill = type)) +
      geom_col(width = 0.7, na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(FST), "", sprintf("%.3f", FST))),
                vjust = -0.4, size = 4) +
      scale_fill_manual(values = bar_fill, guide = "none") +
      scale_x_discrete(labels = c("True", "Called")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "FST", x = NULL, y = NULL) + base_theme
  })
}

shinyApp(ui, server)
