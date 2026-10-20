# Pairwise relatedness explorer (Queller & Goodnight 1989)
# ---------------------------------------------------------
# Put this file in a folder together with TurkeyPO.csv and run:
#   shiny::runApp("path/to/folder")
# Requires the packages shiny and ggplot2.
#
# Data format: column 1 = individual ID, column 2 = group
# ("fem" = female, "off" = offspring, anything else = male),
# then two columns per locus. (A column named "father" or "sire", if present,
# is ignored.) 0 or blank = missing allele.
# Offspring are listed directly after their mother.

library(shiny)
library(ggplot2)

default_file <- "../data/TurkeyPO.csv"
classes <- c("Unrelated", "Half sibs", "Full sibs", "Parent-offspring")
class_r <- c("Unrelated" = 0, "Half sibs" = 0.25, "Full sibs" = 0.5, "Parent-offspring" = 0.5)

## ---- Reading the genotype file --------------------------------------------

read_genotypes <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                colClasses = "character", strip.white = TRUE)

  # drop a father/sire column if there is one; it is not used
  fc <- which(tolower(trimws(names(d))) %in% c("father", "sire"))
  if (length(fc)) d <- d[, -fc, drop = FALSE]

  if (ncol(d) < 4 || (ncol(d) - 2) %% 2 != 0)
    stop("Expected an ID column, a group column, then two columns per locus.")

  id  <- trimws(d[[1]])
  pop <- trimws(d[[2]])
  if (anyDuplicated(id))
    stop("Duplicate IDs in file: ", paste(unique(id[duplicated(id)]), collapse = ", "))

  raw <- as.matrix(d[, -(1:2), drop = FALSE])
  g <- matrix(suppressWarnings(as.numeric(raw)), nrow = nrow(raw))
  g[is.na(g) | g <= 0] <- NA

  odd  <- seq(1, ncol(g), by = 2)
  A1   <- g[, odd, drop = FALSE]
  A2   <- g[, odd + 1, drop = FALSE]
  loci <- sub("[aA]$", "", names(d)[-(1:2)][odd])

  # a genotype with only one allele scored is treated as missing
  miss <- is.na(A1) | is.na(A2)
  A1[miss] <- NA
  A2[miss] <- NA
  dimnames(A1) <- dimnames(A2) <- list(id, loci)

  # roles and mothers (offspring follow their mother in the file)
  role <- ifelse(tolower(pop) == "off", "offspring",
                 ifelse(tolower(pop) == "fem", "female", "male"))
  mother <- rep(NA_character_, length(id))
  cur <- NA_character_
  for (i in seq_along(id)) {
    if (role[i] == "female") cur <- id[i]
    if (role[i] == "offspring") mother[i] <- cur
  }

  info <- data.frame(ID = id, Group = pop, Role = role, Mother = mother,
                     `Loci typed` = rowSums(!is.na(A1)),
                     stringsAsFactors = FALSE, check.names = FALSE)
  list(info = info, A1 = A1, A2 = A2, loci = loci)
}

## ---- Allele frequencies ----------------------------------------------------

allele_freqs <- function(A1, A2, ref) {
  lapply(seq_len(ncol(A1)), function(l) {
    a <- c(A1[ref, l], A2[ref, l])
    a <- a[!is.na(a)]
    if (!length(a)) return(numeric(0))
    tab <- table(a)
    setNames(as.numeric(tab) / sum(tab), names(tab))
  })
}

# matrix of the reference frequency of each allele each bird carries
# (alleles absent from the reference sample get frequency 0)
freq_matrix <- function(A, freqs) {
  P <- matrix(NA_real_, nrow(A), ncol(A), dimnames = dimnames(A))
  for (l in seq_len(ncol(A))) {
    v <- unname(freqs[[l]][as.character(A[, l])])
    v[is.na(v) & !is.na(A[, l])] <- 0
    P[, l] <- v
  }
  P
}

## ---- Queller & Goodnight estimator -----------------------------------------
# Symmetric pairwise version: x as reference + y as reference, summed over
# loci before dividing (ratio of sums), as in Queller & Goodnight (1989).

qg_pair <- function(x, y, A1, A2, P1, P2) {
  x1 <- A1[x, ]; x2 <- A2[x, ]; y1 <- A1[y, ]; y2 <- A2[y, ]
  ok <- !is.na(x1) & !is.na(y1)
  if (!any(ok)) return(c(r = NA_real_, loci = 0))

  x1 <- x1[ok]; x2 <- x2[ok]; y1 <- y1[ok]; y2 <- y2[ok]
  px1 <- P1[x, ok]; px2 <- P2[x, ok]; py1 <- P1[y, ok]; py2 <- P2[y, ok]

  shared <- 0.5 * ((x1 == y1) + (x1 == y2) + (x2 == y1) + (x2 == y2))
  num <- sum(shared - px1 - px2) + sum(shared - py1 - py2)
  den <- sum(1 + (x1 == x2) - px1 - px2) + sum(1 + (y1 == y2) - py1 - py2)

  c(r = if (den != 0) num / den else NA_real_, loci = sum(ok))
}

# r for a vector of pairs (row names = individual IDs)
qg_pairs <- function(id1, id2, A1, A2, P1, P2) {
  if (!length(id1)) return(matrix(numeric(0), 2, 0, dimnames = list(c("r", "loci"), NULL)))
  vapply(seq_along(id1), function(i) qg_pair(id1[i], id2[i], A1, A2, P1, P2),
         c(r = 0, loci = 0))
}

## ---- Simulated pairs from allele frequencies --------------------------------
# Genotypes are drawn from the reference allele frequencies (Hardy-Weinberg,
# unlinked loci, no missing data) and r is computed with the same estimator
# and the same reference frequencies as the real birds.

simulate_pairs <- function(freqs, n) {
  freqs <- freqs[lengths(freqs) > 0]
  L <- length(freqs)
  draw <- function() {
    m <- vapply(seq_len(L), function(l)
      as.numeric(sample(names(freqs[[l]]), n, replace = TRUE, prob = freqs[[l]])),
      numeric(n))
    matrix(m, n, L)
  }
  founder <- function() list(a = draw(), b = draw())
  gamete <- function(P) {
    pick <- matrix(runif(n * L) < 0.5, n, L)
    ifelse(pick, P$a, P$b)
  }
  child <- function(P, Q) list(a = gamete(P), b = gamete(Q))
  pfreq <- function(A) {
    P <- A
    for (l in seq_len(L)) {
      v <- unname(freqs[[l]][as.character(A[, l])])
      v[is.na(v)] <- 0
      P[, l] <- v
    }
    P
  }
  qg_sim <- function(X, Y) {
    pxa <- pfreq(X$a); pxb <- pfreq(X$b); pya <- pfreq(Y$a); pyb <- pfreq(Y$b)
    s <- 0.5 * ((X$a == Y$a) + (X$a == Y$b) + (X$b == Y$a) + (X$b == Y$b))
    num <- rowSums(2 * s - pxa - pxb - pya - pyb)
    den <- rowSums(2 + (X$a == X$b) + (Y$a == Y$b) - pxa - pxb - pya - pyb)
    num / den
  }

  P <- founder(); Q <- founder(); R <- founder()
  out <- list(
    "Unrelated"        = qg_sim(founder(), founder()),
    "Half sibs"        = qg_sim(child(P, Q), child(P, R)),
    "Full sibs"        = qg_sim(child(P, Q), child(P, Q)),
    "Parent-offspring" = qg_sim(P, child(P, Q))
  )
  data.frame(Source = "Simulated",
             Class = rep(names(out), lengths(out)),
             Detail = rep(names(out), lengths(out)),
             r = unlist(out, use.names = FALSE), stringsAsFactors = FALSE)
}

## ---- Parsing pasted text ---------------------------------------------------
# Blank lines separate families; pairs are only formed within a family.
# Each line may hold IDs (separated by spaces, commas, tabs, semicolons)
# or a whole row copied from the CSV, in which case only its first field
# (the ID) is used.

parse_set <- function(txt, known, n_cols) {
  txt <- gsub("\r", "", txt)
  blocks <- strsplit(txt, "\n[ \t]*\n")[[1]]
  fams <- lapply(blocks, function(bl) {
    lines <- strsplit(bl, "\n")[[1]]
    lines <- lines[nzchar(trimws(lines))]
    ids <- unlist(lapply(lines, function(ln) {
      tok <- strsplit(trimws(ln), "[,;\t ]+")[[1]]
      tok <- tok[nzchar(tok)]
      # a CSV row has an ID followed by a group label (fem, off, coal1, ...)
      is_row <- length(tok) >= 3 && tok[1] %in% known && !(tok[2] %in% known) &&
        is.na(suppressWarnings(as.numeric(tok[2])))
      if (is_row) tok[1] else tok
    }))
    unique(ids)
  })
  fams[lengths(fams) > 0]
}

# density curve scaled to histogram counts
scaled_density <- function(r, n_pairs, binwidth, from, to) {
  r <- r[is.finite(r)]
  if (length(r) < 3) return(NULL)
  dd <- density(r, from = from, to = to, n = 300)
  data.frame(x = dd$x, y = dd$y * n_pairs * binwidth)
}

## ---- UI --------------------------------------------------------------------

set_box <- function(n, label, cmp) {
  wellPanel(
    style = "padding: 10px 12px; margin-bottom: 10px;",
    textInput(paste0("name", n), NULL, value = label),
    textAreaInput(paste0("ids", n), NULL, rows = 5, resize = "vertical",
                  placeholder = "Paste IDs or whole rows from the CSV.\nBlank line between families."),
    selectInput(paste0("cmp", n), "Compare with expected distribution for",
                c("(none)" = "none", classes), selected = cmp)
  )
}

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 15px; }
    .well textarea { font-family: monospace; font-size: 13px; }
    .help-note { color: #555; font-size: 13px; margin-top: -4px; }
    table { font-size: 14px; }
  "))),
  titlePanel("Pairwise relatedness from microsatellites",
             windowTitle = "Pairwise relatedness"),
  tags$p(style = "color:#777; margin-top:-8px;",
         "Version 5"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("file", "Genotype file (optional, defaults to TurkeyPO.csv)",
                accept = ".csv"),
      radioButtons("freq_ref", "Allele frequencies estimated from",
                   c("Adults only (unrelated reference)" = "adults",
                     "All individuals" = "all")),
      uiOutput("loci_ui"),
      tags$p(class = "help-note", style = "margin-top:-10px;",
             actionLink("loci_all", "all"), " · ", actionLink("loci_none", "none")),
      hr(),
      tags$p(class = "help-note",
             "Each box is one comparison group. Pairs are formed only within a family; ",
             "leave a blank line between families to pool several families into one distribution."),
      set_box(1, "Full-sib families", "Full sibs"),
      set_box(2, "Half-sib families", "Half sibs"),
      set_box(3, "Unrelated adults", "Unrelated"),
      actionLink("fill_adults", "Fill the third box with all adults"),
      hr(),
      checkboxInput("show_sim", "Overlay simulated expected distributions", TRUE),
      sliderInput("nsim", "Simulated pairs per relationship", min = 500, max = 10000,
                  value = 2000, step = 500),
      checkboxInput("show_expected", "Show theoretical values as lines (0, 0.25, 0.5)", TRUE),
      sliderInput("binwidth", "Histogram bin width", min = 0.025, max = 0.2,
                  value = 0.05, step = 0.025),
      radioButtons("layout", "Layout", c("Separate panels" = "facet",
                                         "Overlaid" = "overlay"), inline = TRUE)
    ),
    mainPanel(
      width = 8,
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Distributions",
                 uiOutput("messages"),
                 plotOutput("hist", height = "auto"),
                 h4("Summary"),
                 tableOutput("summary")),
        tabPanel("Expected",
                 br(),
                 tags$p("Expected distributions of r for each relationship, from pairs simulated ",
                        "with the allele frequencies of the ticked loci. Dotted lines are the theoretical values."),
                 plotOutput("exp_plot", height = "720px"),
                 h4("Summary"),
                 tableOutput("exp_summary")),
        tabPanel("Pairwise values",
                 br(),
                 downloadButton("dl", "Download CSV"),
                 br(), br(),
                 tableOutput("pairs")),
        tabPanel("Data",
                 br(),
                 tags$p("Offspring are assigned to the female listed above them in the file."),
                 tableOutput("info")),
        tabPanel("Method", br(), withMathJax(uiOutput("method")))
      )
    )
  )
)

## ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {

  # default data file: TurkeyPO.csv in the app folder, or a renamed download
  # such as "TurkeyPO (1).csv" (an answer-key file with fathers is skipped)
  find_default <- function() {
    if (file.exists(default_file)) return(default_file)
    cand <- list.files(".", pattern = "^TurkeyPO.*\\.csv$", ignore.case = TRUE)
    cand <- cand[!grepl("father", cand, ignore.case = TRUE)]
    if (length(cand)) cand[1] else NA_character_
  }

  raw_dat <- reactive({
    if (!is.null(input$file)) {
      path <- input$file$datapath; shown <- input$file$name
    } else {
      path <- find_default(); shown <- path
    }
    validate(need(!is.na(path) && file.exists(path),
                  paste0("No data file found. Put TurkeyPO.csv in the same folder as app.R ",
                         "(R is looking in: ", getwd(), "), or upload it with the button on the left.")))
    d <- tryCatch(read_genotypes(path),
                  error = function(e) validate(need(FALSE, paste("Could not read the file:", conditionMessage(e)))))
    d$file <- shown
    d
  })

  output$status <- renderUI({
    d <- raw_dat()
    tags$p(style = "color:#2b6a3a; margin:6px 0 10px;",
           sprintf("Loaded %s: %d individuals (%d adults, %d offspring), %d loci.",
                   basename(d$file), nrow(d$info), sum(d$info$Role != "offspring"),
                   sum(d$info$Role == "offspring"), length(d$loci)))
  })

  ## loci tick boxes -------------------------------------------------------
  # The selection is kept in loci_sel(); it starts as all loci of the file.
  loci_sel <- reactiveVal(NULL)
  observeEvent(raw_dat(), loci_sel(raw_dat()$loci))

  output$loci_ui <- renderUI({
    l <- raw_dat()$loci
    checkboxGroupInput("loci", "Loci used", choices = l, selected = l, inline = TRUE)
  })

  # user ticks/unticks (ignoreInit skips the NULL before the boxes exist)
  observeEvent(input$loci, {
    loci_sel(if (is.null(input$loci)) character(0) else input$loci)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$loci_all,
               updateCheckboxGroupInput(session, "loci", selected = raw_dat()$loci))
  observeEvent(input$loci_none,
               updateCheckboxGroupInput(session, "loci", selected = character(0)))

  # genotype data restricted to the ticked loci; everything downstream uses this
  dat <- reactive({
    d <- raw_dat()
    chosen <- loci_sel()
    if (is.null(chosen)) chosen <- d$loci
    sel <- intersect(d$loci, chosen)
    validate(need(length(sel) > 0, "Tick at least one locus."))
    d$A1 <- d$A1[, sel, drop = FALSE]
    d$A2 <- d$A2[, sel, drop = FALSE]
    d$loci <- sel
    d$info$`Loci typed` <- rowSums(!is.na(d$A1))
    d
  })

  loci_label <- reactive(sprintf("Based on %d of %d loci", length(dat()$loci), length(raw_dat()$loci)))

  adults <- reactive(raw_dat()$info$ID[raw_dat()$info$Role != "offspring"])

  fill_adults <- function() {
    updateTextAreaInput(session, "ids3", value = paste(adults(), collapse = " "))
  }
  observeEvent(raw_dat(), fill_adults())      # prefill on load / new upload
  observeEvent(input$fill_adults, fill_adults())

  freqs <- reactive({
    d <- dat()
    ref <- if (input$freq_ref == "adults") adults() else d$info$ID
    f <- allele_freqs(d$A1, d$A2, ref)
    list(f = f, P1 = freq_matrix(d$A1, f), P2 = freq_matrix(d$A2, f))
  })

  ## reference distributions ----------------------------------------------

  sim_ref <- reactive({
    fq <- freqs()
    set.seed(1)
    simulate_pairs(fq$f, input$nsim)
  })

  ref_all <- reactive(sim_ref())

  ## user sets ------------------------------------------------------------

  sets_text <- debounce(reactive(lapply(1:3, function(i) list(
    name = input[[paste0("name", i)]],
    txt  = input[[paste0("ids", i)]],
    cmp  = input[[paste0("cmp", i)]]
  ))), 600)

  results <- reactive({
    d <- dat(); fq <- freqs(); st <- sets_text()
    known <- d$info$ID
    n_cols <- 2 + 2 * ncol(d$A1)
    out <- list(); msgs <- character(0); used <- character(0); cmp <- character(0)

    for (i in seq_along(st)) {
      txt <- st[[i]]$txt
      if (is.null(txt) || !nzchar(trimws(txt))) next
      nm <- trimws(if (is.null(st[[i]]$name)) "" else st[[i]]$name)
      if (!nzchar(nm)) nm <- paste("Set", i)
      if (nm %in% used) nm <- paste0(nm, " (", i, ")")

      fams <- parse_set(txt, known, n_cols)
      unknown <- setdiff(unlist(fams), known)
      if (length(unknown))
        msgs <- c(msgs, sprintf("%s: ignored IDs not in the file: %s",
                                nm, paste(unknown, collapse = ", ")))
      fams <- lapply(fams, intersect, known)
      fams <- fams[lengths(fams) >= 2]
      if (!length(fams)) {
        msgs <- c(msgs, sprintf("%s: needs at least two individuals in a family.", nm))
        next
      }
      used <- c(used, nm)
      cmp[nm] <- if (is.null(st[[i]]$cmp)) "none" else st[[i]]$cmp

      rows <- lapply(seq_along(fams), function(f) {
        pr <- combn(fams[[f]], 2)
        rl <- qg_pairs(pr[1, ], pr[2, ], d$A1, d$A2, fq$P1, fq$P2)
        data.frame(Set = nm, Family = f, ID1 = pr[1, ], ID2 = pr[2, ],
                   r = rl["r", ], Loci = rl["loci", ], stringsAsFactors = FALSE)
      })
      out[[length(out) + 1]] <- do.call(rbind, rows)
    }

    pairs <- if (length(out)) do.call(rbind, out) else NULL
    if (!is.null(pairs)) {
      pairs$Set <- factor(pairs$Set, levels = used)
      rownames(pairs) <- NULL
    }
    list(pairs = pairs, msgs = msgs, cmp = cmp)
  })

  summary_tab <- reactive({
    res <- results(); p <- res$pairs
    req(p)
    ref <- ref_all()
    do.call(rbind, lapply(split(p, p$Set, drop = TRUE), function(s) {
      r <- s$r[!is.na(s$r)]
      nm <- as.character(s$Set[1])
      cl <- res$cmp[[nm]]
      sim <- if (cl == "none") numeric(0) else ref$r[ref$Class == cl]
      data.frame(Set = nm,
                 Families = length(unique(s$Family)),
                 Individuals = length(unique(c(s$ID1, s$ID2))),
                 Pairs = length(r),
                 Mean = mean(r),
                 SD = if (length(r) > 1) sd(r) else NA_real_,
                 Median = median(r),
                 Min = min(r), Max = max(r),
                 `Compared with` = if (cl == "none") "" else cl,
                 `Simulated mean` = if (length(sim)) mean(sim) else NA_real_,
                 `Simulated SD` = if (length(sim) > 1) sd(sim) else NA_real_,
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  })

  output$messages <- renderUI({
    m <- results()$msgs
    if (!length(m)) return(NULL)
    div(style = "background:#fff4e5; border-left:4px solid #e69500; padding:8px 12px; margin:10px 0;",
        lapply(m, tags$div))
  })

  n_sets <- reactive({
    p <- results()$pairs
    if (is.null(p)) 1 else nlevels(p$Set)
  })

  output$hist <- renderPlot({
    res <- results(); p <- res$pairs
    validate(need(!is.null(p),
                  "Paste at least two individual IDs into one of the boxes to see a distribution."))
    p <- p[!is.na(p$r), ]
    s <- summary_tab()
    s$Set <- factor(s$Set, levels = levels(p$Set))
    ref <- ref_all()
    bw <- input$binwidth

    xl <- range(c(-0.5, 1, p$r))
    # expected curves, scaled to each set's number of pairs
    curves <- if (!isTRUE(input$show_sim)) NULL else
      do.call(rbind, lapply(levels(p$Set), function(nm) {
        cl <- res$cmp[[nm]]
        if (is.null(cl) || cl == "none") return(NULL)
        cd <- scaled_density(ref$r[ref$Class == cl], sum(p$Set == nm), bw, xl[1], xl[2])
        if (is.null(cd)) return(NULL)
        cd$Set <- nm
        cd$Expected <- paste("Expected (simulated):", cl)
        cd
      }))
    if (!is.null(curves)) curves$Set <- factor(curves$Set, levels = levels(p$Set))

    g <- ggplot(p, aes(x = r))
    if (input$show_expected)
      g <- g + geom_vline(xintercept = c(0, 0.25, 0.5),
                          linetype = "dotted", colour = "grey55")

    if (input$layout == "overlay") {
      g <- g +
        geom_histogram(aes(fill = Set), binwidth = bw, boundary = 0,
                       position = "identity", alpha = 0.45, colour = "white") +
        geom_vline(data = s, aes(xintercept = Mean, colour = Set), linewidth = 1.2)
      if (!is.null(curves))
        g <- g + geom_line(data = curves, aes(x = x, y = y, colour = Set, linetype = Expected),
                           linewidth = 0.9)
    } else {
      g <- g +
        geom_histogram(binwidth = bw, boundary = 0, fill = "#4C78A8", colour = "white") +
        geom_vline(data = s, aes(xintercept = Mean), colour = "#C8364A", linewidth = 1.2) +
        geom_text(data = s, aes(x = Mean, y = Inf,
                                label = sprintf("  mean = %.3f  (n = %d pairs)", Mean, Pairs)),
                  vjust = 1.5, hjust = 0, size = 4.3, colour = "#C8364A") +
        facet_wrap(~Set, ncol = 1, scales = "free_y")
      if (!is.null(curves))
        g <- g + geom_line(data = curves, aes(x = x, y = y, linetype = Expected),
                           colour = "grey15", linewidth = 0.9)
    }

    g + coord_cartesian(xlim = xl) +
      scale_linetype_manual(values = c("solid", "longdash", "dotdash", "twodash",
                                       "dashed", "dotted", "11", "41")) +
      labs(x = "Pairwise relatedness r (Queller & Goodnight)", y = "Number of pairs",
           subtitle = loci_label(), fill = NULL, colour = NULL, linetype = NULL) +
      theme_minimal(base_size = 15) +
      theme(legend.position = "top", legend.direction = "vertical",
            strip.text = element_text(face = "bold", hjust = 0, size = 15),
            panel.grid.minor = element_blank())
  }, height = function() {
    if (isTRUE(input$layout == "overlay")) 480 else 160 + 240 * n_sets()
  })

  output$summary <- renderTable(summary_tab(), digits = 3, na = "")

  ## expected tab ---------------------------------------------------------

  output$exp_plot <- renderPlot({
    ref <- ref_all()
    ref$Class <- factor(ref$Class, levels = classes)
    means <- aggregate(r ~ Class, data = ref, FUN = mean)
    theo <- data.frame(Class = factor(names(class_r), levels = classes), r = class_r)

    ggplot(ref, aes(x = r)) +
      geom_histogram(binwidth = input$binwidth, boundary = 0,
                     fill = "#9AA5B1", colour = "white") +
      geom_vline(data = theo, aes(xintercept = r), linetype = "dotted", colour = "grey25") +
      geom_vline(data = means, aes(xintercept = r), colour = "#C8364A", linewidth = 1.1) +
      geom_text(data = means, aes(x = r, y = Inf, label = sprintf("  mean = %.3f", r)),
                vjust = 1.5, hjust = 0, size = 4.3, colour = "#C8364A") +
      facet_wrap(~Class, ncol = 1) +
      coord_cartesian(xlim = range(c(-0.75, 1.1, quantile(ref$r, c(0.001, 0.999))))) +
      labs(x = "Pairwise relatedness r (Queller & Goodnight)", y = "Number of simulated pairs",
           subtitle = loci_label()) +
      theme_minimal(base_size = 15) +
      theme(strip.text = element_text(face = "bold", hjust = 0, size = 15),
            panel.grid.minor = element_blank())
  })

  output$exp_summary <- renderTable({
    ref <- ref_all()
    do.call(rbind, lapply(classes, function(cl) {
      v <- ref$r[ref$Class == cl]
      data.frame(Relationship = cl, `Theoretical r` = class_r[[cl]], Pairs = length(v),
                 Mean = mean(v), SD = sd(v),
                 `2.5%` = unname(quantile(v, 0.025)), `97.5%` = unname(quantile(v, 0.975)),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }, digits = 3)

  ## other tabs -----------------------------------------------------------

  output$pairs <- renderTable({
    p <- results()$pairs
    validate(need(!is.null(p), "No pairs yet."))
    p
  }, digits = 3, na = "")

  output$dl <- downloadHandler(
    filename = function() "pairwise_relatedness.csv",
    content = function(file) write.csv(results()$pairs, file, row.names = FALSE)
  )

  output$info <- renderTable({
    info <- dat()$info
    info
  }, na = "")

  output$method <- renderUI(tagList(
    p("For each pair of individuals x and y, relatedness is estimated with the Queller & Goodnight (1989) estimator, ",
      "summing numerator and denominator over all loci typed in both birds and averaging the two directions ",
      "(x as reference and y as reference):"),
    p("$$r_{xy} = \\frac{\\sum_{l}\\left[\\tfrac12(\\delta_{ac}+\\delta_{ad}+\\delta_{bc}+\\delta_{bd}) - p_a - p_b\\right] + \\sum_{l}\\left[\\tfrac12(\\delta_{ac}+\\delta_{ad}+\\delta_{bc}+\\delta_{bd}) - p_c - p_d\\right]}{\\sum_{l}\\left[1 + \\delta_{ab} - p_a - p_b\\right] + \\sum_{l}\\left[1 + \\delta_{cd} - p_c - p_d\\right]}$$"),
    p("Here x has alleles a, b and y has alleles c, d at locus l; \\(\\delta_{ac}\\) = 1 if a and c are the same allele, otherwise 0; ",
      "\\(p_a\\) is the population frequency of allele a."),
    h4("Expected distributions"),
    tags$ul(
      tags$li("Pairs of each relationship are simulated by drawing founder genotypes from the allele frequencies ",
              "of the ticked loci (Hardy–Weinberg, independent loci, no missing data or scoring errors) and passing ",
              "alleles to offspring by Mendelian inheritance. Their spread shows how much r varies just because ",
              "only a limited number of loci are typed. Untick loci to see the distributions widen."),
      tags$li("Theoretical values: unrelated 0, half sibs 0.25, full sibs and parent–offspring 0.5.")
    ),
    h4("Notes"),
    tags$ul(
      tags$li("Allele frequencies come from the adults by default. Including offspring over-represents ",
              "their parents' alleles. Alleles never seen in the reference sample get frequency 0."),
      tags$li("Because the frequencies are estimated from the same birds, real pairs are shifted slightly below the ",
              "theoretical values (unrelated adults average about −1/(N−1) for N reference individuals). ",
              "Simulated birds are not part of the reference sample, so they do not show this shift."),
      tags$li("A locus is skipped for a pair if either bird is missing a genotype (0 in the file).")
    )
  ))
}

shinyApp(ui, server)
