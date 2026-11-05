# ===========================================================================
# Herring GT-seq: a guided population-genetics analysis
#
# Part A: filter loci and learn the methods on three reference populations
#         (SquaxinPass, Cherry14, ElliottBay).
# Part B: apply the locus set to Oro Bay and work out what is going on there.
#
# Run with:  shiny::runApp("HerringApp")
# Instructor mode (all tabs open + notes): add ?instructor=1 to the app URL.
# Only the 'shiny' package is required.
# ===========================================================================

library(shiny)
# ---------------------------------------------------------------------------
# Population-genetic helper functions for the Herring GT-seq teaching app.
# Base R only. Genotypes are coded as dosage of allele 1 (0/1/2, NA = missing),
# with individuals in rows and loci in columns.
# ---------------------------------------------------------------------------

# Convert a character vector of genotypes ("AG", "TT", NA) to dosage.
to_dosage <- function(g) {
  called <- g[!is.na(g) & nchar(g) == 2]
  alleles <- sort(unique(unlist(strsplit(called, ""))))
  if (length(alleles) == 0) {
    return(list(d = rep(NA_real_, length(g)), alleles = ""))
  }
  a1 <- alleles[1]
  ok <- !is.na(g) & nchar(g) == 2
  d <- rep(NA_real_, length(g))
  d[ok] <- (substr(g[ok], 1, 1) == a1) + (substr(g[ok], 2, 2) == a1)
  list(d = d, alleles = paste(alleles, collapse = "/"))
}

# Per-locus summary for a set of individuals.
locus_stats <- function(X) {
  n  <- colSums(!is.na(X))
  p  <- colMeans(X, na.rm = TRUE) / 2
  Ho <- colMeans(X == 1, na.rm = TRUE)
  He <- 2 * p * (1 - p) * ifelse(n > 0, (2 * n) / (2 * n - 1), NA)
  maf <- pmin(p, 1 - p)
  nAA <- colSums(X == 2, na.rm = TRUE)
  nAB <- colSums(X == 1, na.rm = TRUE)
  nBB <- colSums(X == 0, na.rm = TRUE)
  data.frame(n = n, p = p, maf = maf, Ho = Ho, He = He,
             nAA = nAA, nAB = nAB, nBB = nBB,
             FIS = ifelse(is.finite(He) & He > 0, 1 - Ho / He, NA))
}

# Multilocus FIS (ratio of sums) with a bootstrap CI over loci.
fis_multilocus <- function(X, B = 200) {
  s <- locus_stats(X)
  ok <- is.finite(s$He) & s$He > 0 & s$n >= 5
  Ho <- s$Ho[ok]; He <- s$He[ok]
  if (length(Ho) < 2) return(c(FIS = NA, lo = NA, hi = NA, loci = length(Ho)))
  est <- 1 - sum(Ho) / sum(He)
  bs <- replicate(B, {
    i <- sample.int(length(Ho), replace = TRUE)
    1 - sum(Ho[i]) / sum(He[i])
  })
  q <- stats::quantile(bs, c(0.025, 0.975), names = FALSE)
  c(FIS = est, lo = q[1], hi = q[2], loci = length(Ho))
}

# Exact test for Hardy-Weinberg equilibrium (Wigginton et al. 2005).
hwe_exact <- function(nAB, nAA, nBB) {
  n <- nAA + nAB + nBB
  if (n == 0) return(NA_real_)
  nA <- 2 * nAA + nAB
  nB <- 2 * nBB + nAB
  rare <- min(nA, nB)
  if (rare == 0) return(NA_real_)
  probs <- numeric(rare + 1)
  mid <- floor(rare * (2 * n - rare) / (2 * n))
  if ((rare %% 2) != (mid %% 2)) mid <- mid + 1
  probs[mid + 1] <- 1
  h <- mid; hr <- (rare - mid) / 2; hc <- n - h - hr
  while (h >= 2) {
    probs[h - 1] <- probs[h + 1] * h * (h - 1) / (4 * (hr + 1) * (hc + 1))
    h <- h - 2; hr <- hr + 1; hc <- hc + 1
  }
  h <- mid; hr <- (rare - mid) / 2; hc <- n - h - hr
  while (h <= rare - 2) {
    probs[h + 3] <- probs[h + 1] * 4 * hr * hc / ((h + 2) * (h + 1))
    h <- h + 2; hr <- hr - 1; hc <- hc - 1
  }
  probs <- probs / sum(probs)
  p <- sum(probs[probs <= probs[nAB + 1] * (1 + 1e-7)])
  min(1, p)
}

# HWE p-values for every locus in X.
hwe_table <- function(X) {
  s <- locus_stats(X)
  vapply(seq_len(ncol(X)), function(j) hwe_exact(s$nAB[j], s$nAA[j], s$nBB[j]),
         numeric(1))
}

# Hudson FST (Bhatia et al. 2013), ratio of averages over loci.
fst_hudson <- function(A, B) {
  p1 <- colMeans(A, na.rm = TRUE) / 2; p2 <- colMeans(B, na.rm = TRUE) / 2
  n1 <- 2 * colSums(!is.na(A));        n2 <- 2 * colSums(!is.na(B))
  num <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  ok <- is.finite(num) & is.finite(den) & den > 0 & n1 > 2 & n2 > 2
  sum(num[ok]) / sum(den[ok])
}

# FST with a permutation test (individuals shuffled between the two groups).
fst_perm <- function(A, B, nperm = 999) {
  obs <- fst_hudson(A, B)
  X <- rbind(A, B); na <- nrow(A)
  perm <- replicate(nperm, {
    i <- sample.int(nrow(X))
    fst_hudson(X[i[seq_len(na)], , drop = FALSE], X[i[-seq_len(na)], , drop = FALSE])
  })
  c(FST = obs, p = (sum(perm >= obs) + 1) / (nperm + 1))
}

# Pairwise FST matrix for a list of genotype matrices.
fst_matrix <- function(groups, nperm = 999) {
  g <- names(groups); k <- length(g)
  Fm <- matrix(NA, k, k, dimnames = list(g, g)); Pm <- Fm
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    r <- fst_perm(groups[[i]], groups[[j]], nperm)
    Fm[i, j] <- Fm[j, i] <- r["FST"]; Pm[i, j] <- Pm[j, i] <- r["p"]
  }
  list(FST = Fm, p = Pm, nperm = nperm)
}

# PCA on standardized genotypes (missing set to the mean).
pca_geno <- function(X) {
  p <- colMeans(X, na.rm = TRUE) / 2
  keep <- is.finite(p) & p > 0 & p < 1
  X <- X[, keep, drop = FALSE]; p <- p[keep]
  Z <- sweep(X, 2, 2 * p)
  Z <- sweep(Z, 2, sqrt(2 * p * (1 - p)), "/")
  Z[is.na(Z)] <- 0
  pr <- stats::prcomp(Z, center = FALSE)
  v <- pr$sdev^2 / sum(pr$sdev^2)
  list(x = pr$x[, 1:2, drop = FALSE], var = v[1:2], nloci = ncol(X))
}

# Individual assignment to baseline groups (as in GeneClass / supervised STRUCTURE).
# Allele frequencies get +0.5 pseudocounts; equal priors; missing loci are skipped.
assign_freq <- function(X) (colSums(X, na.rm = TRUE) + 0.5) / (2 * colSums(!is.na(X)) + 1)
assign_loglik <- function(X, p) {
  n <- nrow(X)
  l2 <- matrix(2 * log(p), n, length(p), byrow = TRUE)
  l1 <- matrix(log(2 * p * (1 - p)), n, length(p), byrow = TRUE)
  l0 <- matrix(2 * log(1 - p), n, length(p), byrow = TRUE)
  rowSums(ifelse(X == 2, l2, ifelse(X == 1, l1, l0)), na.rm = TRUE)
}
to_post <- function(LL) { LL <- LL - apply(LL, 1, max); E <- exp(LL); E / rowSums(E) }
# baseline: named list of genotype matrices; query: genotype matrix of fish to assign.
# Baseline fish themselves are assigned leave-one-out.
assign_fish <- function(baseline, query = NULL) {
  K <- names(baseline)
  P <- lapply(baseline, assign_freq)
  out <- list()
  for (k in K) {
    X <- baseline[[k]]
    cnt <- colSums(X, na.rm = TRUE); nn <- colSums(!is.na(X))
    LL <- sapply(K, function(j) assign_loglik(X, P[[j]]))
    LL <- matrix(LL, nrow = nrow(X))
    for (i in seq_len(nrow(X))) {             # leave fish i out of its own baseline
      xi <- X[i, ]; ok <- !is.na(xi)
      p_loo <- (cnt - ifelse(ok, xi, 0) + 0.5) / (2 * (nn - ok) + 1)
      LL[i, match(k, K)] <- assign_loglik(X[i, , drop = FALSE], p_loo)
    }
    out[[k]] <- to_post(LL)
  }
  q <- NULL
  if (!is.null(query) && nrow(query) > 0)
    q <- to_post(matrix(sapply(K, function(j) assign_loglik(query, P[[j]])), nrow = nrow(query)))
  list(base = out, query = q, K = K)
}

# Key results for a set of reference fish and loci (used by the comparison panels).
key_results <- function(G, pops, inds, loci, ref_pops) {
  out <- c(Fish = length(inds), Loci = length(loci))
  for (pp in ref_pops) {
    X <- G[inds[pops[inds] == pp], loci, drop = FALSE]
    out[paste("FIS", pp)] <- fis_multilocus(X, B = 2)["FIS"]
  }
  for (pr in utils::combn(ref_pops, 2, simplify = FALSE)) {
    out[paste("FST", pr[1], "\u2013", pr[2])] <-
      fst_hudson(G[inds[pops[inds] == pr[1]], loci, drop = FALSE], G[inds[pops[inds] == pr[2]], loci, drop = FALSE])
  }
  out
}
compare_table <- function(a, b, name_a, name_b) {
  df <- data.frame(Result = names(a), A = a, B = b, check.names = FALSE)
  df$A <- ifelse(df$Result %in% c("Fish", "Loci"), sprintf("%.0f", df$A), sprintf("%.3f", df$A))
  df$B <- ifelse(df$Result %in% c("Fish", "Loci"), sprintf("%.0f", df$B), sprintf("%.3f", df$B))
  names(df) <- c("Result", name_a, name_b)
  df
}

# Small HTML table for the report.
html_table <- function(df, digits = 3) {
  if (is.null(df) || nrow(df) == 0) return("<p><em>None.</em></p>")
  fmt <- function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "fg") else as.character(v)
  head <- paste0("<tr>", paste0("<th>", htmltools::htmlEscape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(df[i, , drop = FALSE], fmt, "")
    paste0("<tr>", paste0("<td>", htmltools::htmlEscape(cells), "</td>", collapse = ""), "</tr>")
  }, "")
  paste0("<table>", head, paste(body, collapse = ""), "</table>")
}

# ---- Locate the data folder ---------------------------------------------------
DATA_DIR <- "data"
if (!file.exists(file.path(DATA_DIR, "genotypes.csv"))) {
  stop("Cannot find the 'data' folder. Start the app with ",
       "shiny::runApp(\"<path to>/HerringApp\"), or open app.R in RStudio and click 'Run App'.",
       call. = FALSE)
}

# ---- Data ------------------------------------------------------------------
geno_raw  <- read.csv(file.path(DATA_DIR, "genotypes.csv"), check.names = FALSE,
                      colClasses = "character", na.strings = c("", "NA"))
samples   <- read.csv(file.path(DATA_DIR, "samples.csv"), stringsAsFactors = FALSE,
                      na.strings = c("", "NA"))
loci_info <- read.csv(file.path(DATA_DIR, "loci.csv"), stringsAsFactors = FALSE)

samples <- samples[match(geno_raw$Sample, samples$Sample), ]
samples$CaptureDate <- as.Date(samples$CaptureDate)
samples$Month <- as.integer(format(samples$CaptureDate, "%m"))

LOCI <- loci_info$Locus
conv <- lapply(LOCI, function(l) to_dosage(geno_raw[[l]]))
G <- sapply(conv, function(x) x$d)
colnames(G) <- LOCI
rownames(G) <- geno_raw$Sample
loci_info$Alleles <- vapply(conv, function(x) x$alleles, character(1))
CHROM <- setNames(loci_info$Chromosome, LOCI)
POS   <- setNames(loci_info$Position, LOCI)

REF_POPS <- c("SquaxinPass", "Cherry14", "ElliottBay")
POP      <- samples$Population
IS_REF   <- samples$Source == "Reference"
IS_ORO   <- samples$Source == "OroBay"
CALLRATE <- rowMeans(!is.na(G))

STAGE_ORDER <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII")
STAGES <- STAGE_ORDER[STAGE_ORDER %in% samples$Stage[IS_ORO]]

# Okabe-Ito colours; references = filled circles, Oro Bay groups = other symbols.
COLS <- c(SquaxinPass = "#0072B2", Cherry14 = "#E69F00", ElliottBay = "#009E73",
          Immature = "#CC79A7", `Early spawners` = "#D55E00", `Late spawners` = "#000000")
PCH  <- c(SquaxinPass = 16, Cherry14 = 16, ElliottBay = 16,
          Immature = 17, `Early spawners` = 15, `Late spawners` = 18)
ORO_GROUPS <- c("Immature", "Early spawners", "Late spawners")

# Genome-wide x axis: chromosomes laid end to end (positions in Mb).
CHR_LEVELS <- unique(loci_info$Chromosome)
CHR_MAX <- tapply(loci_info$Position, factor(loci_info$Chromosome, levels = CHR_LEVELS), max)
CHR_OFF <- setNames(c(0, cumsum(as.numeric(CHR_MAX) + 3e6))[seq_along(CHR_LEVELS)], CHR_LEVELS)
GPOS <- setNames((loci_info$Position + CHR_OFF[loci_info$Chromosome]) / 1e6, LOCI)
CHR_LAB <- setNames(sub("^NC_0451?", "", sub("\\.1$", "", CHR_LEVELS)), CHR_LEVELS)
STAGE_COLS <- setNames(hcl.colors(length(STAGES), "Dark 3"), STAGES)

STEPS <- c("A1", "A2", "A3", "B1", "B2", "B3", "B4")

ref_idx <- function(inds, pop) inds[POP[inds] == pop]
fmt3 <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "f"))

# ---- UI helpers ------------------------------------------------------------
decision_box <- function(id, prompt) {
  wellPanel(class = "decision",
    h4("Your decision"),
    textAreaInput(paste0(id, "_why"), prompt, width = "100%", rows = 3),
    actionButton(paste0(id, "_commit"), "Commit decision and continue",
                 class = "btn-primary"))
}
note <- function(id) uiOutput(paste0("note_", id))
intro <- function(...) div(class = "intro", ...)

# ---- UI --------------------------------------------------------------------
ui <- navbarPage(
  title = "Herring GT-seq", id = "nav",
  header = tags$head(tags$style(HTML("
    body { font-size: 17px; }
    .navbar, .navbar-nav > li > a { font-size: 17px; }
    h2 { font-size: 30px; } h3 { font-size: 25px; } h4 { font-size: 20px; }
    .help-block, .form-group label, .control-label { font-size: 16px; }
    .intro { max-width: 1000px; margin-bottom: 12px; }
    .decision { background: #f4f8fb; border-color: #c9dbe8; }
    .instructor { background: #fff8e6; border: 1px solid #f0d68a; padding: 10px 14px;
                  border-radius: 6px; margin-top: 12px; }
    .hint { color: #8a5a00; }
    table.table { font-size: 16px; }
  "))),

  # ---- Start ----
  tabPanel("Start", value = "start",
    fluidPage(intro(
      h2("Herring GT-seq: a guided analysis"),
      p("You will analyse GT-seq genotypes from Pacific herring in two parts."),
      tags$ol(
        tags$li(strong("Part A: the reference populations."),
                " SquaxinPass, Cherry14 and ElliottBay were genotyped well. You will use them to decide which loci to keep, ",
                "test Hardy\u2013Weinberg equilibrium, and estimate FST and a PCA."),
        tags$li(strong("Part B: Oro Bay."),
                " You will apply your locus set to fish from Oro Bay, decide how to split them into groups, ",
                "and work out what is going on there.")),
      p("Each step ends with a decision. Write a short justification and commit it to unlock the next step. ",
        "Your choices carry forward, so a decision made early shows up in later results. ",
        "You can go back and change a decision, but later steps will then have to be redone."),
      p("At the end you download a report with every decision, justification and key result."),
      textInput("student", "Your name", width = "300px"),
      actionButton("start_go", "Begin Part A", class = "btn-primary")
    ))),

  # ---- A1 ----
  tabPanel("A1 Data QC", value = "A1",
    sidebarLayout(
      sidebarPanel(width = 4,
        helpText("Filters are applied to the 96 reference fish only."),
        sliderInput("a1_ind", "Minimum call rate per individual", 0, 1, 0.5, 0.05),
        sliderInput("a1_loc", "Minimum call rate per locus", 0, 1, 0.5, 0.05),
        sliderInput("a1_maf", "Minimum minor allele frequency (MAF)", 0, 0.2, 0, 0.01),
        radioButtons("a1_mafbasis", "Apply the MAF filter to",
          c("All reference fish pooled" = "pooled",
            "Each population separately (all must pass)" = "each")),
        checkboxInput("a1_mono", "Remove loci that are monomorphic in the pooled reference fish", FALSE),
        tableOutput("a1_counts"),
        decision_box("a1", "Why did you choose these thresholds?")),
      mainPanel(width = 8,
        intro(p("Look at missing data per fish and per locus, and at how common the rarer allele is. ",
                "Rare alleles carry little information and can make some estimators fail; ",
                "strict filters, on the other hand, throw data away.")),
        plotOutput("a1_plots", height = "740px"),
        checkboxInput("a1_cmp", "Compare key results with no filtering at all", FALSE),
        conditionalPanel("input.a1_cmp", tableOutput("a1_cmp_table"),
          helpText("FIS and FST computed with all reference fish and all loci (left) and with your filters (right).")),
        note("A1")))),

  # ---- A2 ----
  tabPanel("A2 HWE & FIS", value = "A2",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("a2_corr", "Multiple-testing correction",
          c("None (p < 0.05)" = "none", "Bonferroni" = "bonferroni",
            "False discovery rate (Benjamini\u2013Hochberg)" = "BH"),
          selected = character(0)),
        radioButtons("a2_rule", "Remove a locus if it is out of HWE in",
          c("All three populations" = "3", "At least two populations" = "2",
            "Any population" = "1", "Never (keep all loci)" = "0"),
          selected = character(0)),
        tableOutput("a2_counts"),
        decision_box("a2", "Why did you choose this correction and this rule?")),
      mainPanel(width = 8,
        intro(p("Each locus is tested for Hardy\u2013Weinberg equilibrium (HWE) in each population with an exact test. ",
                "FIS measures the heterozygote deficit (positive) or excess (negative).")),
        tabsetPanel(
          tabPanel("Populations", tableOutput("a2_summary"),
                   helpText("Multilocus FIS = 1 \u2212 \u03a3Ho / \u03a3He, with He corrected for sample size (2n / (2n \u2212 1)) and a 95% bootstrap interval over loci. ",
                            "GenAlEx reports F = mean over loci of 1 \u2212 Ho / He with uncorrected He; with ~30 fish that value is about 0.015 lower.")),
          tabPanel("Loci out of HWE", tableOutput("a2_loci")),
          tabPanel("Compare", p("Key results before and after removing loci that are out of HWE."),
                   tableOutput("a2_cmp_table")),
          tabPanel("Per-locus FIS", plotOutput("a2_plot", height = "580px"))),
        note("A2")))),

  # ---- A4 ----
  tabPanel("A3 FST & PCA", value = "A3",
    sidebarLayout(
      sidebarPanel(width = 4,
        helpText("FST is Hudson's estimator; p-values come from 999 permutations of fish between populations."),
        decision_box("a4", "Describe how the reference populations relate to each other and why you trust your locus set.")),
      mainPanel(width = 8,
        uiOutput("a4_h_final"), tableOutput("a4_fst"),
        helpText("Upper triangle: FST. Lower triangle: permutation p-values (999 permutations). ",
                 "'p \u2264 0.001' means that none of the permuted data sets reached the observed FST."),
        plotOutput("a4_pca", height = "600px"),
        note("A3")))),

  # ---- B1 ----
  tabPanel("B1 Split Oro Bay", value = "B1",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Spawners (ripe fish)"),
        helpText("Fish that pass the GSI and stage criteria are late spawners if caught on or after the date, early spawners if caught before it."),
        sliderInput("b1_gsi", "Minimum gonadosomatic index (GSI)", 0, 0.25, 0, 0.01),
        sliderInput("b1_date", "Late spawners: caught on or after",
                    min = min(samples$CaptureDate[IS_ORO], na.rm = TRUE),
                    max = max(samples$CaptureDate[IS_ORO], na.rm = TRUE),
                    value = min(samples$CaptureDate[IS_ORO], na.rm = TRUE), step = 1, timeFormat = "%d %b"),
        checkboxGroupInput("b1_late_st", "Spawning stages allowed", STAGES, selected = STAGES, inline = TRUE),
        h4("Immature fish"),
        checkboxGroupInput("b1_imm_st", "Spawning stages counted as immature", STAGES, inline = TRUE),
        checkboxInput("b1_imm_gsi", "Also require GSI below the late-spawner cutoff", FALSE),
        h4("Genotype quality"),
        sliderInput("b1_cr", "Minimum call rate per fish (share of loci genotyped)", 0, 1, 0, 0.05),
        sliderInput("b1_loc", "Minimum call rate per locus (share of your Oro Bay fish genotyped)", 0, 1, 0, 0.05),
        tableOutput("b1_counts"),
        decision_box("b1", "Explain how you defined the two groups and which fish you left out.")),
      mainPanel(width = 8,
        intro(p("The Oro Bay fish come with metadata: capture date, spawning stage and GSI (gonad weight / body weight). ",
                "Use the plots to decide which fish are late spawners and which are immature. ",
                "Fish that fit neither definition are left out.")),
        plotOutput("b1_plots", height = "1100px"),
        h4("Fish by spawning stage and month"), tableOutput("b1_table"),
        note("B1")))),

  # ---- B2 ----
  tabPanel("B2 FIS in Oro Bay", value = "B2",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("b2_imm", "What best explains FIS in the immature fish?",
          c("A Wahlund effect (mixed stocks)" = "wahlund", "Genotyping artefact (allelic dropout)" = "dropout",
            "Inbreeding" = "inbreeding", "Both a Wahlund effect and dropout" = "both",
            "No real deficit" = "none"), selected = character(0)),
        radioButtons("b2_early", "What best explains FIS in the early spawners?",
          c("A Wahlund effect (mixed stocks)" = "wahlund", "Genotyping artefact (allelic dropout)" = "dropout",
            "Inbreeding" = "inbreeding", "No real deficit" = "none"), selected = character(0)),
        radioButtons("b2_late", "What best explains FIS in the late spawners?",
          c("A Wahlund effect (mixed stocks)" = "wahlund", "Genotyping artefact (allelic dropout)" = "dropout",
            "Inbreeding" = "inbreeding", "No real deficit" = "none"), selected = character(0)),
        decision_box("b2", "Justify your answers with the numbers.")),
      mainPanel(width = 8,
        intro(p("Your locus set from Part A is applied to the Oro Bay groups you defined in B1. ",
                "The plots recalculate FIS while the call-rate filters on fish (left) and on loci (right) get stricter. ",
                "If missing data biased the estimates, FIS would change as the filter tightens; ",
                "a Wahlund effect is a property of the group and should stay.")),
        tableOutput("b2_table"),
        plotOutput("b2_plot", height = "520px"),
        note("B2")))),

  # ---- B3 ----
  tabPanel("B3 FST, PCA & assignment", value = "B3",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("b3_late", "Are the late spawners a distinct population?",
          c("Yes" = "yes", "No, they belong to one of the references" = "no", "Cannot tell" = "unsure"),
          selected = character(0)),
        radioButtons("b3_early", "What are the early spawners?",
          c("A mixture of stocks" = "mixture", "A single distinct population" = "single",
            "Mostly one of the reference populations" = "ref", "Cannot tell" = "unsure"),
          selected = character(0)),
        radioButtons("b3_assign", "With the late spawners as a baseline, the early spawners are mostly",
          c("Late-spawner fish, plus some fish from other stocks" = "late_plus", "SquaxinPass or ElliottBay fish" = "south",
            "Cherry Point fish" = "cherry", "Impossible to say" = "unsure"), selected = character(0)),
        radioButtons("b3_imm", "What are the immature fish?",
          c("A mixture of stocks" = "mixture", "A single distinct population" = "single",
            "Mostly one of the reference populations" = "ref", "Cannot tell" = "unsure"),
          selected = character(0)),
        decision_box("b3", "Justify your answers with FST, the PCA and the assignment test.")),
      mainPanel(width = 8,
        tabsetPanel(
          tabPanel("FST & PCA",
            uiOutput("b3_h"),
            tableOutput("b3_fst"),
            helpText("p-values come from 999 permutations of fish between the two groups. ",
                     "'p \u2264 0.001' means that none of the permuted data sets reached the observed FST."),
            plotOutput("b3_pca", height = "600px")),
          tabPanel("Assignment",
            p("Each fish is assigned to the baseline groups by how likely its genotype is under each group's allele frequencies ",
              "(the same idea as STRUCTURE run with known groups, or GeneClass). One bar per fish; colours show assignment probabilities. ",
              "Baseline fish are assigned with themselves left out, which shows how well the baselines can be told apart."),
            fluidRow(
              column(7, checkboxGroupInput("b3_base", "Baseline groups", c(REF_POPS, ORO_GROUPS),
                                           selected = REF_POPS, inline = TRUE)),
              column(5, checkboxInput("b3_other", "Also assign Oro Bay fish left out of the groups in B1 (spent fish, intermediate GSI)", FALSE))),
            plotOutput("b3_assign_plot", height = "470px"),
            h4("Fish assigned to each baseline (highest probability)"),
            tableOutput("b3_assign_table"))),
        note("B3")))),

  # ---- B4 ----
  tabPanel("B4 Synthesis", value = "B4",
    fluidPage(intro(
      h3("What is going on in Oro Bay?"),
      textAreaInput("b4_text", "Write your conclusion (a short paragraph).", width = "100%", rows = 8),
      downloadButton("report", "Download report (HTML)"),
      hr(), h4("Your decision log"), tableOutput("log_table"),
      note("B4"))))
)

# ---- Server ----------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(log = data.frame(Step = character(), Time = character(),
                                        Choices = character(), Justification = character(),
                                        stringsAsFactors = FALSE))

  instructor <- reactive({
    q <- parseQueryString(session$clientData$url_search)
    !is.null(q$instructor)
  })

  # Gating: only A1 is open at the start (all tabs in instructor mode).
  for (s in STEPS[-1]) hideTab("nav", s)
  observeEvent(instructor(), {
    if (instructor()) for (s in STEPS) showTab("nav", s)
  }, once = TRUE)

  observeEvent(input$start_go, updateNavbarPage(session, "nav", selected = "A1"))

  why_ok <- function(id) {
    txt <- input[[paste0(id, "_why")]]
    if (is.null(txt) || nchar(trimws(txt)) < 20) {
      showNotification("Please write a short justification (at least one sentence) before committing.",
                       type = "warning")
      return(FALSE)
    }
    TRUE
  }
  add_log <- function(step, choices, why) {
    rv$log <- rbind(rv$log, data.frame(Step = step, Time = format(Sys.time(), "%H:%M"),
                                       Choices = choices, Justification = why,
                                       stringsAsFactors = FALSE))
  }
  # Committing a step clears every later step, then opens the next one.
  advance <- function(step) {
    i <- match(step, STEPS)
    later <- STEPS[seq_along(STEPS) > i]
    for (s in later) { rv[[s]] <- NULL; rv[[paste0("res_", s)]] <- NULL }
    if ("B3" %in% later) rv$res_B3b <- NULL
    if (!instructor()) for (s in later[-1]) hideTab("nav", s)
    if (length(later)) {
      showTab("nav", later[1])
      updateNavbarPage(session, "nav", selected = later[1])
    }
    showNotification(paste("Step", step, "committed."), type = "message")
  }
  need_step <- function(step) validate(need(!is.null(rv[[step]]), paste("Commit step", step, "first.")))

  # ======================= A1: data QC ===================================
  a1_res <- reactive({
    idx  <- which(IS_REF)
    cr_i <- rowMeans(!is.na(G[idx, , drop = FALSE]))
    inds <- idx[cr_i >= input$a1_ind]
    Xk   <- G[inds, , drop = FALSE]
    cr_l <- colMeans(!is.na(Xk))
    maf  <- locus_stats(Xk)$maf
    maf_pop <- sapply(REF_POPS, function(pp) locus_stats(G[ref_idx(inds, pp), , drop = FALSE])$maf)
    mono <- !is.finite(maf) | maf == 0
    maf_ok <- if (input$a1_mafbasis == "pooled") maf >= input$a1_maf else
      apply(maf_pop, 1, function(v) all(is.finite(v) & v >= input$a1_maf))
    keep <- (cr_l >= input$a1_loc) & (maf_ok %in% TRUE) & !(input$a1_mono & mono)
    keep[is.na(keep)] <- FALSE
    list(inds = inds, loci = LOCI[keep], cr_i = cr_i, cr_l = cr_l, maf = maf,
         n_mono = sum(mono), n_ref = length(idx))
  })

  output$a1_counts <- renderTable({
    r <- a1_res()
    data.frame(Item = c("Reference fish kept", "Loci kept", "Monomorphic loci (pooled)"),
               Value = c(sprintf("%d of %d", length(r$inds), r$n_ref),
                         sprintf("%d of %d", length(r$loci), length(LOCI)),
                         as.character(r$n_mono)))
  }, colnames = FALSE)

  output$a1_cmp_table <- renderTable({
    r <- a1_res()
    compare_table(key_results(G, POP, which(IS_REF), LOCI, REF_POPS),
                  key_results(G, POP, r$inds, r$loci, REF_POPS), "No filtering", "Your filters")
  })
  output$a1_plots <- renderPlot(res = 110, {
    r <- a1_res()
    op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    hist(r$cr_i, breaks = seq(0, 1, 0.025), col = "grey80", border = "white",
         main = "Call rate per reference fish", xlab = "Share of loci genotyped", ylab = "Fish")
    abline(v = input$a1_ind, col = "red", lwd = 2)
    hist(r$cr_l, breaks = seq(0, 1, 0.025), col = "grey80", border = "white",
         main = "Call rate per locus (fish that passed)", xlab = "Share of fish genotyped", ylab = "Loci")
    abline(v = input$a1_loc, col = "red", lwd = 2)
    m <- r$maf; m[!is.finite(m)] <- NA
    hist(m, breaks = seq(0, 0.5, 0.02), col = "grey80", border = "white",
         main = "Minor allele frequency (pooled)", xlab = "MAF", ylab = "Loci")
    abline(v = input$a1_maf, col = "red", lwd = 2)
    plot(r$cr_l, m, pch = 19, col = adjustcolor("grey30", 0.6), xlim = c(0, 1), ylim = c(0, 0.5),
         xlab = "Locus call rate", ylab = "MAF", main = "Loci: call rate vs MAF")
    abline(v = input$a1_loc, h = input$a1_maf, col = "red", lty = 2)
  })

  observeEvent(input$a1_commit, {
    if (!why_ok("a1")) return()
    r <- a1_res()
    if (length(r$loci) < 10) {
      showNotification("Fewer than 10 loci are left. Relax your filters.", type = "error"); return()
    }
    rv$A1 <- list(inds = r$inds, loci = r$loci)
    add_log("A1 Data QC",
      sprintf("Fish call rate \u2265 %.2f; locus call rate \u2265 %.2f; MAF \u2265 %.2f (%s); monomorphic removed: %s. Kept %d fish and %d loci.",
              input$a1_ind, input$a1_loc, input$a1_maf, input$a1_mafbasis,
              ifelse(input$a1_mono, "yes", "no"), length(r$inds), length(r$loci)),
      input$a1_why)
    advance("A1")
  })

  # ======================= A2: HWE and FIS ===============================
  a2_base <- reactive({
    need_step("A1")
    inds <- rv$A1$inds; loci <- rv$A1$loci
    pv  <- sapply(REF_POPS, function(pp) hwe_table(G[ref_idx(inds, pp), loci, drop = FALSE]))
    fis <- sapply(REF_POPS, function(pp) locus_stats(G[ref_idx(inds, pp), loci, drop = FALSE])$FIS)
    rownames(pv) <- rownames(fis) <- loci
    set.seed(42)
    summ <- do.call(rbind, lapply(REF_POPS, function(pp) {
      X <- G[ref_idx(inds, pp), loci, drop = FALSE]
      s <- locus_stats(X); f <- fis_multilocus(X)
      heb <- 2 * s$p * (1 - s$p); okb <- is.finite(heb) & heb > 0 & s$n >= 5
      data.frame(Population = pp, Fish = nrow(X), Ho = mean(s$Ho, na.rm = TRUE),
                 He = mean(s$He, na.rm = TRUE), FIS = unname(f["FIS"]),
                 `95% CI` = sprintf("%.3f to %.3f", f["lo"], f["hi"]),
                 `F (GenAlEx style)` = mean(1 - s$Ho[okb] / heb[okb]), check.names = FALSE)
    }))
    list(pv = pv, fis = fis, summ = summ, loci = loci)
  })
  a2_sig <- reactive({
    b <- a2_base()
    validate(need(length(input$a2_corr) == 1, "Choose a multiple-testing correction."))
    adj <- if (input$a2_corr == "none") b$pv else apply(b$pv, 2, p.adjust, method = input$a2_corr)
    adj <- matrix(adj, nrow = nrow(b$pv), dimnames = dimnames(b$pv))
    sig <- adj < 0.05; sig[is.na(sig)] <- FALSE
    list(adj = adj, sig = sig, nsig = rowSums(sig))
  })
  a2_keep <- reactive({
    b <- a2_base(); s <- a2_sig()
    validate(need(length(input$a2_rule) == 1, "Choose a rule for removing loci."))
    k <- as.integer(input$a2_rule)
    if (k == 0) b$loci else b$loci[s$nsig < k]
  })

  output$a2_summary <- renderTable(a2_base()$summ, digits = 3)
  output$a2_cmp_table <- renderTable({
    need_step("A1")
    keep <- a2_keep()
    compare_table(key_results(G, POP, rv$A1$inds, rv$A1$loci, REF_POPS),
                  key_results(G, POP, rv$A1$inds, keep, REF_POPS), "After A1 only", "After your HWE rule")
  })
  output$a2_counts <- renderTable({
    b <- a2_base(); s <- a2_sig()
    removed <- if (length(input$a2_rule) == 1) length(b$loci) - length(a2_keep()) else NA
    data.frame(Item = c(paste("Out of HWE in", REF_POPS), "Out of HWE in all three", "Loci removed by your rule"),
               Value = c(colSums(s$sig), sum(s$nsig == 3), removed))
  }, colnames = FALSE, digits = 0)
  output$a2_loci <- renderTable({
    b <- a2_base(); s <- a2_sig()
    w <- which(s$nsig > 0)
    validate(need(length(w) > 0, "No locus is out of HWE under this correction."))
    df <- data.frame(Locus = b$loci[w], `Populations out of HWE` = s$nsig[w], check.names = FALSE)
    for (pp in REF_POPS) {
      df[[paste("p", pp)]] <- formatC(s$adj[w, pp], digits = 2, format = "g")
      df[[paste("FIS", pp)]] <- fmt3(b$fis[w, pp])
    }
    df[order(-df[[2]]), ]
  })
  output$a2_plot <- renderPlot(res = 110, {
    b <- a2_base(); s <- a2_sig()
    x <- GPOS[b$loci]
    op <- par(mar = c(5, 4.5, 3, 1)); on.exit(par(op))
    plot(NA, xlim = range(GPOS), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "FIS",
         main = "Per-locus FIS along the genome; filled symbols are out of HWE")
    for (k in seq_along(CHR_LEVELS)) {
      ch <- CHR_LEVELS[k]
      x0 <- CHR_OFF[ch] / 1e6; x1 <- (CHR_OFF[ch] + CHR_MAX[ch]) / 1e6
      if (k %% 2 == 0) rect(x0, -1.2, x1, 1.2, col = "grey93", border = NA)
      axis(1, at = (x0 + x1) / 2, labels = CHR_LAB[ch], tick = FALSE, cex.axis = 0.7, las = 2)
    }
    mtext("Chromosome (RefSeq NC_0451xx.1); loci at their physical positions", side = 1, line = 3.5)
    abline(h = 0, col = "grey60")
    for (pp in REF_POPS)
      points(x, b$fis[, pp], col = COLS[pp], pch = ifelse(s$sig[, pp], 19, 1), cex = 0.9)
    legend("bottomright", REF_POPS, col = COLS[REF_POPS], pch = 19, bty = "n", bg = "white")
  })

  observeEvent(input$a2_commit, {
    if (length(input$a2_corr) == 0 || length(input$a2_rule) == 0) {
      showNotification("Choose a correction and a rule first.", type = "warning"); return()
    }
    if (!why_ok("a2")) return()
    keep <- a2_keep()
    rv$A2 <- list(loci = keep)
    rv$res_A2 <- a2_base()$summ
    add_log("A2 HWE & FIS",
      sprintf("Correction: %s; remove loci out of HWE in \u2265 %s populations (0 = never). Removed %d, kept %d loci.",
              input$a2_corr, input$a2_rule, length(rv$A1$loci) - length(keep), length(keep)),
      input$a2_why)
    advance("A2")
  })

  # ======================= A4: FST and PCA (references) ==================
  ref_groups <- function(loci) {
    inds <- rv$A1$inds
    setNames(lapply(REF_POPS, function(pp) G[ref_idx(inds, pp), loci, drop = FALSE]), REF_POPS)
  }
  fst_display <- function(fm) {
    k <- nrow(fm$FST)
    out <- matrix("", k, k, dimnames = dimnames(fm$FST))
    for (i in seq_len(k)) for (j in seq_len(k)) {
      if (i < j) out[i, j] <- sprintf("%.3f", fm$FST[i, j])
      if (i > j) out[i, j] <- if (fm$p[i, j] <= 1 / (fm$nperm + 1) + 1e-12)
        sprintf("p \u2264 %.3f", 1 / (fm$nperm + 1)) else sprintf("p = %.3f", fm$p[i, j])
    }
    data.frame(` ` = rownames(out), out, check.names = FALSE)
  }
  a4_fst <- reactive({
    need_step("A2")
    withProgress(message = "Permutation tests", value = 0.5, {
      set.seed(7)
      list(final = fst_matrix(ref_groups(rv$A2$loci)))
    })
  })
  output$a4_h_final <- renderUI({ need_step("A2"); h4(sprintf("Pairwise FST among the reference populations (%d loci)", length(rv$A2$loci))) })
  output$a4_fst <- renderTable(fst_display(a4_fst()$final))

  plot_pca <- function(X, grp, main) {
    p <- pca_geno(X)
    is_ref <- grp %in% REF_POPS
    o <- c(which(is_ref), which(!is_ref))          # Oro Bay fish drawn on top
    pcol <- ifelse(is_ref[o], adjustcolor(COLS[grp[o]], 0.45), adjustcolor(COLS[grp[o]], 0.9))
    plot(p$x[o, , drop = FALSE], col = pcol,
         pch = PCH[grp[o]], cex = ifelse(is_ref[o], 1, 1.35),
         xlab = sprintf("PC1 (%.1f%%)", 100 * p$var[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * p$var[2]),
         main = sprintf("%s (%d loci)", main, p$nloci))
    u <- names(COLS)[names(COLS) %in% grp]
    legend("topright", u, col = COLS[u], pch = PCH[u], pt.cex = ifelse(u %in% REF_POPS, 1, 1.35),
           bty = "n", cex = 0.85)
  }
  output$a4_pca <- renderPlot(res = 110, {
    need_step("A2")
    inds <- rv$A1$inds
    op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    plot_pca(G[inds, rv$A2$loci, drop = FALSE], POP[inds], "Reference populations")
  })

  observeEvent(input$a4_commit, {
    if (!why_ok("a4")) return()
    rv$A3 <- list(loci = rv$A2$loci)
    rv$res_A3 <- fst_display(a4_fst()$final)
    add_log("A3 FST & PCA", sprintf("Final locus set: %d loci.", length(rv$A2$loci)), input$a4_why)
    advance("A3")
  })

  # ======================= B1: split Oro Bay =============================
  oro <- which(IS_ORO)
  b1_def <- reactive({
    st <- samples$Stage[oro]; gsi <- samples$GSI[oro]; dt <- samples$CaptureDate[oro]
    ripe  <- (gsi >= input$b1_gsi) %in% TRUE & st %in% input$b1_late_st
    late  <- ripe & (dt >= input$b1_date) %in% TRUE
    early <- ripe & (dt < input$b1_date) %in% TRUE
    imm   <- st %in% input$b1_imm_st
    if (isTRUE(input$b1_imm_gsi)) imm <- imm & (gsi < input$b1_gsi) %in% TRUE
    both <- ripe & imm
    late[both] <- FALSE; early[both] <- FALSE; imm[both] <- FALSE
    grp <- ifelse(late, "Late spawners", ifelse(early, "Early spawners", ifelse(imm, "Immature", NA)))
    q <- CALLRATE[oro] >= input$b1_cr
    base_loci <- if (is.null(rv$A3)) LOCI else rv$A3$loci
    sel <- oro[!is.na(grp) & q]
    loc_cr <- if (length(sel)) colMeans(!is.na(G[sel, base_loci, drop = FALSE])) else setNames(rep(NA_real_, length(base_loci)), base_loci)
    loci <- base_loci[(loc_cr >= input$b1_loc) %in% TRUE]
    list(grp = grp, q = q, both = sum(both), loc_cr = loc_cr, loci = loci, base_loci = base_loci)
  })
  output$b1_counts <- renderTable({
    d <- b1_def()
    cnt <- function(g, q = TRUE) sum(d$grp == g & q, na.rm = TRUE)
    data.frame(Group = c(ORO_GROUPS, "Fish fitting two groups (left out)", "Loci kept (of your Part A set)"),
               `Meet definition` = c(sapply(ORO_GROUPS, cnt), d$both, length(d$base_loci)),
               `After call-rate filters` = c(sapply(ORO_GROUPS, cnt, q = d$q), NA, length(d$loci)),
               check.names = FALSE)
  }, digits = 0)
  output$b1_plots <- renderPlot(res = 110, {
    st <- samples$Stage[oro]; gsi <- samples$GSI[oro]; dt <- samples$CaptureDate[oro]
    op <- par(mfrow = c(3, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    hist(gsi, breaks = seq(0, 0.3, 0.01), col = "grey80", border = "white",
         main = "GSI of Oro Bay fish", xlab = "GSI", ylab = "Fish")
    abline(v = input$b1_gsi, col = "red", lwd = 2)
    plot(dt, gsi, pch = 19, col = adjustcolor(STAGE_COLS[st], 0.8), xlab = "Capture date", ylab = "GSI",
         main = "GSI by capture date (colour = spawning stage)")
    abline(h = input$b1_gsi, v = input$b1_date, col = "red", lty = 2)
    legend("topright", STAGES, col = STAGE_COLS, pch = 19, bty = "n", cex = 0.8, title = "Stage")
    stripchart(gsi ~ factor(st, levels = STAGES), vertical = TRUE, method = "jitter", jitter = 0.25,
               pch = 19, col = adjustcolor(STAGE_COLS, 0.7), xlab = "Spawning stage", ylab = "GSI",
               main = "GSI by spawning stage")
    abline(h = input$b1_gsi, col = "red", lty = 2)
    hist(CALLRATE[oro], breaks = seq(0, 1, 0.025), col = "grey80", border = "white",
         main = "Call rate per Oro Bay fish", xlab = "Share of loci genotyped", ylab = "Fish")
    abline(v = input$b1_cr, col = "red", lwd = 2)
    d <- b1_def()
    lc <- d$loc_cr
    if (all(is.na(lc))) {
      plot.new(); text(0.5, 0.5, "Define the groups to see locus call rates."); plot.new()
    } else {
      hist(lc, breaks = seq(0, 1, 0.025), col = "grey80", border = "white",
           main = "Call rate per locus in your Oro Bay fish", xlab = "Share of fish genotyped", ylab = "Loci")
      abline(v = input$b1_loc, col = "red", lwd = 2)
      ref_lc <- colMeans(!is.na(G[IS_REF, d$base_loci, drop = FALSE]))
      plot(ref_lc, lc, pch = 19, col = adjustcolor(ifelse(lc >= input$b1_loc, "grey30", "red"), 0.7),
           xlim = c(0, 1), ylim = c(0, 1), xlab = "Locus call rate in reference fish",
           ylab = "Locus call rate in your Oro Bay fish", main = "Same loci, two sequencing runs")
      abline(0, 1, col = "grey60", lty = 2); abline(h = input$b1_loc, col = "red", lty = 2)
    }
  })
  output$b1_table <- renderTable({
    st <- factor(samples$Stage[oro], levels = STAGES)
    mo <- factor(month.abb[samples$Month[oro]], levels = month.abb[sort(unique(samples$Month[oro]))])
    tab <- table(st, mo)
    data.frame(Stage = rownames(tab), as.data.frame.matrix(tab), check.names = FALSE)
  }, digits = 0)

  b1_issues <- reactive({
    d <- b1_def(); iss <- character()
    if ("VIII" %in% input$b1_imm_st)
      iss <- c(iss, "You count stage VIII fish as immature. Look at GSI by spawning stage and by date: are these fish immature, or have they already spawned?")
    if (sum(d$grp == "Early spawners", na.rm = TRUE) < 5)
      iss <- c(iss, "Your early-spawner group is (almost) empty. Look at GSI by capture date: where does the spawning season split?")
    if (input$b1_cr < 0.6)
      iss <- c(iss, "You keep fish with fewer than 60% of loci genotyped. Do these fish add information, or mostly missing data?")
    if (input$b1_loc < 0.5)
      iss <- c(iss, "Some loci are genotyped in only a small share of your Oro Bay fish. Compare locus call rates between the two sequencing runs: did your Part A filter catch these loci?")
    iss
  })
  b1_commit <- function() {
    d <- b1_def()
    all_g <- setNames(lapply(ORO_GROUPS, function(g) oro[d$grp %in% g]), ORO_GROUPS)
    q_g   <- setNames(lapply(ORO_GROUPS, function(g) oro[d$grp %in% g & d$q]), ORO_GROUPS)
    rv$B1 <- list(all = all_g, groups = q_g, loci = d$loci, other = oro[is.na(d$grp) & d$q])
    rv$res_B1 <- data.frame(Group = ORO_GROUPS, `Meet definition` = lengths(all_g),
                            `After call-rate filter` = lengths(q_g),
                            Loci = length(d$loci), check.names = FALSE)
    add_log("B1 Split Oro Bay",
      sprintf("Spawners: GSI \u2265 %.2f, stages %s; late = caught on or after %s, early = before. Immature: stages %s%s. Fish call rate \u2265 %.2f; locus call rate in Oro Bay \u2265 %.2f. Result: %d immature, %d early spawners, %d late spawners, %d loci.",
              input$b1_gsi, paste(input$b1_late_st, collapse = ","), format(input$b1_date),
              paste(input$b1_imm_st, collapse = ","), ifelse(isTRUE(input$b1_imm_gsi), " and GSI below cutoff", ""),
              input$b1_cr, input$b1_loc, length(q_g[[1]]), length(q_g[[2]]), length(q_g[[3]]), length(d$loci)),
      input$b1_why)
    advance("B1")
  }
  observeEvent(input$b1_commit, {
    need_A3 <- is.null(rv$A3)
    if (need_A3) { showNotification("Commit step A3 first.", type = "warning"); return() }
    if (!why_ok("b1")) return()
    d <- b1_def()
    n_g <- sapply(ORO_GROUPS, function(g) sum(d$grp == g & d$q, na.rm = TRUE))
    if (any(n_g < 5)) {
      showNotification("Each of the three groups needs at least 5 fish after the call-rate filter.", type = "error"); return()
    }
    if (length(d$loci) < 10) {
      showNotification("Fewer than 10 loci pass your locus filter. Relax it.", type = "error"); return()
    }
    iss <- b1_issues()
    if (length(iss)) {
      showModal(modalDialog(title = "Before you commit, consider:",
        tags$ul(lapply(iss, tags$li)),
        footer = tagList(modalButton("Go back and reconsider"),
                         actionButton("b1_force", "Commit anyway", class = "btn-warning"))))
    } else b1_commit()
  })
  observeEvent(input$b1_force, { removeModal(); b1_commit() })

  # ======================= B2: FIS in Oro Bay ============================
  b2_res <- reactive({
    need_step("B1")
    loci <- rv$B1$loci; inds <- rv$A1$inds
    sets <- c(setNames(lapply(REF_POPS, function(pp) ref_idx(inds, pp)), REF_POPS), rv$B1$groups)
    set.seed(11)
    tab <- do.call(rbind, lapply(names(sets), function(g) {
      X <- G[sets[[g]], loci, drop = FALSE]; s <- locus_stats(X); f <- fis_multilocus(X)
      pv <- hwe_table(X)
      data.frame(Group = g, Fish = nrow(X), Ho = mean(s$Ho, na.rm = TRUE), He = mean(s$He, na.rm = TRUE),
                 FIS = unname(f["FIS"]), `95% CI` = sprintf("%.3f to %.3f", f["lo"], f["hi"]),
                 `Loci out of HWE (p < 0.05)` = sum(pv < 0.05, na.rm = TRUE), check.names = FALSE)
    }))
    # (a) stricter call-rate filter on fish, with the B1 locus set
    thr <- seq(0.5, 0.95, 0.05)
    sens_fish <- lapply(ORO_GROUPS, function(g) {
      all <- rv$B1$all[[g]]
      t(vapply(thr, function(t) {
        ii <- all[CALLRATE[all] >= t]
        if (length(ii) < 4) return(c(NA_real_, as.numeric(length(ii))))
        c(unname(fis_multilocus(G[ii, loci, drop = FALSE], B = 2)["FIS"]), as.numeric(length(ii)))
      }, numeric(2)))
    })
    # (b) stricter call-rate filter on loci (within each group), starting from the Part A loci
    thr_l <- seq(0, 0.95, 0.05)
    sens_loc <- lapply(ORO_GROUPS, function(g) {
      ii <- rv$B1$groups[[g]]
      X <- G[ii, rv$A3$loci, drop = FALSE]; cr <- colMeans(!is.na(X))
      t(vapply(thr_l, function(t) {
        k <- cr >= t
        if (sum(k) < 5) return(c(NA_real_, as.numeric(sum(k))))
        c(unname(fis_multilocus(X[, k, drop = FALSE], B = 2)["FIS"]), as.numeric(sum(k)))
      }, numeric(2)))
    })
    names(sens_fish) <- names(sens_loc) <- ORO_GROUPS
    list(tab = tab, sens_fish = sens_fish, thr = thr, sens_loc = sens_loc, thr_l = thr_l)
  })
  output$b2_table <- renderTable(b2_res()$tab, digits = 3)
  sens_panel <- function(thr, sens, xlab, main, what) {
    rng <- range(c(0, unlist(lapply(sens, function(m) m[, 1]))), na.rm = TRUE)
    plot(NA, xlim = range(thr), ylim = rng + c(-0.05, 0.08), xlab = xlab, ylab = "Multilocus FIS", main = main)
    abline(h = 0, col = "grey60")
    for (g in names(sens)) {
      m <- sens[[g]]
      lines(thr, m[, 1], col = COLS[g], lwd = 2, type = "b", pch = PCH[g], cex = 1.2)
      text(thr, m[, 1], labels = ifelse(is.na(m[, 1]), "", m[, 2]), pos = 3, cex = 0.65, col = COLS[g])
    }
    legend("topright", names(sens), col = COLS[names(sens)], lwd = 2, pch = PCH[names(sens)], bty = "n",
           title = paste("Numbers =", what), cex = 0.85)
  }
  output$b2_plot <- renderPlot(res = 110, {
    r <- b2_res()
    op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    sens_panel(r$thr, r$sens_fish, "Minimum call rate of fish included", "Filtering fish", "fish included")
    sens_panel(r$thr_l, r$sens_loc, "Minimum call rate of loci (within each group)", "Filtering loci", "loci included")
  })
  observeEvent(input$b2_commit, {
    if (length(input$b2_imm) == 0 || length(input$b2_early) == 0 || length(input$b2_late) == 0) {
      showNotification("Answer all three questions first.", type = "warning"); return()
    }
    if (!why_ok("b2")) return()
    rv$B2 <- TRUE
    rv$res_B2 <- b2_res()$tab
    add_log("B2 FIS in Oro Bay", sprintf("Immature: %s. Early spawners: %s. Late spawners: %s.",
                                         input$b2_imm, input$b2_early, input$b2_late),
            input$b2_why)
    advance("B2")
  })

  # ======================= B3: FST and PCA with Oro Bay ==================
  b3_sets <- reactive({
    need_step("B2")
    inds <- rv$A1$inds
    c(setNames(lapply(REF_POPS, function(pp) ref_idx(inds, pp)), REF_POPS), rv$B1$groups)
  })
  b3_fst <- reactive({
    sets <- b3_sets(); loci <- rv$B1$loci
    withProgress(message = "Permutation tests", value = 0.5, {
      set.seed(13)
      fst_matrix(lapply(sets, function(ii) G[ii, loci, drop = FALSE]))
    })
  })
  output$b3_h <- renderUI({ need_step("B2"); h4(sprintf("Pairwise FST (upper triangle) and permutation p-values (lower triangle), %d loci", length(rv$B1$loci))) })
  output$b3_fst <- renderTable(fst_display(b3_fst()))
  b3_assign <- reactive({
    sets <- b3_sets(); loci <- rv$B1$loci
    validate(need(length(input$b3_base) >= 2, "Choose at least two baseline groups."))
    base <- lapply(sets[input$b3_base], function(ii) G[ii, loci, drop = FALSE])
    qnames <- setdiff(names(sets), input$b3_base)
    qsets <- sets[qnames]
    if (isTRUE(input$b3_other) && length(rv$B1$other)) qsets[["Not grouped"]] <- rv$B1$other
    qi <- unlist(qsets); qg <- rep(names(qsets), lengths(qsets))
    res <- assign_fish(base, if (length(qi)) G[qi, loci, drop = FALSE] else NULL)
    rows <- c(lapply(names(res$base), function(k) data.frame(group = k, res$base[[k]], check.names = FALSE)),
              if (!is.null(res$query)) list(data.frame(group = qg, res$query, check.names = FALSE)))
    df <- do.call(rbind, rows)
    names(df) <- c("group", res$K)
    df$group <- factor(df$group, levels = unique(c(res$K, names(qsets))))
    df
  })
  output$b3_assign_plot <- renderPlot(res = 110, {
    df <- b3_assign(); K <- setdiff(names(df), "group")
    cols <- c(COLS, `Not grouped` = "grey60")[K]
    groups <- levels(df$group)[levels(df$group) %in% df$group]
    sizes <- as.numeric(table(df$group)[groups])
    gap <- max(4, round(sum(sizes) / 40))
    op <- par(mar = c(3.8, 4.5, 1, 1)); on.exit(par(op))
    M <- NULL; mids <- numeric()
    for (g in groups) {
      sub <- as.matrix(df[df$group == g, K, drop = FALSE])
      top <- max.col(sub, ties.method = "first")
      sub <- sub[order(top, -apply(sub, 1, max)), , drop = FALSE]
      start <- if (is.null(M)) 0 else ncol(M)
      M <- cbind(M, t(sub), matrix(0, length(K), gap))
      mids <- c(mids, start + nrow(sub) / 2)
    }
    barplot(M, col = cols, border = NA, space = 0, ylab = "Assignment probability", axes = TRUE,
            ylim = c(0, 1.32), yaxt = "n")
    axis(2, at = seq(0, 1, 0.25))
    # group labels above the bars, alternating heights so narrow groups do not overlap
    yl <- ifelse(seq_along(mids) %% 2 == 1, 1.06, 1.18)
    text(mids, yl, sprintf("%s (%d)", groups, sizes), cex = 0.95, font = 2)
    legend("bottom", inset = c(0, -0.1), xpd = TRUE, K, fill = cols, horiz = TRUE, bty = "n", cex = 0.9)
  })
  output$b3_assign_table <- renderTable({
    df <- b3_assign(); K <- setdiff(names(df), "group")
    top <- factor(K[max.col(as.matrix(df[, K]), ties.method = "first")], levels = K)
    tab <- as.data.frame.matrix(table(df$group, top))
    conf <- tapply(apply(as.matrix(df[, K]), 1, max) > 0.9, df$group, mean)
    data.frame(Group = rownames(tab), tab, `Share with p > 0.9` = sprintf("%.0f%%", 100 * conf[rownames(tab)]),
               check.names = FALSE)
  }, digits = 0)
  output$b3_pca <- renderPlot(res = 110, {
    sets <- b3_sets()
    ii <- unlist(sets); grp <- rep(names(sets), lengths(sets))
    op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    plot_pca(G[ii, rv$B1$loci, drop = FALSE], grp, "Reference populations and Oro Bay groups")
  })
  observeEvent(input$b3_commit, {
    if (length(input$b3_late) == 0 || length(input$b3_early) == 0 || length(input$b3_imm) == 0 ||
        length(input$b3_assign) == 0) {
      showNotification("Answer all three questions first.", type = "warning"); return()
    }
    if (!why_ok("b3")) return()
    rv$B3 <- TRUE
    rv$res_B3 <- fst_display(b3_fst())
    rv$res_B3b <- isolate(tryCatch({
      df <- b3_assign(); K <- setdiff(names(df), "group")
      top <- factor(K[max.col(as.matrix(df[, K]), ties.method = "first")], levels = K)
      tab <- as.data.frame.matrix(table(df$group, top)); data.frame(Group = rownames(tab), tab, check.names = FALSE)
    }, error = function(e) NULL))
    add_log("B3 FST, PCA & assignment",
            sprintf("Late spawners distinct: %s. Early spawners: %s (assignment: %s). Immature fish: %s. Assignment baselines: %s.",
                    input$b3_late, input$b3_early, input$b3_assign, input$b3_imm, paste(input$b3_base, collapse = ", ")),
            input$b3_why)
    advance("B3")
  })

  # ======================= B4: synthesis and report ======================
  output$log_table <- renderTable(rv$log)
  output$report <- downloadHandler(
    filename = function() {
      nm <- gsub("[^A-Za-z0-9]+", "_", trimws(input$student))
      paste0("herring_report_", ifelse(nchar(nm) > 0, nm, "student"), ".html")
    },
    content = function(file) {
      esc <- htmltools::htmlEscape
      sec <- function(title, df) paste0("<h2>", esc(title), "</h2>", html_table(df))
      html <- paste0(
        "<!doctype html><html><head><meta charset='utf-8'><title>Herring GT-seq report</title>",
        "<style>body{font-family:sans-serif;max-width:960px;margin:24px auto;padding:0 16px;color:#222}",
        "table{border-collapse:collapse;margin:8px 0 20px}td,th{border:1px solid #ccc;padding:4px 8px;",
        "font-size:13px;vertical-align:top;text-align:left}th{background:#f2f2f2}</style></head><body>",
        "<h1>Herring GT-seq analysis report</h1>",
        "<p><strong>Student:</strong> ", esc(input$student), "<br><strong>Date:</strong> ",
        format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>",
        sec("Decision log", rv$log),
        sec("A2: reference populations (HWE & FIS)", rv$res_A2),
        sec("A3: pairwise FST among references (upper) and p-values (lower)", rv$res_A3),
        sec("B1: Oro Bay groups", rv$res_B1),
        sec("B2: FIS by group", rv$res_B2),
        sec("B3: pairwise FST with Oro Bay groups (upper) and p-values (lower)", rv$res_B3),
        sec("B3: assignment (number of fish assigned to each baseline)", rv$res_B3b),
        "<h2>Conclusion</h2><p>", gsub("\n", "<br>", esc(input$b4_text)), "</p>",
        "</body></html>")
      writeLines(html, file, useBytes = TRUE)
    })

  # ======================= Instructor notes ==============================
  notes <- list(
    A1 = "Reference fish are well genotyped (median ~96% of loci). With fish and locus call rate \u2265 0.8 and MAF \u2265 0.01 (pooled), 92 fish and 120 of 133 loci remain; about nine loci are monomorphic in the reference fish. Note: two loci that pass here (NC_045152.1:21070884, NC_045163.1:4795125) are rarely genotyped in the Oro Bay run; students should catch them in B1.",
    A2 = "Reference populations are close to HWE: multilocus FIS about 0.02\u20130.06 (SquaxinPass and ElliottBay highest; their 95% CIs should be checked in the table). The app uses He corrected for sample size; GenAlEx's F (mean over loci of 1 \u2212 Ho/He, uncorrected He) is about 0.015 lower, i.e. \u22120.002 to 0.04 (column 'F (GenAlEx style)'). SquaxinPass (0.06, 95% CI 0.04\u20130.09) and ElliottBay (0.06, 0.03\u20130.10) are slightly but significantly positive; Cherry14 (0.02, \u22120.02\u20130.05) is not. The per-locus median is near 0 (0.01 and \u22120.02); the mean is pulled up by a tail of loci with F 0.3\u20130.8, so this is a good place to discuss null alleles or slight substructure within a sample. Without correction 2\u20134 loci per population fail HWE; with Bonferroni or FDR only one locus fails, in one population. No locus fails in all three, so the removal rule matters little here.",
    A3 = "With the intended A1/A2 choices (120 loci), FST among references is 0.15 (SquaxinPass\u2013Cherry14), 0.20 (SquaxinPass\u2013ElliottBay) and 0.32 (Cherry14\u2013ElliottBay); all p \u2264 0.001. Values are high because the panel was designed for stock identification. The PCA separates the three populations clearly; ElliottBay is the most distinct.",
    B1 = "Intended groups: GSI \u2265 0.05; late spawners caught from 1 April, early spawners before; fish call rate \u2265 0.7. Late spawners: 43 meet the definition, 30 after the call-rate filter (31 if one fish at 69.9% is kept). Early spawners: 47 meet the definition but only 19 pass (early-season fish were genotyped poorly). Immature = stage II only (15; 12 pass). Locus call rate \u2265 0.5 in these fish leaves 111 of 120 loci and removes two loci that are rarely genotyped in the Oro Bay run but pass in the references. Decisions here matter most: with GSI cutoff 0 (spent fish included), the late group (71 fish) has FIS \u2248 0.11, FST to SquaxinPass roughly halves, and in the assignment test about a third of its fish assign elsewhere. Trap: stage VIII fish are spent and have low GSI like immature fish.",
    B2 = "Expected with the intended choices: late spawners FIS \u2248 0 (\u22120.01), consistent with one population (and with GenAlEx on the earlier lab file). Immature fish FIS \u2248 0.24 and early spawners \u2248 0.22, both stable whether fish or loci are filtered more strictly: a Wahlund effect from mixed stocks, not a missing-data artefact. Heterozygosity is unrelated to call rate (r = \u22120.09). Reference FIS 0.01\u20130.08.",
    B3 = "Expected (111 loci): late spawners differ from all references (FST \u2248 0.23 vs SquaxinPass, 0.22 vs Cherry14, 0.54 vs ElliottBay) and form their own PCA cluster: a new population, and not Cherry Point fish. Early spawners sit between SquaxinPass (\u2248 0.04) and late spawners (\u2248 0.10); immature fish are closest to SquaxinPass (\u2248 0.01) and early spawners (\u2248 0.03). Assignment with the three references only forces late spawners into a reference group. With late spawners added as a baseline: baseline fish self-assign well (SquaxinPass 28/32, Cherry14 27/29, ElliottBay 29/31, late 30/30). Early spawners: 12 of 19 assign to late spawners, 4 SquaxinPass, 3 ElliottBay; immature: 6 SquaxinPass, 4 ElliottBay, 2 late. No Oro Bay fish assigns to Cherry Point. The result hardly depends on the date cutoff, because the 7 non-late fish were all caught before mid-February. Spent (stage VIII) fish left out in B1 are about half SquaxinPass-like.",
    B4 = "A good conclusion: Oro Bay late spawners are a distinct population in Hardy\u2013Weinberg equilibrium; immature fish are a mixed aggregation showing a Wahlund effect. Good answers also mention the quality control: loci that work in one sequencing run can fail in another.")
  for (s in names(notes)) local({
    ss <- s
    output[[paste0("note_", ss)]] <- renderUI({
      if (!instructor()) return(NULL)
      div(class = "instructor", strong("Instructor note: "), notes[[ss]])
    })
  })
}

shinyApp(ui, server)
