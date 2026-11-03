# Herring GT-seq: guided analysis (Shiny app)

A teaching app that walks students through a population-genetic analysis of
Pacific herring GT-seq data in two parts:

- **Part A (reference populations: SquaxinPass (Case Inlet), Cherry14, ElliottBay):** data QC,
  Hardy–Weinberg tests and FIS, FST and PCA.
  Students build the locus set they will use in Part B.
- **Part B (Oro Bay):** students split the Oro Bay fish into immature fish,
  early spawners and late spawners using GSI, spawning stage, capture date and
  call rate (per fish and per locus), then
  estimate FIS, FST, a PCA and an individual assignment test (a STRUCTURE-like
  plot) together with the references.

Each step ends with a decision and a written justification. Committing a step
unlocks the next one; recommitting an earlier step clears the later ones.
At the end, students download an HTML report with every decision,
justification and key result.

## Running it

Requires R (4.1 or later) and the `shiny` package only.

```r
install.packages("shiny")
shiny::runApp("path/to/HerringApp")   # or open app.R in RStudio and click 'Run App'
```

Do not run app.R with `source()` or line by line: the app must be started from its folder.

**Instructor mode:** add `?instructor=1` to the app URL
(for example `http://127.0.0.1:1234/?instructor=1`). All tabs open at once and
each step shows an instructor note with the expected results.

**Hosting for a class:** with an account on shinyapps.io,
`rsconnect::deployApp("HerringApp")`. A university Posit Connect server works
the same way.

## Files

- `app.R`: everything in one file: population-genetic functions in base R (exact HWE test,
  multilocus FIS with bootstrap, Hudson FST with permutation tests, PCA,
  individual assignment with leave-one-out for baseline fish).
- `data/genotypes.csv`: 288 fish × 133 loci (96 reference fish, 192 Oro Bay fish);
  genotypes as letter pairs, empty = missing.
- `data/samples.csv`: one row per fish: source, population, lab-reported call
  rate, capture date, sex, age, spawning stage and GSI (Oro Bay fish only).
- `data/loci.csv`: locus name, chromosome (RefSeq) and position.

## Data notes

- Loci are the 133 shared by the reference (VAL) and Oro Bay (Prod1) runs.
  Two of them (NC_045152.1:21070884, NC_045163.1:4795125) are well genotyped
  in the reference fish but called in fewer than 20% of Oro Bay fish. Filtering
  on the reference fish does not catch them; the locus call-rate filter for
  Oro Bay fish in step B1 does.
- Oro Bay fish are linked to the metadata through their sample number
  (`orobay_259` = `22OB- 259`). Ripe fish caught before the date cutoff form a third group, early spawners,
  which students should recognise as a mixture.

## Expected results with the intended choices

| Step | Result |
| --- | --- |
| A1 | Call rate ≥ 0.8 (fish and loci), MAF ≥ 0.01: 92 fish, 120 loci. The comparison panel shows that filtering barely changes FIS or FST. |
| A2 | Few loci out of HWE after correction; reference FIS 0.02–0.06 (GenAlEx-style F about 0.015 lower) |
| A3 | FST 0.15 (SquaxinPass–Cherry14), 0.20 (SquaxinPass–ElliottBay), 0.32 (Cherry14–ElliottBay) |
| B1 | GSI ≥ 0.05 and fish call rate ≥ 0.7: late spawners caught from 1 April (30 fish), early spawners before (19 fish); immature: stage II (12 fish); locus call rate ≥ 0.5 in these fish: 111 loci |
| B2 | FIS: late spawners ≈ 0; immature ≈ 0.24 and early spawners ≈ 0.22, stable across fish and locus filters (Wahlund effect) |
| B3 | Late spawners FST ≈ 0.23 vs SquaxinPass, 0.22 vs Cherry14, 0.54 vs ElliottBay; early spawners between SquaxinPass (0.04) and late spawners (0.10); immature closest to SquaxinPass (0.01). Assignment with late spawners as a baseline: early spawners 12 of 19 late-type, rest SquaxinPass/ElliottBay-like; immature mostly SquaxinPass/ElliottBay-like; no Oro Bay fish assigns to Cherry Point |

## Status

The calculations were checked against an independent implementation, but the
app was first run in R by the author of the exercise; please run it again after
every change.
