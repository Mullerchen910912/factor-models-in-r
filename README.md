# Factor models in R — hands-on

**Practical R walkthroughs for factor models, from statistical factor extraction
(PCA) to a modern causal factor model — all on public data.**

Factor models capture the idea that many observed series are driven by a few
unobserved common forces. This repo works through the main flavours by hand.

## Contents

### 1. Statistical factor models via PCA
- **`hands-on_R_factors/sp500/pca_sp500.R`** — extract latent factors from S&P 500
  stock returns with principal components; interpret how many factors matter and
  what they represent ([`sp500_prices.csv`](hands-on_R_factors/sp500/sp500_prices.csv) included).
- **`hands-on_R_factors/str_pca/str_pca.R`** — PCA on the California schools
  dataset (a classic teaching dataset) to illustrate dimension reduction.
- **`hands-on_R_factors/R_script_for_factor_analysis.R`** — general factor-analysis
  routines.

### 2. Causal factor model — Bai & Wang (2024)
- **`bai_wang_2024_replication.R`** — replicates the Bai & Wang (2024) approach to
  treatment-effect estimation under a factor structure.
- **`bai_wang_2024_walkthrough.Rmd` / `.pdf`** — an annotated walkthrough of the
  method with worked output (`fig_ca_*` show the estimated counterfactual and
  treatment effect).

## Run it

Open any script in RStudio and run top-to-bottom, or:

```r
source("hands-on_R_factors/sp500/pca_sp500.R")   # PCA factors from S&P 500 returns
source("bai_wang_2024_replication.R")            # Bai & Wang (2024) causal factor model
rmarkdown::render("bai_wang_2024_walkthrough.Rmd")
```

All data here is public (S&P 500 prices; California schools). Typical packages:
`tidyverse`, `stats` (`prcomp`), and base R.

## References

- Bai, J. & Wang, P. (2024). *Causal inference with factor models.*
- Tsay, R. (2010). *Analysis of Financial Time Series*, Ch. 9 (statistical factor models).

---

*By Sheng-Yuan Chen — MS Economics, National Taiwan University.*
