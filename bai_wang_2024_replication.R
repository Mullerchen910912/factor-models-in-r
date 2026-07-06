# =============================================================================
# Replication: Bai & Wang (2024) "Causal Inference Using Factor Models"
# Application 1: Abadie et al. (2010) — California Prop 99
#
# Algorithm (Section 5.1 — no covariates):
#   Step 1: PCA on demeaned control units -> F_hat (T x r)
#   Step 2a: regress Y_treated[pre]  ~ 1 + F_hat[pre]  -> (mu_0, lambda_0)
#   Step 2b: regress Y_treated[post] ~ 1 + F_hat[post] -> (mu_1, lambda_1)
#   Step 3:  counterfactual Y_CF = mu_0 + lambda_0' F_hat_t (all t)
#            causal effect tau_t = Y_treated,t - Y_CF_t  (post only)
#
# NOTE: Column-demeaning before PCA extracts "common trend" factors.
#       The intercept in the loading regression absorbs each unit's level.
#       This gives ATT = -22.2, closely matching SCM (-22.0).
# =============================================================================

library(tidyverse)
library(tidysynth)
library(strucchange)

set.seed(2024)

fig_dir <- "/Users/chenshengyuan/Desktop/NTU_ECON/Course/Advanced econometric topic/factor_models"

# =============================================================================
# Helper: Bai-Ng (2002) IC1 criterion to select number of factors
# =============================================================================
bai_ng_ic1 <- function(Y, r_max = 8) {
  T_n  <- nrow(Y)
  N_n  <- ncol(Y)
  C2   <- min(T_n, N_n)
  Y_dm <- scale(Y, center = TRUE, scale = FALSE)
  sv   <- svd(Y_dm)

  ic1 <- numeric(r_max)
  for (r in 1:r_max) {
    F_r      <- sv$u[, 1:r, drop = FALSE] * sqrt(T_n)
    Lambda_r <- t(Y_dm) %*% F_r / T_n
    E_r      <- Y_dm - F_r %*% t(Lambda_r)
    V_r      <- sum(E_r^2) / (T_n * N_n)
    ic1[r]   <- log(V_r) + r * (T_n + N_n) / (T_n * N_n) * log(C2)
  }
  which.min(ic1)
}

# =============================================================================
# Helper: variance estimator for causal effect (Proposition 1 in paper)
#   var(tau_hat_it) = z_it' var(delta_hat_i) z_it + a_hat_i' var(f_hat_t) a_hat_i
#   where z_it = [f_hat_t; X_it],  a_hat_i = lambda_1 - lambda_0
#   (no covariates: z_it = f_hat_t)
# =============================================================================
tau_variance <- function(Y_treated_post, F_post, lambda_0, lambda_1,
                         Y_treated_pre,  F_pre,
                         F_hat, V_hat_inv, Gamma_hat) {
  # delta_hat_i = lambda_1 - lambda_0
  delta <- lambda_1 - lambda_0
  T_post <- length(Y_treated_post)
  T_pre  <- length(Y_treated_pre)

  # var(delta_hat_i) = var(lambda_hat_1) + var(lambda_hat_0)
  # Using HC0 sandwich for each period
  var_lambda <- function(Y_vec, F_mat) {
    X  <- cbind(1, F_mat)
    Xe <- X * residuals(lm(Y_vec ~ F_mat))
    bread  <- solve(t(X) %*% X / length(Y_vec))
    meat   <- t(Xe) %*% Xe / length(Y_vec)^2
    bread %*% meat %*% bread
  }
  V0 <- var_lambda(Y_treated_pre,  F_pre)[-1, -1]   # drop intercept row/col
  V1 <- var_lambda(Y_treated_post, F_post)[-1, -1]
  V_delta <- V0 + V1

  # var(tau_hat_it) for each post-period t
  var_tau <- numeric(T_post)
  for (t in 1:T_post) {
    f_t <- F_hat[nrow(F_hat) - T_post + t, ]
    var_tau[t] <- t(f_t) %*% V_delta %*% f_t
  }
  sqrt(var_tau)
}

# =============================================================================
# APPLICATION 1: Abadie et al. (2010) — California Prop 99
# =============================================================================
cat(strrep("=", 60), "\n")
cat("APPLICATION 1: California Prop 99 (Abadie et al. 2010)\n")
cat(strrep("=", 60), "\n\n")

data(smoking)

# ---- Reshape to T x N wide matrix ----
Y_wide <- smoking |>
  select(state, year, cigsale) |>
  pivot_wider(names_from = state, values_from = cigsale) |>
  arrange(year)

years_ca  <- Y_wide$year
CA_col    <- which(names(Y_wide) == "California")
col_order <- c(CA_col, setdiff(2:ncol(Y_wide), CA_col))
Y_ca      <- as.matrix(Y_wide[, col_order])
Y_control <- Y_ca[, -1]        # T x 38
Y_treated <- Y_ca[,  1]        # T x 1 (California)

T_n      <- nrow(Y_ca)
T0_idx   <- which(years_ca == 1989)
pre_idx  <- 1:(T0_idx - 1)     # 1970–1988 (19 years)
post_idx <- T0_idx:T_n          # 1989–2000 (12 years)

cat(sprintf("Units: %d treated (CA) + %d controls\n", 1, ncol(Y_control)))
cat(sprintf("T: %d (%d–%d) | T0: %d | pre: %d yrs | post: %d yrs\n\n",
            T_n, min(years_ca), max(years_ca), 1989,
            length(pre_idx), length(post_idx)))

# ---- Bai-Ng IC1 ----
r_ic1 <- bai_ng_ic1(Y_control)
r     <- 3   # paper uses r = 3
cat(sprintf("Bai-Ng IC1 selects r = %d  |  Using r = %d (as in paper)\n\n",
            r_ic1, r))

# ---- Step 1: PCA on demeaned control units ----
Y_c_dm <- scale(Y_control, center = TRUE, scale = FALSE)
sv_c   <- svd(Y_c_dm)
F_hat  <- sv_c$u[, 1:r] * sqrt(T_n)   # T x r, F_hat'F_hat = T * I_r

# ---- Step 2: Regress CA on factors (with intercept) ----
fit_pre  <- lm(Y_treated[pre_idx]  ~ F_hat[pre_idx, ])
fit_post <- lm(Y_treated[post_idx] ~ F_hat[post_idx, ])

lambda_0 <- coef(fit_pre)[-1]   # r-vector (slope only)
mu_0     <- coef(fit_pre)[1]    # intercept (level)
lambda_1 <- coef(fit_post)[-1]
mu_1     <- coef(fit_post)[1]

cat(sprintf("mu_0 = %.2f  |  lambda_0 = [%.3f, %.3f]\n", mu_0, lambda_0[1], lambda_0[2]))
cat(sprintf("mu_1 = %.2f  |  lambda_1 = [%.3f, %.3f]\n", mu_1, lambda_1[1], lambda_1[2]))
cat(sprintf("delta = lambda_1 - lambda_0 = [%.3f, %.3f]\n\n",
            (lambda_1 - lambda_0)[1], (lambda_1 - lambda_0)[2]))

# ---- Step 3: Counterfactual and causal effect ----
Y_CF    <- as.vector(mu_0 + F_hat %*% lambda_0)   # T x 1
tau_all <- Y_treated - Y_CF                        # gap for all t

cat("Post-treatment causal effects:\n")
for (i in seq_along(post_idx)) {
  cat(sprintf("  %d: %+.2f\n", years_ca[post_idx[i]], tau_all[post_idx[i]]))
}
cat(sprintf("\nFactor model ATT (1989–2000): %+.2f\n\n", mean(tau_all[post_idx])))

# =============================================================================
# STRUCTURAL BREAK TEST
# H0: tau_it = 0  <=>  H0: lambda(1) = lambda(0)
# Regress Y_CA on F_hat (full sample), test for break at T0 = 1989
# =============================================================================
cat(strrep("=", 60), "\n")
cat("STRUCTURAL BREAK TEST\n")
cat(strrep("=", 60), "\n\n")

df_sb <- data.frame(y = Y_treated, f1 = F_hat[, 1], f2 = F_hat[, 2])

# Chow test at known break T0 = 1989
chow <- sctest(y ~ f1 + f2, data = df_sb, type = "Chow", point = T0_idx)
cat(sprintf("Chow test at T0 = 1989:\n"))
cat(sprintf("  F-statistic = %.2f\n", chow$statistic))
cat(sprintf("  p-value     = %.4f\n\n", chow$p.value))

# QLR (sup-F) test: unknown break, 15% trimming
fs  <- Fstats(y ~ f1 + f2, data = df_sb, from = 0.15)
qlr <- sctest(fs, type = "supF")
cat(sprintf("QLR (sup-F) test with 15%% trimming:\n"))
cat(sprintf("  sup-F statistic = %.2f\n", qlr$statistic))
cat(sprintf("  p-value         = %.4f\n",  qlr$p.value))
cat(sprintf("  Max F at year   = %d\n\n",  years_ca[fs$breakpoint]))

# ---- Synthetic control for comparison ----
sc_ca <- smoking |>
  synthetic_control(
    outcome           = cigsale,
    unit              = state,
    time              = year,
    i_unit            = "California",
    i_time            = 1988,
    generate_placebos = FALSE
  ) |>
  generate_predictor(time_window = 1970:1988,
                     cigsale_avg = mean(cigsale, na.rm = TRUE)) |>
  generate_predictor(time_window = 1975, cigsale_1975 = mean(cigsale)) |>
  generate_predictor(time_window = 1980, cigsale_1980 = mean(cigsale)) |>
  generate_predictor(time_window = 1988, cigsale_1988 = mean(cigsale)) |>
  generate_weights(optimization_window = 1970:1988) |>
  generate_control()

synth_tbl <- sc_ca |>
  grab_synthetic_control() |>
  rename(year = time_unit, Y_actual = real_y, Y_synth = synth_y) |>
  mutate(tau_synth = Y_actual - Y_synth)

synth_att <- mean(synth_tbl$tau_synth[synth_tbl$year >= 1989])
cat(sprintf("Synthetic control ATT (1989–2000): %+.2f\n\n", synth_att))

# ---- Approximate SE from sandwich variance ----
se_tau <- tau_variance(
  Y_treated_post = Y_treated[post_idx],
  F_post         = F_hat[post_idx, ],
  lambda_0       = lambda_0,
  lambda_1       = lambda_1,
  Y_treated_pre  = Y_treated[pre_idx],
  F_pre          = F_hat[pre_idx, ],
  F_hat          = F_hat,
  V_hat_inv      = NULL,
  Gamma_hat      = NULL
)

# ---- Build results tibble ----
results_ca <- tibble(
  year       = years_ca,
  Y_actual   = Y_treated,
  Y_factor   = Y_CF,
  tau_factor = tau_all,
  se_factor  = c(rep(NA, length(pre_idx)), se_tau),
  post       = year >= 1989
) |>
  left_join(synth_tbl |> select(year, Y_synth, tau_synth), by = "year")

# =============================================================================
# FIGURE A: Counterfactual comparison
# =============================================================================
p_cf <- ggplot(results_ca, aes(x = year)) +
  geom_line(aes(y = Y_actual, color = "Actual CA"),   linewidth = 0.9) +
  geom_line(aes(y = Y_factor, color = "Factor model (counterfactual)"),
            linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = Y_synth, color = "Synthetic CA"), shape = 3, size = 2.5) +
  geom_vline(xintercept = 1988.5, color = "steelblue", linewidth = 0.7) +
  annotate("text", x = 1985.5, y = 125,
           label = "Passage of\nProposition 99 →",
           hjust = 1, size = 3.2, color = "grey40") +
  scale_color_manual(
    values = c("Actual CA"                     = "#2171b5",
               "Factor model (counterfactual)" = "#d95f0e",
               "Synthetic CA"                  = "#d95f0e"),
    name = NULL
  ) +
  labs(
    title    = "Causal Factor Model vs. Synthetic Control: Counterfactual California",
    subtitle = sprintf("r = %d factors  |  Factor ATT = %.1f  |  Synth ATT = %.1f packs",
                       r, mean(tau_all[post_idx]), synth_att),
    x = "Year", y = "Per-capita cigarette sales (packs)",
    caption = "Source: Abadie et al. (2010). Method: Bai & Wang (2024)."
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"),
        legend.position  = "bottom")

ggsave(file.path(fig_dir, "fig_ca_counterfactual.png"),
       p_cf, width = 9, height = 5.5, dpi = 150, bg = "white")

# =============================================================================
# FIGURE B: Causal effects with 95% CI
# =============================================================================
p_gap <- ggplot(results_ca, aes(x = year)) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 1988.5, color = "steelblue", linewidth = 0.7) +
  # 95% CI (factor model)
  geom_ribbon(data = . %>% filter(post),
              aes(ymin = tau_factor - 1.96 * se_factor,
                  ymax = tau_factor + 1.96 * se_factor),
              fill = "#d95f0e", alpha = 0.2) +
  # Factor model causal effect
  geom_line(aes(y = tau_factor, color = "Factor model"), linewidth = 0.9) +
  # Synthetic control gap
  geom_point(aes(y = tau_synth, color = "Synthetic control"), shape = 3, size = 2.5) +
  annotate("text", x = 1985.5, y = -2,
           label = "Passage of\nProposition 99 →",
           hjust = 1, size = 3.2, color = "grey40") +
  scale_color_manual(
    values = c("Factor model"       = "#d95f0e",
               "Synthetic control"  = "#2171b5"),
    name = NULL
  ) +
  labs(
    title    = "Causal Factor Model vs. Synthetic Control: Causal Effects",
    subtitle = "Shaded region: 95% CI (sandwich SE, Proposition 1)",
    x = "Year", y = "Per-capita cigarette sales (packs)",
    caption = "Negative = Prop 99 reduced cigarette consumption. Source: Abadie et al. (2010)."
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"),
        legend.position  = "bottom")

ggsave(file.path(fig_dir, "fig_ca_causal_effect.png"),
       p_gap, width = 9, height = 5, dpi = 150, bg = "white")

cat("Figures saved:\n")
cat("  fig_ca_counterfactual.png\n")
cat("  fig_ca_causal_effect.png\n")
cat("\n✅ Application 1 complete.\n")
