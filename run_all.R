# Reproduce the results of run_all.py with the base-R implementation (R/dnsmodel.R).
#
#   Rscript run_all.R           # full run, including the Kalman backtest
#   Rscript run_all.R --fast    # skip the Kalman out-of-sample backtest
#
# Outputs go to results/r/, in the same format as the Python outputs in results/,
# so that compare_r_python.py can compare them.

source("R/dnsmodel.R")
args <- commandArgs(trailingOnly = TRUE)
fast <- "--fast" %in% args
t_start <- Sys.time()

LAMBDA <- 0.0609
LAMBDA2 <- 2 * LAMBDA
HORIZONS <- c(1, 6, 12)
OUT <- "results/r"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

Y <- load_yields("data/bonds.csv")
Zns <- nelson_siegel_loadings(LAMBDA, MATURITIES)
Zsv <- svensson_loadings(LAMBDA, LAMBDA2, MATURITIES)

# ---- In-sample ------------------------------------------------------------------
fns <- fit_factors(Y, Zns)
fsv <- fit_factors(Y, Zsv)
t_k <- Sys.time()
kfit <- kalman_fit(Y)
kalman_seconds <- as.numeric(difftime(Sys.time(), t_k, units = "secs"))
kp <- kalman_params(kfit)
kstates <- kalman_filtered(kfit)
kfitted <- kstates %*% t(nelson_siegel_loadings(kp$lam, MATURITIES))

fit_tbl <- data.frame(maturity_months = MATURITIES,
                      `Nelson-Siegel (two-step)` = rmse(Y - fns$fitted) * 100,
                      `Svensson (two-step)` = rmse(Y - fsv$fitted) * 100,
                      `Nelson-Siegel (Kalman)` = rmse(Y - kfitted) * 100, check.names = FALSE)
write.csv(fit_tbl, file.path(OUT, "in_sample_rmse_bp.csv"), row.names = FALSE)
write.csv(data.frame(date = rownames(Y), level = fns$betas[, 1], slope = fns$betas[, 2], curvature = fns$betas[, 3]),
          file.path(OUT, "two_step_factors.csv"), row.names = FALSE)
write.csv(data.frame(date = rownames(Y), level = kstates[, 1], slope = kstates[, 2], curvature = kstates[, 3]),
          file.path(OUT, "kalman_factors.csv"), row.names = FALSE)

num <- function(x) paste(sprintf("%.10g", x), collapse = ", ")
writeLines(c("{",
  sprintf('  "lambda": %.10g,', kp$lam),
  sprintf('  "loglik": %.10g,', kfit$loglik),
  sprintf('  "mu": [%s],', num(kp$mu)),
  sprintf('  "A": [[%s], [%s], [%s]],', num(kp$A[1, ]), num(kp$A[2, ]), num(kp$A[3, ])),
  sprintf('  "measurement_sd_bp": [%s],', num(sqrt(kp$h) * 100)),
  sprintf('  "optim_convergence_code": %d,', kfit$convergence),
  sprintf('  "fit_seconds": %.1f', kalman_seconds),
  "}"), file.path(OUT, "kalman_params.json"))

# ---- Out-of-sample backtest -------------------------------------------------------
first <- which(rownames(Y) == "1994-01-31")
models <- list(
  "Random walk" = random_walk,
  "AR(1) on each yield" = ar1_yields,
  "DNS two-step (AR1)" = two_step_forecast(Zns, "ar1"),
  "DNS two-step (VAR1)" = two_step_forecast(Zns, "var1"),
  "Svensson two-step (AR1)" = two_step_forecast(Zsv, "ar1"))
if (!fast) models[["DNS Kalman (VAR1)"]] <- "kalman"
t_b <- Sys.time()
errors <- expanding_window_errors(Y, models, HORIZONS, first)
backtest_seconds <- as.numeric(difftime(Sys.time(), t_b, units = "secs"))

rows <- list(); dm_rows <- list()
for (name in names(errors)) for (h in HORIZONS) {
  e <- errors[[name]][[as.character(h)]]
  rows[[length(rows) + 1]] <- data.frame(model = name, horizon = h, t(rmse(e) * 100), check.names = FALSE)
  if (name != "Random walk") for (j in seq_along(MATURITIES)) {
    d <- diebold_mariano(e[, j], errors[["Random walk"]][[as.character(h)]][, j], h)
    dm_rows[[length(dm_rows) + 1]] <- data.frame(model = name, horizon = h, maturity_months = MATURITIES[j],
                                                 dm_stat = d[["stat"]], p_value = d[["p"]])
  }
}
rm_tbl <- do.call(rbind, rows)
names(rm_tbl)[-(1:2)] <- MATURITIES
write.csv(rm_tbl, file.path(OUT, "oos_rmse_bp.csv"), row.names = FALSE)
write.csv(do.call(rbind, dm_rows), file.path(OUT, "dm_tests_vs_random_walk.csv"), row.names = FALSE)

total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
writeLines(sprintf('{"total_seconds": %.1f, "kalman_fit_seconds": %.1f, "backtest_seconds": %.1f, "fast": %s, "r_version": "%s"}',
                   total, kalman_seconds, backtest_seconds, tolower(fast), R.version.string),
           file.path(OUT, "runtime.json"))
cat(sprintf("Kalman lambda = %.4f, loglik = %.2f, convergence code %d\n", kp$lam, kfit$loglik, kfit$convergence))
cat(sprintf("Total %.0f s\n", total))
